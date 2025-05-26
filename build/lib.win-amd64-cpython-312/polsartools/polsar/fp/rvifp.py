import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it
from polsartools.utils.convert_matrices import C3_T3_mat
from .fp_infiles import fp_c3t3files

@time_it
def rvifp(infolder, outname=None, window_size=1,write_flag=True,max_workers=None):

    input_filepaths = fp_c3t3files(infolder)
    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "rvifp.tif"))
    
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
            processing_func=process_chunk_rvifp,block_size=(512, 512), max_workers=max_workers,  num_outputs=1)


def process_chunk_rvifp(chunks, window_size,input_filepaths,*args):

    t11_T1 = np.array(chunks[0])
    t12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    t13_T1 = np.array(chunks[3])+1j*np.array(chunks[4])
    t21_T1 = np.conj(t12_T1)
    t22_T1 = np.array(chunks[5])
    t23_T1 = np.array(chunks[6])+1j*np.array(chunks[7])
    t31_T1 = np.conj(t13_T1)
    t32_T1 = np.conj(t23_T1)
    t33_T1 = np.array(chunks[8])

    T_T1 = np.array([[t11_T1, t12_T1, t13_T1], 
                     [t21_T1, t22_T1, t23_T1], 
                     [t31_T1, t32_T1, t33_T1]])

    if window_size>1:
        kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

        t11f = conv2d(T_T1[0,0,:,:],kernel)
        t12f = conv2d(np.real(T_T1[0,1,:,:]),kernel)+1j*conv2d(np.imag(T_T1[0,1,:,:]),kernel)
        t13f = conv2d(np.real(T_T1[0,2,:,:]),kernel)+1j*conv2d(np.imag(T_T1[0,2,:,:]),kernel)
        
        t21f = np.conj(t12f) 
        t22f = conv2d(T_T1[1,1,:,:],kernel)
        t23f = conv2d(np.real(T_T1[1,2,:,:]),kernel)+1j*conv2d(np.imag(T_T1[1,2,:,:]),kernel)

        t31f = np.conj(t13f) 
        t32f = np.conj(t23f) 
        t33f = conv2d(T_T1[2,2,:,:],kernel)

        T_T1 = np.array([[t11f, t12f, t13f], [t21f, t22f, t23f], [t31f, t32f, t33f]])


    reshaped_arr = T_T1.reshape(3, 3, -1).transpose(2, 0, 1)
    eigenvalues = np.linalg.eigvals(reshaped_arr)
    sorted_eigenvalues = np.sort(eigenvalues, axis=1)[:, ::-1]
    sorted_eigenvalues = sorted_eigenvalues.reshape(T_T1.shape[2], T_T1.shape[3], 3)
    
    p1 = sorted_eigenvalues[:,:,0]/(sorted_eigenvalues[:,:,0] + sorted_eigenvalues[:,:,1] + sorted_eigenvalues[:,:,2])
    p2 = sorted_eigenvalues[:,:,1]/(sorted_eigenvalues[:,:,0] + sorted_eigenvalues[:,:,1] + sorted_eigenvalues[:,:,2])
    p3 = sorted_eigenvalues[:,:,2]/(sorted_eigenvalues[:,:,0] + sorted_eigenvalues[:,:,1] + sorted_eigenvalues[:,:,2])
    
    p1[p1<0] = 0
    p2[p2<0] = 0
    p3[p3<0] = 0
        
    p1[p1>1] = 1
    p2[p2>1] = 1
    p3[p3>1] = 1

    
    rvi = np.real((4*p3)/(p1 + p2 + p3))

    idx = np.argwhere(rvi>1)

    rvi[idx] = (3/4)*rvi[idx]
    rvi[~idx] = rvi[~idx]
    rvi[rvi==0] = np.nan

    rvi = np.real(rvi)

    return rvi
