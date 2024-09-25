import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel

def rvifp(infolder, outname=None, window_size=1,write_flag=True,max_workers=None):

    if os.path.isfile(os.path.join(infolder,"T11.bin")):
        input_filepaths = [
        os.path.join(infolder,"T11.bin"),
        
        os.path.join(infolder,'T12_real.bin'), os.path.join(infolder,'T12_imag.bin'),  
        os.path.join(infolder,'T13_real.bin'), os.path.join(infolder,'T13_imag.bin'),
        os.path.join(infolder,"T22.bin"),
        os.path.join(infolder,'T23_real.bin'), os.path.join(infolder,'T23_imag.bin'),  
     
        os.path.join(infolder,"T33.bin"),
        ]
    elif os.path.isfile(os.path.join(infolder,"C11.bin")):

        input_filepaths = [
        os.path.join(infolder,"C11.bin"),
        
        os.path.join(infolder,'C12_real.bin'), os.path.join(infolder,'C12_imag.bin'),  
        os.path.join(infolder,'C13_real.bin'), os.path.join(infolder,'C13_imag.bin'),
        os.path.join(infolder,"C22.bin"),
        os.path.join(infolder,'C23_real.bin'), os.path.join(infolder,'C23_imag.bin'),  
     
        os.path.join(infolder,"C33.bin"),
        ]



    else:

    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "rvifp.tif"))
    
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
            processing_func=process_chunk_rvifp,block_size=(512, 512), max_workers=max_workers,  num_outputs=1)


def process_chunk_rvifp(chunks, window_size):

    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

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

    # if window_size>1:
    #     c11_T1r = conv2d(np.real(c11_T1),kernel)
    #     c11_T1i = conv2d(np.imag(c11_T1),kernel)
    #     c11s = c11_T1r+1j*c11_T1i

    #     c12_T1r = conv2d(np.real(c12_T1),kernel)
    #     c12_T1i = conv2d(np.imag(c12_T1),kernel)
    #     c12s = c12_T1r+1j*c12_T1i


    #     c21_T1r = conv2d(np.real(c21_T1),kernel)
    #     c21_T1i = conv2d(np.imag(c21_T1),kernel)
    #     c21s = c21_T1r+1j*c21_T1i


    #     c22_T1r = conv2d(np.real(c22_T1),kernel)
    #     c22_T1i = conv2d(np.imag(c22_T1),kernel)
    #     c22s = c22_T1r+1j*c22_T1i
