import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d, eig22


"""
normlized shannon entropy parameters are not agreeing with polsarpro
others are fine

"""
@time_it
def shannondp(infolder, outname=None,  chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    input_filepaths = [
        os.path.join(infolder, "C11.bin"), 
        os.path.join(infolder, "C12_real.bin"),
        os.path.join(infolder, "C12_imag.bin"),
        os.path.join(infolder, "C22.bin")
    ]
    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "HS.tif"))
        output_filepaths.append(os.path.join(infolder, "HSI.tif"))
        output_filepaths.append(os.path.join(infolder, "HSP.tif"))
        
        output_filepaths.append(os.path.join(infolder, "HS_norm.tif"))
        output_filepaths.append(os.path.join(infolder, "HSI_norm.tif"))
        output_filepaths.append(os.path.join(infolder, "HSP_norm.tif"))
        

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,processing_func=process_chunk_shannondp,block_size=(512, 512), max_workers=max_workers,  num_outputs=6)

def process_chunk_shannondp(chunks, window_size,*args):
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    # C2_stack = np.zeros((np.shape(c11_T1)[0],np.shape(c11_T1)[1],4))
    C2_stack = np.dstack((c11_T1,c12_T1,np.conj(c12_T1),c22_T1)).astype(np.complex64)

    if window_size>1:
        C2_stack[:,:,0] = conv2d(np.real(c11_T1),kernel)+1j*conv2d(np.imag(c11_T1),kernel)
        C2_stack[:,:,1] = conv2d(np.real(c12_T1),kernel)+1j*conv2d(np.imag(c12_T1),kernel)
        C2_stack[:,:,2] = conv2d(np.real(c21_T1),kernel)+1j*conv2d(np.imag(c21_T1),kernel)
        C2_stack[:,:,3] = conv2d(np.real(c22_T1),kernel)+1j*conv2d(np.imag(c22_T1),kernel)

    data = C2_stack.reshape( C2_stack.shape[0]*C2_stack.shape[1], C2_stack.shape[2] ).reshape((-1,2,2))
    rows, cols,_ = C2_stack.shape
    
    
    evals, evecs = np.linalg.eig(data)
    
    
    evals[:,0][evals[:,0] <0] = 0
    evals[:,1][evals[:,1] >1] = 1
    
    eval_norm1 = (evals[:,0])/(evals[:,0] + evals[:,1])
    # eval_norm1[eval_norm1<0]=0
    # eval_norm1[eval_norm1>1]=1
    # eval_norm1 = (evals[:,0])/(evals[:,0] + evals[:,1])
    
    eval_norm2 = (evals[:,1])/(evals[:,0] + evals[:,1])
    
    # eval_norm2[eval_norm2<0]=0
    # eval_norm2[eval_norm2>1]=1
    
    # # %Alpha 1
    eig_vec_r1 = np.real(evecs[:,0,0])
    eig_vec_c1 = np.imag(evecs[:,0,0])
    # alpha1 = np.arccos(np.sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1))*180/np.pi
    
    # # %Alpha 2
    eig_vec_r2 = np.real(evecs[:,0,1])
    eig_vec_c2 = np.imag(evecs[:,0,1])
    # alpha2 = np.arccos(np.sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2))*180/np.pi
    eps  = 1e-8
    D = evals[:,0]*evals[:,1]
    I = evals[:,0]+evals[:,1]
    DegPol = np.ones(rows*cols).astype(np.float32) - 4* D / (I*I + eps)

    HSP = np.zeros(rows*cols).astype(np.float32)
    HSI = np.zeros(rows*cols).astype(np.float32)
    HS = np.zeros(rows*cols).astype(np.float32)

    condition = (np.ones(rows*cols) - DegPol) < eps
    HSP = np.where(condition, 0, np.log(np.abs(np.ones(rows*cols) - DegPol)))
    
    HSI= 2 * np.log(np.exp(1)*np.pi*I/2)
    HS = HSP + HSI
    

    D_norm = eval_norm1*eval_norm2
    I_norm = eval_norm1+eval_norm2
    DegPol_norm = np.ones(rows*cols).astype(np.float32) - 4* D_norm / (I_norm*I_norm + eps)
    HSP_norm = np.zeros(rows*cols).astype(np.float32)
    HSI_norm = np.zeros(rows*cols).astype(np.float32)
    HS_norm = np.zeros(rows*cols).astype(np.float32)
    
    
    condition = (np.ones(rows*cols) - DegPol_norm) < eps
    HSP_norm = np.where(condition, 0, np.log(np.abs(np.ones(rows*cols) - DegPol_norm)))
    HSI_norm= 2 * np.log(np.exp(1)*np.pi*I_norm/2)
    
    
    # HSP_norm = (HSP_norm - np.nanmin(HSP_norm)) / (np.nanmax(HSP_norm) - np.nanmin(HSP_norm))
    
    HS_norm = HSP_norm + HSI_norm
    
    
    

    return np.real(HS).reshape(rows,cols),np.real(HSI).reshape(rows,cols),np.real(HSP).reshape(rows,cols),np.real(HS_norm).reshape(rows,cols),np.real(HSI_norm).reshape(rows,cols),np.real(HSP_norm).reshape(rows,cols) 