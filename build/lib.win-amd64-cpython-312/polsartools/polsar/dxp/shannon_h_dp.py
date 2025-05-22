import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d, eig22


"""
normlized shannon entropy parameters are not agreeing with polsarpro
others are fine

"""
@time_it
def shannon_h_dp(infolder, outname=None,  chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    input_filepaths = [
        os.path.join(infolder, "C11.bin"), 
        os.path.join(infolder, "C12_real.bin"),
        os.path.join(infolder, "C12_imag.bin"),
        os.path.join(infolder, "C22.bin")
    ]
    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "H_Shannon.tif"))
        output_filepaths.append(os.path.join(infolder, "HI_Shannon.tif"))
        output_filepaths.append(os.path.join(infolder, "HP_Shannon.tif"))
        
        # output_filepaths.append(os.path.join(infolder, "HS_norm.tif"))
        # output_filepaths.append(os.path.join(infolder, "HSI_norm.tif"))
        # output_filepaths.append(os.path.join(infolder, "HSP_norm.tif"))
        

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,processing_func=process_chunk_shannondp,block_size=(512, 512), max_workers=max_workers,  num_outputs=3)

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
  
    eps  = 1e-8
    D = evals[:,0]*evals[:,1]
    I = evals[:,0]+evals[:,1]
    
    # Barakat degree of polarization
    DoP = np.ones(rows*cols).astype(np.float32) - 4* D / (I*I + eps)

    HSP = np.zeros(rows*cols).astype(np.float32)
    # HSI = np.zeros(rows*cols).astype(np.float32)
    # HS = np.zeros(rows*cols).astype(np.float32)

    condition = (np.ones(rows*cols) - DoP) < eps
    HSP = np.where(condition, 0, np.log(np.abs(np.ones(rows*cols) - DoP)))
    HSP[np.isinf(HSP)] = np.nan
    HSP[HSP==0] = np.nan
    
    with np.errstate(divide='ignore', invalid='ignore'):
        HSI = 2 * np.log(np.exp(1) * np.pi * I / 2)
        HSI[np.isinf(HSI)] = np.nan
    
    HS = np.nansum(np.dstack((HSP, HSI)), 2)

    """ Normalization will not not work as expected if we are processing individual blocks of data. 
    Therefore we will normalize the whole image at the end.
    """
    # HSP_norm = (HSP - np.nanmin(HSP)) / (np.nanmax(HSP) - np.nanmin(HSP))
    # HSI_norm = (HSI - np.nanmin(HSI)) / (np.nanmax(HSI) - np.nanmin(HSI))
    # HS_norm = (HS - np.nanmin(HS)) / (np.nanmax(HS) - np.nanmin(HS))

    

    return np.real(HS).reshape(rows,cols),np.real(HSI).reshape(rows,cols),np.real(HSP).reshape(rows,cols) #np.real(HS_norm).reshape(rows,cols),np.real(HSI_norm).reshape(rows,cols),np.real(HSP_norm).reshape(rows,cols) 