import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d, eig22

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
    print(np.shape(C2_stack),np.shape(data))
    rows, cols,_ = C2_stack.shape
    
    
    evals, evecs = np.linalg.eig(data)
    
    
    evals[:,0][evals[:,0] <0] = 0
    evals[:,1][evals[:,1] >1] = 1
    
    eval_norm1 = (evals[:,1])/(evals[:,0] + evals[:,1])
    eval_norm1[eval_norm1<0]=0
    eval_norm1[eval_norm1>1]=1
    
    
    eval_norm2 = (evals[:,0])/(evals[:,0] + evals[:,1])
    eval_norm2[eval_norm2<0]=0
    eval_norm2[eval_norm2>1]=1
    
    # # %Alpha 1
    eig_vec_r1 = np.real(evecs[:,0,1])
    eig_vec_c1 = np.imag(evecs[:,0,1])
    # alpha1 = np.arccos(np.sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1))*180/np.pi
    
    # # %Alpha 2
    eig_vec_r2 = np.real(evecs[:,0,0])
    eig_vec_c2 = np.imag(evecs[:,0,0])
    # alpha2 = np.arccos(np.sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2))*180/np.pi
    print('here...1 ')
    eps  = 1e-8
    D = evals[:,0]*evals[:,1]
    I = evals[:,0]+evals[:,1]
    print('here...2 ',np.shape(I),np.shape(D))
    DegPol = np.ones(rows*cols).astype(np.float32) - 4* D / (I*I + eps)
    print('here...2.5 ')
    HSP = np.zeros(rows*cols).astype(np.float32)
    HSI = np.zeros(rows*cols).astype(np.float32)
    HS = np.zeros(rows*cols).astype(np.float32)
    HSPN = np.zeros(rows*cols).astype(np.float32)
    HSIN = np.zeros(rows*cols).astype(np.float32)
    HSN = np.zeros(rows*cols).astype(np.float32)
    
    print('here...3')
    
    # if ((1. - DegPol) < eps) M_out[Flag[HSP]][lig][col] = 0.;
    # else M_out[Flag[HSP]][lig][col] = log(fabs(1. - DegPol));

    condition = (np.ones(rows*cols) - DegPol) < eps
    HSP = np.where(condition, 0, np.log(np.abs(np.ones(rows*cols) - DegPol)))
    
    HSI= 2 * np.log(np.exp(1)*np.pi*I/2)
    HS = HSP + HSI
    
    # # # %Cloude Alpha
    # alpha_ = (eval_norm1*alpha1 + eval_norm2*alpha2)
    # alpha_ = alpha_.reshape(C2_stack.shape[0],C2_stack.shape[1])
    # # # %Entropy
    # H = - eval_norm1*np.log10(eval_norm1)/np.log10(2) - eval_norm2*np.log10(eval_norm2)/np.log10(2)
    # H = H.reshape(C2_stack.shape[0],C2_stack.shape[1])

    # # alpha1 = alpha1.reshape(C2_stack.shape[0],C2_stack.shape[1])
    # # alpha2 = alpha2.reshape(C2_stack.shape[0],C2_stack.shape[1])

    # # print(np.nanmean(H),np.nanmean(alpha_))

    return np.real(HS).reshape(rows,cols),np.real(HSI).reshape(rows,cols),np.real(HSP).reshape(rows,cols)