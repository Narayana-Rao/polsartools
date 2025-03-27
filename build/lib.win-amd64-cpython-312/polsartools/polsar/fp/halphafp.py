import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d
from polsartools.utils.convert_matrices import T3_C3_mat

@time_it
def halphafp(infolder, outname=None, chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):

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
        print(f"Invalid C3 or T3 folder!!")

    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "H_fp.tif"))
        output_filepaths.append(os.path.join(infolder, "alpha_fp.tif"))
        output_filepaths.append(os.path.join(infolder, "anisotropy_fp.tif"))
        output_filepaths.append(os.path.join(infolder, "e1_norm.tif"))
        output_filepaths.append(os.path.join(infolder, "e2_norm.tif"))
        output_filepaths.append(os.path.join(infolder, "e3_norm.tif"))
    
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
            processing_func=process_chunk_halphafp,
            block_size=(512, 512), max_workers=max_workers, 
            num_outputs=6)

def process_chunk_halphafp(chunks, window_size, input_filepaths, *args):

    # additional_arg1 = args[0] if len(args) > 0 else None
    # additional_arg2 = args[1] if len(args) > 1 else None

    if 'T11' in input_filepaths[0] and 'T22' in input_filepaths[5] and 'T33' in input_filepaths[8]:
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
        # T_T1 = T3_C3_mat(T3)


    elif 'C11' in input_filepaths[0] and 'C22' in input_filepaths[5] and 'C33' in input_filepaths[8]:
        C11 = np.array(chunks[0])
        C12 = np.array(chunks[1])+1j*np.array(chunks[2])
        C13 = np.array(chunks[3])+1j*np.array(chunks[4])
        C21 = np.conj(C12)
        C22 = np.array(chunks[5])
        C23 = np.array(chunks[6])+1j*np.array(chunks[7])
        C31 = np.conj(C13)
        C32 = np.conj(C23)
        C33 = np.array(chunks[8])
        T_T1 = np.array([[C11, C12, C13], 
                         [C21, C22, C23], 
                         [C31, C32, C33]])

    else:
        raise ValueError("Invalid input matrices. Ensure the input is either T3 or C3 matrix foolder.")


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

    _,_,rows,cols = np.shape(T_T1)
    T_T1 = T_T1.reshape(9, rows, cols)
    
    # Indices for vectorized access
    i, j = np.indices((rows, cols))
    
    T_T1 = np.dstack((T_T1[0,:,:],T_T1[1,:,:],T_T1[2,:,:],
                      T_T1[3,:,:],T_T1[4,:,:],T_T1[5,:,:],T_T1[6,:,:],T_T1[7,:,:],T_T1[8,:,:]))
    
    data = T_T1.reshape( T_T1.shape[0]*T_T1.shape[1], T_T1.shape[2]).reshape((-1,3,3))

    evals_, evecs_ = np.linalg.eig(data.reshape(-1, 3, 3))

    # Sort eigenvalues for each pixel in descending order; 
    sorted_indices = np.argsort(evals_, axis=-1)[:, ::-1] 
    
    # Reorder eigenvalues and eigenvectors based on sorted indices
    evals = np.take_along_axis(evals_, sorted_indices, axis=-1)  # Reorder eigenvalues
    
    # To reorder the eigenvectors, we use  indexing
    # Use `sorted_indices` to index along the second axis (the 3 components of the eigenvector)
    evecs = np.array([evecs_[i, :, sorted_indices[i]] for i in range(evecs_.shape[0])])
    
    # print('Eigen!')
    # Not sure if this is required; if enabled entropy values are different between polsarpro and this code
    # evals[:,0][evals[:,0] <0] = 0
    # evals[:,1][evals[:,1] <0] = 0
    # evals[:,2][evals[:,2] <0] = 0
    # evals[:,0][evals[:,0] >1] = 1
    # evals[:,1][evals[:,1] >1] = 1
    # evals[:,2][evals[:,2] >1] = 1
    
    eval_norm1 = (evals[:,0])/(evals[:,0] + evals[:,1]+ evals[:,2])
    eval_norm1[eval_norm1<0]=0
    # eval_norm1[eval_norm1>1]=1
    
    
    eval_norm2 = (evals[:,1])/(evals[:,0] + evals[:,1]+ evals[:,2])
    eval_norm2[eval_norm2<0]=0
    # eval_norm2[eval_norm2>1]=1
    
    
    eval_norm3 = (evals[:,2])/(evals[:,0] + evals[:,1]+evals[:,2])
    eval_norm3[eval_norm3<0]=0
    # eval_norm3[eval_norm3>1]=1
      
    
    # # %Alpha 1
    eig_vec_r1 = np.real(evecs[:,0,0])
    eig_vec_c1 = np.imag(evecs[:,0,0])
    alpha1 = np.arccos(np.sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1))*180/np.pi
    
    # # %Alpha 2
    eig_vec_r2 = np.real(evecs[:,0,1])
    eig_vec_c2 = np.imag(evecs[:,0,1])
    alpha2 = np.arccos(np.sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2))*180/np.pi
    
    # # %Alpha 3
    eig_vec_r3 = np.real(evecs[:,0,2])
    eig_vec_c3 = np.imag(evecs[:,0,2])
    alpha3 = np.arccos(np.sqrt(eig_vec_r3*eig_vec_r3 + eig_vec_c3*eig_vec_c3))*180/np.pi
        
    # # %Cloude Alpha
    alpha_ = (eval_norm1*alpha1 + eval_norm2*alpha2+ eval_norm3*alpha3)
    alpha_ = alpha_.reshape(rows,cols)
    
    # print('Alpha!')
    
    # # %Entropy
    H = - eval_norm1*np.log10(eval_norm1)/np.log10(3) - eval_norm2*np.log10(eval_norm2)/np.log10(3) - eval_norm3*np.log10(eval_norm3)/np.log10(3)
    H = H.reshape(rows,cols)

    # alpha1 = alpha1.reshape(rows,cols)
    # alpha2 = alpha2.reshape(rows,cols)
    # alpha3 = alpha3.reshape(rows,cols)
    
    ## POLARIMETRIC SCATTERING ANISOTROPY (A)
    Anisotropy = (eval_norm2-eval_norm3)/(eval_norm2+eval_norm3)
    
    return H,alpha_,Anisotropy.reshape(rows,cols),eval_norm1.reshape(rows,cols),eval_norm2.reshape(rows,cols),eval_norm3.reshape(rows,cols)