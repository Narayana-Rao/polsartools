import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d
from polsartools.utils.convert_matrices import C3_T3_mat

@time_it
def neufp(infolder, outname=None, chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):

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
        output_filepaths.append(os.path.join(infolder, "Neumann_psi.tif"))
        output_filepaths.append(os.path.join(infolder, "Neumann_delta_mod.tif"))
        output_filepaths.append(os.path.join(infolder, "Neumann_delta_pha.tif"))
        output_filepaths.append(os.path.join(infolder, "Neumann_tau.tif"))
    
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
            processing_func=process_chunk_neufp,
            block_size=(512, 512), max_workers=max_workers, 
            num_outputs=4)

def process_chunk_neufp(chunks, window_size, input_filepaths, *args):

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


    if 'C11' in input_filepaths[0] and 'C22' in input_filepaths[5] and 'C33' in input_filepaths[8]:
        C11 = np.array(chunks[0])
        C12 = np.array(chunks[1])+1j*np.array(chunks[2])
        C13 = np.array(chunks[3])+1j*np.array(chunks[4])
        C21 = np.conj(C12)
        C22 = np.array(chunks[5])
        C23 = np.array(chunks[6])+1j*np.array(chunks[7])
        C31 = np.conj(C13)
        C32 = np.conj(C23)
        C33 = np.array(chunks[8])
        C3 = np.array([[C11, C12, C13], 
                         [C21, C22, C23], 
                         [C31, C32, C33]])

        T_T1 = C3_T3_mat(C3)


    if window_size>1:
        kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

        t11f = conv2d(T_T1[0,0,:,:],kernel)
        t12f = conv2d(T_T1[0,1,:,:],kernel)
        t13f = conv2d(T_T1[0,2,:,:],kernel)
        
        t21f = conv2d(T_T1[1,0,:,:],kernel)
        t22f = conv2d(T_T1[1,1,:,:],kernel)
        t23f = conv2d(T_T1[1,2,:,:],kernel)

        t31f = conv2d(T_T1[2,0,:,:],kernel)
        t32f = conv2d(T_T1[2,1,:,:],kernel)
        t33f = conv2d(T_T1[2,2,:,:],kernel)

        T_T1 = np.array([[t11f, t12f, t13f], [t21f, t22f, t23f], [t31f, t32f, t33f]])

    

    _,_,rows,cols = np.shape(T_T1)
    T_T1 = T_T1.reshape(9, rows, cols)
    
    # Indices for vectorized access
    i, j = np.indices((rows, cols))

    # Extract components of T3_stack once
    T3_stack_real = np.real(T_T1[:,i, j ])
    T3_stack_imag = np.imag(T_T1[:, i, j])

    T_00_0 = T3_stack_real[0, :, :]
    T_01_0 = T3_stack_real[1, :, :]
    T_01_1 = T3_stack_imag[1, :, :]
    T_02_0 = T3_stack_real[2, :, :]
    T_02_1 = T3_stack_imag[2, :, :]
    T_11_0 = T3_stack_real[4, :, :]
    T_12_0 = T3_stack_real[5, :, :]
    T_22_0 = T3_stack_real[8, :, :]

    # Compute the Phi value using the given formula
    Phi = 0.25 * (np.pi + np.arctan2(-2. * T_12_0, T_22_0 - T_11_0))

    # Adjust Phi for the correct orientation
    Phi[Phi <= np.pi / 4.] = Phi[Phi <= np.pi / 4.]
    Phi[Phi > np.pi / 4.] -= np.pi / 2.

    # Convert Phi to degrees
    Neumann_psi = Phi * 180. / np.pi

    # Pre-compute trigonometric terms
    cos_2Phi = np.cos(2 * Phi)
    sin_2Phi = np.sin(2 * Phi)
    sin_4Phi = np.sin(4 * Phi)

    # Coherency Matrix de-orientation
    T110 = T_00_0
    T12re0 = T_01_0 * cos_2Phi - T_02_0 * sin_2Phi
    T12im0 = T_01_1 * cos_2Phi - T_02_1 * sin_2Phi
    T220 = T_11_0 * cos_2Phi**2 - T_12_0 * sin_4Phi + T_22_0 * sin_2Phi**2
    T330 = T_11_0 * sin_2Phi**2 - T_12_0 * sin_4Phi + T_22_0 * cos_2Phi**2

    # Compute Neumann_delta_mod, Neumann_delta_pha, and Neumann_tau
    Neumann_delta_mod = np.sqrt((T220 + T330) / (T110 + np.finfo(float).eps))
    Neumann_delta_pha = np.arctan2(T12im0, T12re0) * 180. / np.pi
    Neumann_tau = 1. - ((np.sqrt(T12re0 ** 2 + T12im0 ** 2) / T110) / (Neumann_delta_mod + np.finfo(float).eps))

    return Neumann_psi, Neumann_delta_mod, Neumann_delta_pha,Neumann_tau