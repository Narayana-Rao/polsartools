import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d
from polsartools.utils.convert_matrices import C3_T3_mat

@time_it
def mf4cf(infolder, outname=None, window_size=1,write_flag=True,max_workers=None):

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
        output_filepaths.append(os.path.join(infolder, "Ps_mf4cf.tif"))
        output_filepaths.append(os.path.join(infolder, "Pd_mf4cf.tif"))
        output_filepaths.append(os.path.join(infolder, "Pv_mf4cf.tif"))
        output_filepaths.append(os.path.join(infolder, "Pc_mf4cf.tif"))
        output_filepaths.append(os.path.join(infolder, "Theta_FP_mf4cf.tif"))
        output_filepaths.append(os.path.join(infolder, "Tau_FP_mf4cf.tif"))
    
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
            processing_func=process_chunk_mf4cf,
            block_size=(512, 512), max_workers=max_workers, 
            num_outputs=6)

def process_chunk_mf4cf(chunks, window_size, input_filepaths):

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

    # det_T3 = t11s*(t22s*t33s-t23s*t32s)-t12s*(t21s*t33s-t23s*t31s)+t13s*(t21s*t32s-t22s*t31s)
    # trace_T3 = t11s + t22s + t33s
    # m1 = np.real(np.sqrt(1-(27*(det_T3/(trace_T3**3)))))

    reshaped_arr = T_T1.reshape(3, 3, -1).transpose(2, 0, 1)
    det_T3 = np.linalg.det(reshaped_arr)
    # del reshaped_arr
    det_T3 = det_T3.reshape(T_T1.shape[2], T_T1.shape[3])

    s0_f = T_T1[0,0,:,:] + T_T1[1,1,:,:] + T_T1[2,2,:,:] #trace_T3
    dop_f = np.real(np.sqrt(1-(27*(det_T3/(s0_f**3)))))

    k11_f = (T_T1[0,0,:,:] + T_T1[1,1,:,:] + T_T1[2,2,:,:])/2
    k44_f = (-T_T1[0,0,:,:] + T_T1[1,1,:,:] + T_T1[2,2,:,:])/2
    k14_f = np.imag(T_T1[1,2,:,:])

    
    # s0_f = trace_T3
    # dop_f = m1

    val1 = (4*dop_f*k11_f*k44_f)/(k44_f**2 - (1 + 4*dop_f**2)*k11_f**2)
    val2 = np.abs(k14_f)/(k11_f)
    


    theta_f = np.real(np.arctan(val1)) # separation for surface and dbl
    tau_f = np.real(np.arctan(val2)) # separation for helix
    # thet = np.rad2deg(thet)
    theta_FP = np.rad2deg(theta_f).astype(np.float32)
    tau_FP = np.rad2deg(tau_f).astype(np.float32)

    pc_f = (dop_f*s0_f*(np.sin(2*tau_f))).astype(np.float32)
    pv_f = ((1-dop_f)*s0_f).astype(np.float32)
    res_pow = s0_f - (pc_f + pv_f)
    ps_f = ((res_pow/2)*(1+np.sin((2*theta_f)))).astype(np.float32)
    pd_f = ((res_pow/2)*(1-np.sin((2*theta_f)))).astype(np.float32)

    return ps_f, pd_f, pv_f,pc_f,theta_FP,tau_FP