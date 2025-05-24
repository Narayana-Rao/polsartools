import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it

@time_it
def mf3cd(infolder, outname=None,  chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    if os.path.isfile(os.path.join(infolder,"T11.bin")):
        input_filepaths = [
            os.path.join(infolder, "T11.bin"), 
            os.path.join(infolder, "T12_real.bin"),
            os.path.join(infolder, "T12_imag.bin"),
            os.path.join(infolder, "T22.bin")
        ]
    else:
        raise(f"Invalid T2 folder!!")
        
    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "Ps_mf3cd.tif"))
        output_filepaths.append(os.path.join(infolder, "Pd_mf3cd.tif"))
        output_filepaths.append(os.path.join(infolder, "Pv_mf3cd.tif"))
        output_filepaths.append(os.path.join(infolder, "Theta_DP_mf3cd.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,processing_func=process_chunk_mf3cd,block_size=(512, 512), max_workers=max_workers,  num_outputs=4)

def process_chunk_mf3cd(chunks, window_size,*args):
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    
    t11_T1 = np.array(chunks[0])
    t12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    t21_T1 = np.conj(t12_T1)
    t22_T1 = np.array(chunks[3])


    if window_size>1:
        t11s = conv2d(np.real(t11_T1),kernel)+1j*conv2d(np.imag(t11_T1),kernel)
        t12s = conv2d(np.real(t12_T1),kernel)+1j*conv2d(np.imag(t12_T1),kernel)
        t21s = np.conj(t12_T1)
        t22s = conv2d(np.real(t22_T1),kernel)+1j*conv2d(np.imag(t22_T1),kernel)
    
    else:
        t11s = t11_T1
        t12s = t12_T1
        t21s = t21_T1
        t22s = t22_T1
    
    det_T2 = t11s*t22s-t12s*t21s
    trace_T2 = t11s + t22s

    m1 = np.real(np.sqrt(1-(4*(det_T2/(trace_T2**2)))))
    h = (t11s - t22s)
    g = t22s
    span = t11s + t22s

    
    val = (m1*span*h)/(t11s*g+m1**2*span**2);
    thet = np.real(np.arctan(val))
    # thet = np.rad2deg(thet)
    theta_DP = np.rad2deg(thet)

    
    Ps_DP = (((m1*(span)*(1+np.sin(2*thet))/2)))
    Pd_DP = (((m1*(span)*(1-np.sin(2*thet))/2)))
    Pv_DP = (span*(1-m1))


    return np.real(Ps_DP),np.real(Pd_DP),np.real(Pv_DP),np.real(theta_DP)