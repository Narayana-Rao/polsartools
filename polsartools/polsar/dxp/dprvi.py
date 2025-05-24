import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it

@time_it
def dprvi(infolder, outname=None, chi_in=0, psi_in=0,window_size=1,write_flag=True,max_workers=None):
    input_filepaths = [
        os.path.join(infolder, "C11.bin"), 
        os.path.join(infolder, "C12_real.bin"),
        os.path.join(infolder, "C12_imag.bin"),
        os.path.join(infolder, "C22.bin")
    ]
    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "dprvi.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, 
        write_flag=write_flag,processing_func=process_chunk_dprvi,block_size=(512, 512), max_workers=max_workers,  num_outputs=1)

def process_chunk_dprvi(chunks, window_size,*args):
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    if window_size>1:
        c11s = conv2d(np.real(c11_T1),kernel)+1j*conv2d(np.imag(c11_T1),kernel)
        c12s = conv2d(np.real(c12_T1),kernel)+1j*conv2d(np.imag(c12_T1),kernel)
        c21s = conv2d(np.real(c21_T1),kernel)+1j*conv2d(np.imag(c21_T1),kernel)
        c22s = conv2d(np.real(c22_T1),kernel)+1j*conv2d(np.imag(c22_T1),kernel)

        c2_det = (c11s*c22s-c12s*c21s)
        c2_trace = c11s+c22s
        # t2_span = t11s*t22s
        m = (np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))
        egv1,egv2 = eig22(np.dstack([c11s,c12s,c21s,c22s]))
        egf = np.vstack([egv1,egv2])
        egfmax = egf.max(axis=0)#.reshape(np.shape(C2_stack[:,:,0]))
        beta = (egfmax/(egv1+egv2)).reshape(np.shape(c11s))
        dprvi = np.real(1-(m*beta))
    else:
        c2_det = (c11_T1*c22_T1-c12_T1*c21_T1)
        c2_trace = c11_T1+c22_T1
        m = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))

        egv1,egv2 = eig22(np.dstack([c11_T1,c12_T1,c21_T1,c22_T1]))
        egf = np.vstack([egv1,egv2])
        egfmax = egf.max(axis=0)#.reshape(np.shape(C2_stack[:,:,0]))
        beta = (egfmax/(egv1+egv2)).reshape(np.shape(c11_T1))
        dprvi = np.real(1-(m*beta))


    return dprvi