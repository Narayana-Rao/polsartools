import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d, eig22

@time_it
def prvidp(infolder, outname=None, chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    input_filepaths = [
        os.path.join(infolder, "C11.bin"), 
        os.path.join(infolder, "C12_real.bin"),
        os.path.join(infolder, "C12_imag.bin"),
        os.path.join(infolder, "C22.bin")
    ]
    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "prvidp.tif"))
    
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,processing_func=process_chunk_prvidp,block_size=(512, 512), max_workers=max_workers, num_outputs=1)

def process_chunk_prvidp(chunks, window_size,*args):
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    if window_size>1:
        c11_T1r = conv2d(np.real(c11_T1),kernel)
        c11_T1i = conv2d(np.imag(c11_T1),kernel)
        c11s = c11_T1r+1j*c11_T1i

        c12_T1r = conv2d(np.real(c12_T1),kernel)
        c12_T1i = conv2d(np.imag(c12_T1),kernel)
        c12s = c12_T1r+1j*c12_T1i


        c21_T1r = conv2d(np.real(c21_T1),kernel)
        c21_T1i = conv2d(np.imag(c21_T1),kernel)
        c21s = c21_T1r+1j*c21_T1i


        c22_T1r = conv2d(np.real(c22_T1),kernel)
        c22_T1i = conv2d(np.imag(c22_T1),kernel)
        c22s = c22_T1r+1j*c22_T1i

        
        c2_det = (c11s*c22s-c12s*c21s)
        c2_trace = c11s+c22s
        dopdp = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))
        prvidp = np.real((1-dopdp)*c22s)
    else:
        c2_det = (c11_T1*c22_T1-c12_T1*c21_T1)
        c2_trace = c11_T1+c22_T1
        dopdp = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))
        prvidp = np.real((1-dopdp)*c22_T1)  

    return prvidp