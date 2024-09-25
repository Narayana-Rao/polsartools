import os
import numpy as np
# from ..utils.utils import process_chunks_parallel,conv2d,eig22
from polsartools.utils.utils import process_chunks_parallel,conv2d,eig22

def rvidp(infolder, outname=None, window_size=1,write_flag=True,max_workers=None):
    input_filepaths = [
        os.path.join(infolder, "C11.bin"), 
        os.path.join(infolder, "C12_real.bin"),
        os.path.join(infolder, "C12_imag.bin"),
        os.path.join(infolder, "C22.bin")
    ]
    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "rvidp.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,processing_func=process_chunk_rvidp,block_size=(512, 512), max_workers=max_workers,  num_outputs=1)

def process_chunk_rvidp(chunks, window_size):
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
        rvidp = np.real(4*c22s/c2_trace)
    else:
        c2_det = (c11_T1*c22_T1-c12_T1*c21_T1)
        c2_trace = c11_T1+c22_T1
        rvidp = np.real(4*c22_T1/c2_trace)

    return rvidp