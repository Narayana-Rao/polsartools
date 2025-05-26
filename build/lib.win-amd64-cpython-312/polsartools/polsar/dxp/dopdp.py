import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it
from .dxp_infiles import dxpc2files
@time_it
def dopdp(infolder, outname=None, window_size=1,write_flag=True,max_workers=None):
    
    input_filepaths = dxpc2files(infolder)
    
    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "dopdp.tif"))
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, 
            write_flag=write_flag,processing_func=process_chunk_dopdp,block_size=(512, 512), 
            max_workers=max_workers,  num_outputs=1)

def process_chunk_dopdp(chunks, window_size,*args):
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    if window_size>1:
        c11_T1 = conv2d(np.real(c11_T1),kernel)+1j*conv2d(np.imag(c11_T1),kernel)
        c12_T1 = conv2d(np.real(c12_T1),kernel)+1j*conv2d(np.imag(c12_T1),kernel)
        c21_T1 = conv2d(np.real(c21_T1),kernel)+1j*conv2d(np.imag(c21_T1),kernel)
        c22_T1 = conv2d(np.real(c22_T1),kernel)+1j*conv2d(np.imag(c22_T1),kernel)
        
    c2_det = (c11_T1*c22_T1-c12_T1*c21_T1)
    c2_trace = c11_T1+c22_T1
    dopdp = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))

    return dopdp