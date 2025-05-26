import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it
from .cp_infiles import cpc2files
@time_it
def dopcp(infolder, outname=None, chi_in=45, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    input_filepaths = cpc2files(infolder)

    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "dopcp.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, 
                write_flag=write_flag,processing_func=process_chunk_dopcp,block_size=(512, 512), 
                max_workers=max_workers,  num_outputs=1, chi_in=chi_in,psi_in=psi_in)

def process_chunk_dopcp(chunks, window_size,chi_in,psi_in):

    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    ncols,nrows = np.shape(c11_T1)

    if window_size>1:
        c11_T1 = conv2d(np.real(c11_T1),kernel)+1j*conv2d(np.imag(c11_T1),kernel)
        c12_T1 = conv2d(np.real(c12_T1),kernel)+1j*conv2d(np.imag(c12_T1),kernel)
        c21_T1 = conv2d(np.real(c21_T1),kernel)+1j*conv2d(np.imag(c21_T1),kernel)
        c22_T1 = conv2d(np.real(c22_T1),kernel)+1j*conv2d(np.imag(c22_T1),kernel)
        
    # Compute Stokes parameters
    s0 = c11_T1 + c22_T1
    s1 = c11_T1 - c22_T1
    s2 = c12_T1 + c21_T1
    s3 = np.where(chi_in >= 0, 1j * (c12_T1 - c21_T1), -1j * (c12_T1 - c21_T1))

    dop= np.sqrt(np.power(s1,2) + np.power(s2,2) + np.power(s3,2))/(s0);   

    return dop