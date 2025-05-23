import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it

@time_it
def mf3cc(infolder, outname=None, chi_in=45, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    input_filepaths = [
        os.path.join(infolder, "C11.bin"), 
        os.path.join(infolder, "C12_real.bin"),
        os.path.join(infolder, "C12_imag.bin"),
        os.path.join(infolder, "C22.bin")
    ]

    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "Ps_mf3cc.tif"))
        output_filepaths.append(os.path.join(infolder, "Pd_mf3cc.tif"))
        output_filepaths.append(os.path.join(infolder, "Pv_mf3cc.tif"))
        output_filepaths.append(os.path.join(infolder, "Theta_CP_mf3cc.tif"))
        

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, 
                write_flag=write_flag,processing_func=process_chunk_mf3cc,block_size=(512, 512), 
                max_workers=max_workers,  num_outputs=4, chi_in=chi_in,psi_in=psi_in)

def process_chunk_mf3cc(chunks, window_size,input_filepaths,chi_in,psi_in):

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
    
    c2_det = (c11_T1*c22_T1-c12_T1*c21_T1)
    c2_trace = c11_T1+c22_T1
    # t2_span = t11s*t22s
    m1 = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))

    # Compute Stokes parameters
    s0 = c11_T1 + c22_T1
    s1 = c11_T1 - c22_T1
    s2 = np.real(c12_T1 + c21_T1)
    s3 = np.where(chi_in >= 0, 1j * (c12_T1 - c21_T1), -1j * (c12_T1 - c21_T1))
    s3 = np.real(s3)

    SC = ((s0)-(s3))/2;
    OC = ((s0)+(s3))/2;

    h = (OC-SC)
    # span = c11s + c22s

    val = ((m1*s0*h))/((SC*OC + (m1**2)*(s0**2)))
    thet = np.real(np.arctan(val))
    theta_CP = np.rad2deg(thet)

    Ps_CP= (((m1*(c2_trace)*(1.0+np.sin(2*thet))/2)))
    Pd_CP= (((m1*(c2_trace)*(1.0-np.sin(2*thet))/2)))
    Pv_CP= (c2_trace*(1.0-m1))
    

    return Ps_CP, Pd_CP, Pv_CP, theta_CP