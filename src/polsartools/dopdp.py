from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d

def load_C2(folder):

    C11 = read_bin(folder+"/C11.bin")
    C22 = read_bin(folder+"/C22.bin")

    C12_i = read_bin(folder+'/C12_imag.bin')
    C12_r = read_bin(folder+'/C12_real.bin')

    C12 = C12_r + 1j*C12_i

    return np.dstack((C11,C12,np.conj(C12),C22))

def dopdp(C2_folder,window_size=1,write_flag=None):

    C2_stack = load_C2(C2_folder)

    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = C2_stack[:,:,0]
    c12_T1 = C2_stack[:,:,1]
    c21_T1 = C2_stack[:,:,2]
    c22_T1 = C2_stack[:,:,3]

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
    # t2_span = t11s*t22s
    dop = (np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))


    if write_flag:
        infile = C2_folder+'/C11.bin'
        """Write files to disk"""
        if os.path.exists(C2_folder+'/C11.bin'):
            infile = C2_folder+'/C11.bin'

        ofile = C2_folder+'/DOP_DP.bin'
        write_bin(ofile,dop,infile)
                    

    return dop
