from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d,load_C2


def eig22(c2):
    c11 = c2[:,:,0].flatten()
    c12 = c2[:,:,1].flatten()
    c21 = c2[:,:,2].flatten()
    c22 = c2[:,:,3].flatten()
    trace = -(c11+c22)
    det = c11*c22-c12*c21
    # const= 1
    sqdiscr = np.sqrt(trace*trace - 4*det);
    lambda1 = -(trace + sqdiscr)*0.5;
    lambda2 = -(trace - sqdiscr)*0.5;
    
    return lambda1,lambda2

def prvidp(C2_folder,window_size=1,write_flag=None):

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
    m = (np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))
    prvi = (1-m)*c22s

    if write_flag:
        infile = C2_folder+'/C11.bin'
        """Write files to disk"""
        if os.path.exists(C2_folder+'/C11.bin'):
            infile = C2_folder+'/C11.bin'

        ofile = C2_folder+'/PRVI.bin'
        write_bin(ofile,dprvi,infile)
                    
    return prvi
