from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d, load_C2


def dopcp(C2_folder,chi_in=45,window_size=1,write_flag=None):

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

    # Stokes Parameter
    s0 = c11s + c22s;
    s1 = c11s - c22s;
    s2 = (c12s + c21s);

    if (chi_in >= 0):
        s3 = (1j*(c12s - c21s)); # The sign is according to RC or LC sign !!
    if (chi_in < 0):
        s3 = -(1j*(c12s - c21s)); # The sign is according to RC or LC sign !!
    

    dop= np.sqrt(np.power(s1,2) + np.power(s2,2) + np.power(s3,2))/(s0);   

    if write_flag:
        infile = os.path.join(C2_folder,'C11.bin')
        """Write files to disk"""
        if os.path.exists(os.path.join(C2_folder,'C11.bin')):
            infile = os.path.join(C2_folder,'C11.bin')

        ofile = os.path.join(C2_folder,'DOP_CP.bin')
        write_bin(ofile,dop,infile)
                    

    return np.real(dop)
