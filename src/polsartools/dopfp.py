from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d

def dopfp(T3_folder,window_size=1,dop_type=1,write_flag=None):
	T11 = read_bin(T3_folder+"/T11.bin")
	T22 = read_bin(T3_folder+"/T22.bin")
	T33 = read_bin(T3_folder+"/T33.bin")

	T12_i = read_bin(T3_folder+'/T12_imag.bin')
	T12_r = read_bin(T3_folder+'/T12_real.bin')
	T13_i = read_bin(T3_folder+'/T13_imag.bin')
	T13_r = read_bin(T3_folder+'/T13_real.bin')
	T23_i = read_bin(T3_folder+'/T23_imag.bin')
	T23_r = read_bin(T3_folder+'/T23_real.bin')

	T12 = T12_r + 1j*T12_i
	T13 = T13_r + 1j*T13_i
	T23 = T23_r + 1j*T23_i
	    
	t11_T1 = T11
	t12_T1 = T12
	t13_T1 = T13
	t21_T1 = np.conj(T12)
	t22_T1 = T22
	t23_T1 = T23
	t31_T1 = np.conj(T13)
	t32_T1 = np.conj(T23)
	t33_T1 = T33
	del T11, T22, T33,T12, T13, T23, T12_i,T12_r, T13_i, T13_r, T23_i, T23_r
	            
	kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

	        
    t11_T1r = conv2d(np.real(t11_T1),kernel)
    t11_T1i = conv2d(np.imag(t11_T1),kernel)
    t11s = t11_T1r+1j*t11_T1i
    del t11_T1r,t11_T1i

    t12_T1r = conv2d(np.real(t12_T1),kernel)
    t12_T1i = conv2d(np.imag(t12_T1),kernel)
    t12s = t12_T1r+1j*t12_T1i
    del t12_T1r,t12_T1i

    t13_T1r = conv2d(np.real(t13_T1),kernel)
    t13_T1i = conv2d(np.imag(t13_T1),kernel)
    t13s = t13_T1r+1j*t13_T1i
    del t13_T1r,t13_T1i

    t21_T1r = conv2d(np.real(t21_T1),kernel)
    t21_T1i = conv2d(np.imag(t21_T1),kernel)
    t21s = t21_T1r+1j*t21_T1i
    del t21_T1r,t21_T1i

    t22_T1r = conv2d(np.real(t22_T1),kernel)
    t22_T1i = conv2d(np.imag(t22_T1),kernel)
    t22s = t22_T1r+1j*t22_T1i
    del t22_T1r,t22_T1i

    t23_T1r = conv2d(np.real(t23_T1),kernel)
    t23_T1i = conv2d(np.imag(t23_T1),kernel)
    t23s = t23_T1r+1j*t23_T1i
    del t23_T1r,t23_T1i

    t31_T1r = conv2d(np.real(t31_T1),kernel)
    t31_T1i = conv2d(np.imag(t31_T1),kernel)
    t31s = t31_T1r+1j*t31_T1i
    del t31_T1r,t31_T1i

    t32_T1r = conv2d(np.real(t32_T1),kernel)
    t32_T1i = conv2d(np.imag(t32_T1),kernel)
    t32s = t32_T1r+1j*t32_T1i
    del t32_T1r,t32_T1i

    t33_T1r = conv2d(np.real(t33_T1),kernel)
    t33_T1i = conv2d(np.imag(t33_T1),kernel)
    t33s = t33_T1r+1j*t33_T1i
    del t33_T1r,t33_T1i

	if dop_type==0:
		dop_fp=0

	# Barakat DOP
	if dop_type==1:            
	    det_T3 = t11s*(t22s*t33s-t23s*t32s)-t12s*(t21s*t33s-t23s*t31s)+t13s*(t21s*t32s-t22s*t31s)
	    trace_T3 = t11s + t22s + t33s
	    dop_fp = np.sqrt(1-((27*det_T3)/trace_T3**3))

	if write_flag:
	    infile = T3_folder+'/T11.bin'
	    """Write files to disk"""
	    if os.path.exists(T3_folder+'/T11.bin'):
	        infile = iFolder+'/T11.bin'
	    elif os.path.exists(T3_folder+'/C11.bin'):
	        infile = T3_folder+'/C11.bin'

	    ofilegrvi = T3_folder+'/dop_fp.bin'
	    write_bin(ofilegrvi,dop_fp,infile)

	return dop_fp
