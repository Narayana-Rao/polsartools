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

def mf3cc(C2_folder,chi_in=45,window_size=1,write_flag=None):

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
    m1 = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))

    # Stokes Parameter
    s0 = c11s + c22s;
    s1 = c11s - c22s;
    s2 = (c12s + c21s);

    if (chi_in >= 0):
        s3 = (1j*(c12s - c21s)); # The sign is according to RC or LC sign !!
    if (chi_in < 0):
        s3 = -(1j*(c12s - c21s)); # The sign is according to RC or LC sign !!
    
    SC = ((s0)-(s3))/2;
    OC = ((s0)+(s3))/2;

    h = (OC-SC)
    span = c11s + c22s

    val = ((m1*s0*h))/((SC*OC + (m1**2)*(s0**2)))
    thet = np.real(np.arctan(val))
    theta_CP = np.rad2deg(thet)

    Ps_CP= (((m1*(span)*(1.0+np.sin(2*thet))/2)))
    Pd_CP= (((m1*(span)*(1.0-np.sin(2*thet))/2)))
    Pv_CP= (span*(1.0-m1))


    if write_flag:
        infile = C2_folder+'/C11.bin'
        """Write files to disk"""
        if os.path.exists(C2_folder+'/C11.bin'):
            infile = C2_folder+'/C11.bin'

        ofile = C2_folder+'/Theta_CP.bin'
        write_bin(ofile,theta_CP,infile)
        ofile1 = C2_folder+'/Pd_CP.bin'
        write_bin(ofile1,Pd_CP,infile)
        ofile2 = C2_folder+'/Ps_CP.bin'
        write_bin(ofile2,Ps_CP,infile)
        ofile3 = C2_folder+'/Pv_CP.bin'
        write_bin(ofile3,Pv_CP,infile)   
                    

    return Ps_CP, Pd_CP, Pv_CP, theta_CP
