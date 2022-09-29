from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d, load_T2



def mf3cd(T2_folder,window_size=1,write_flag=None):

    T2_stack = load_T2(T2_folder)

    t11_T1 = T2_stack[:,:,0]
    t12_T1 = T2_stack[:,:,1]
    t21_T1 = T2_stack[:,:,2]
    t22_T1 = T2_stack[:,:,3]
    
    
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

    t11_T1r = conv2d(np.real(t11_T1),kernel)
    t11_T1i = conv2d(np.imag(t11_T1),kernel)
    t11s = t11_T1r+1j*t11_T1i

    t12_T1r = conv2d(np.real(t12_T1),kernel)
    t12_T1i = conv2d(np.imag(t12_T1),kernel)
    t12s = t12_T1r+1j*t12_T1i

    t21_T1r = conv2d(np.real(t21_T1),kernel)
    t21_T1i = conv2d(np.imag(t21_T1),kernel)
    t21s = t21_T1r+1j*t21_T1i

    t22_T1r = conv2d(np.real(t22_T1),kernel)
    t22_T1i = conv2d(np.imag(t22_T1),kernel)
    t22s = t22_T1r+1j*t22_T1i
    
    det_T2 = t11s*t22s-t12s*t21s
    trace_T2 = t11s + t22s

    m1 = np.real(np.sqrt(1-(4*(det_T2/(trace_T2**2)))))
    h = (t11s - t22s)
    g = t22s
    span = t11s + t22s
    
    val = (m1*span*h)/(t11s*g+m1**2*span**2);
    thet = np.real(np.arctan(val))
    # thet = np.rad2deg(thet)
    theta_DP = np.rad2deg(thet)
    
    Ps_DP = np.real((((m1*(span)*(1+np.sin(2*thet))/2))))
    Pd_DP = np.real((((m1*(span)*(1-np.sin(2*thet))/2))))
    Pv_DP = np.real((span*(1-m1))) 

    if write_flag:
        infile = T2_folder+'/T11.bin'
        """Write files to disk"""
        if os.path.exists(T2_folder+'/T11.bin'):
            infile = T2_folder+'/T11.bin'

        ofile = T2_folder+'/Theta_DP.bin'
        write_bin(ofile,theta_DP,infile)
        ofile1 = T2_folder+'/Pd_DP.bin'
        write_bin(ofile1,Pd_DP,infile)
        
        ofile2 = T2_folder+'/Ps_DP.bin'
        write_bin(ofile2,Ps_DP,infile)
        
        ofile3 = T2_folder+'/Pv_DP.bin'
        write_bin(ofile3,Pv_DP,infile)      

    return Ps_DP, Pd_DP, Pv_DP, theta_DP
