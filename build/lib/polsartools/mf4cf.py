from osgeo import gdal
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d

def mf4cf(T3_folder,window_size=1,write_flag=None):
        
    T11 = read_bin(T3_folder+"/T11.bin").astype(np.float32)
    T22 = read_bin(T3_folder+"/T22.bin").astype(np.float32)
    T33 = read_bin(T3_folder+"/T33.bin").astype(np.float32)

    T12_i = read_bin(T3_folder+'/T12_imag.bin').astype(np.float32)
    T12_r = read_bin(T3_folder+'/T12_real.bin').astype(np.float32)
    T13_i = read_bin(T3_folder+'/T13_imag.bin').astype(np.float32)
    T13_r = read_bin(T3_folder+'/T13_real.bin').astype(np.float32)
    T23_i = read_bin(T3_folder+'/T23_imag.bin').astype(np.float32)
    T23_r = read_bin(T3_folder+'/T23_real.bin').astype(np.float32)

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

    t12_T1r = conv2d(np.real(t12_T1),kernel)
    t12_T1i = conv2d(np.imag(t12_T1),kernel)
    t12s = t12_T1r+1j*t12_T1i

    t13_T1r = conv2d(np.real(t13_T1),kernel)
    t13_T1i = conv2d(np.imag(t13_T1),kernel)
    t13s = t13_T1r+1j*t13_T1i


    t21_T1r = conv2d(np.real(t21_T1),kernel)
    t21_T1i = conv2d(np.imag(t21_T1),kernel)
    t21s = t21_T1r+1j*t21_T1i

    t22_T1r = conv2d(np.real(t22_T1),kernel)
    t22_T1i = conv2d(np.imag(t22_T1),kernel)
    t22s = t22_T1r+1j*t22_T1i

    t23_T1r = conv2d(np.real(t23_T1),kernel)
    t23_T1i = conv2d(np.imag(t23_T1),kernel)
    t23s = t23_T1r+1j*t23_T1i
    

    
    t31_T1r = conv2d(np.real(t31_T1),kernel)
    t31_T1i = conv2d(np.imag(t31_T1),kernel)
    t31s = t31_T1r+1j*t31_T1i

    t32_T1r = conv2d(np.real(t32_T1),kernel)
    t32_T1i = conv2d(np.imag(t32_T1),kernel)
    t32s = t32_T1r+1j*t32_T1i

    t33_T1r = conv2d(np.real(t33_T1),kernel)
    t33_T1i = conv2d(np.imag(t33_T1),kernel)
    t33s = t33_T1r+1j*t33_T1i

    det_T3 = t11s*(t22s*t33s-t23s*t32s)-t12s*(t21s*t33s-t23s*t31s)+t13s*(t21s*t32s-t22s*t31s)
    trace_T3 = t11s + t22s + t33s
    m1 = np.real(np.sqrt(1-(27*(det_T3/(trace_T3**3)))))

    k11_f = (t11s + t22s + t33s)/2
    k44_f = (-t11s + t22s + t33s)/2
    k14_f = np.imag(t23s)

    trace_T3 = t11s + t22s + t33s
 

    s0_f = trace_T3
    dop_f = m1

    val1 = (4*dop_f*k11_f*k44_f)/(k44_f**2 - (1 + 4*dop_f**2)*k11_f**2)
    val2 = np.abs(k14_f)/(k11_f)
    


    theta_f = np.real(np.arctan(val1)) # separation for surface and dbl
    tau_f = np.real(np.arctan(val2)) # separation for helix
    # thet = np.rad2deg(thet)
    theta_FP = np.rad2deg(theta_f).astype(np.float32)
    tau_FP = np.rad2deg(tau_f).astype(np.float32)

    

    pc_f = (dop_f*s0_f*(np.sin(2*tau_f))).astype(np.float32)
    pv_f = ((1-dop_f)*s0_f).astype(np.float32)
    res_pow = s0_f - (pc_f + pv_f)
    ps_f = ((res_pow/2)*(1+np.sin((2*theta_f)))).astype(np.float32)
    pd_f = ((res_pow/2)*(1-np.sin((2*theta_f)))).astype(np.float32)
    
    if write_flag:
        infile = T3_folder+'/T11.bin'
        """Write files to disk"""
        if os.path.exists(T3_folder+'/T11.bin'):
            infile = iFolder+'/T11.bin'
        elif os.path.exists(T3_folder+'/C11.bin'):
            infile = T3_folder+'/C11.bin'

        ofilegrvi0 = T3_folder+'/Tau_FP_4c.bin'
        write_bin(ofilegrvi0,tau_FP,infile)

        ofilegrvi = T3_folder+'/Theta_FP_4c.bin'
        write_bin(ofilegrvi,theta_FP,infile)
        
        ofilegrvi1 = T3_folder+'/Pd_FP_4c.bin'
        write_bin(ofilegrvi1,pd_f,infile)
        
        ofilegrvi2 = T3_folder+'/Ps_FP_4c.bin'
        write_bin(ofilegrvi2,ps_f,infile)
        
        ofilegrvi3 = T3_folder+'/Pv_FP_4c.bin'
        write_bin(ofilegrvi3,pv_f,infile)

        ofilegrvi4 = T3_folder+'/Pc_FP_4c.bin'
        write_bin(ofilegrvi4,pc_f,infile)     

    return ps_f, pd_f, pv_f,pc_f,theta_FP,tau_FP