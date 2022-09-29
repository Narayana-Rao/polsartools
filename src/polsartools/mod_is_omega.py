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

def mod_is_omega(C2_folder,chi_in=45,psi_in=0,window_size=1,write_flag=None):

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
    s0 = np.float32(np.real(c11s + c22s))
    s1 = np.float32(np.real(c11s - c22s))
    s2 = np.float32(np.real(c12s + c21s))

    if (chi_in >= 0):
        s3 = np.float32(np.real(1j*(c12s - c21s))) # The sign is according to RC or LC sign !!
    if (chi_in < 0):
        s3 = np.float32(np.real(-(1j*(c12s - c21s)))) # The sign is according to RC or LC sign !!
    
    ## Stokes child parameters
    SC = ((s0)-(s3))/2;
    OC = ((s0)+(s3))/2;
    #old_err_state = np.seterr(divide='raise')
    #ignored_states = np.seterr(**old_err_state)
    CPR = np.divide(SC,OC)  ##SC/OC
    # CPR = SC/OC
    ##scattered fields    
    dop= np.sqrt(np.power(s1,2) + np.power(s2,2) + np.power(s3,2))/(s0)
    Psi = 0.5*((180/np.pi)*np.arctan2(s2,s1))
    DOCP = (-s3)/(dop*s0);
    Chi = 0.5*((180/np.pi)*np.arcsin(DOCP))
    ##---------------------------------
    #psi_in = self.tau
    
    ##---------------------------------
    # Calculating Omega from S-Omega decomposition        
    x1 = np.cos(2*chi_in*np.pi/180)*np.cos(2*psi_in*np.pi/180)*np.cos(2*Chi*np.pi/180)*np.cos(2*Psi*np.pi/180)
    x2 = np.cos(2*chi_in*np.pi/180)*np.sin(2*psi_in*np.pi/180)*np.cos(2*Chi*np.pi/180)*np.sin(2*Psi*np.pi/180)
    x3 = np.abs(np.sin(2*chi_in*np.pi/180)*np.sin(2*Chi*np.pi/180))
    Prec =  dop*(1 + x1 + x2 + x3)
    Prec1 = (1 - dop) + dop*(1 + x1 + x2 + x3)
    omega = (Prec/Prec1)
                      
    
    # ## Improved S-Omega (i-SOmega powers
    ind_g1 = (CPR>1).astype(int)
    s_new_g1 = omega*(1 - omega)*OC
    db_new_g1 = omega*s0 - omega*(1 - omega)*OC   ##depolarized of OC x polarized of SC

    ind_l1 = (CPR<1).astype(int)
    s_new_l1 = omega*s0 - omega*(1 - omega)*SC
    db_new_l1 = omega*(1 - omega)*SC   ##depolarized of OC x polarized of SC

    ind_e1 = (CPR==1).astype(int)
    s_new_e1 = omega*OC
    db_new_e1 = omega*SC

    surface_new = s_new_g1*ind_g1+s_new_l1*ind_l1+s_new_e1*ind_e1
    double_bounce_new = db_new_g1*ind_g1+db_new_l1*ind_l1+db_new_e1*ind_e1

    
    diffused_new = (1 - omega)*s0; ##diffused scattering                       
    
    surface_new[surface_new==0] = np.nan
    double_bounce_new[double_bounce_new==0] = np.nan
    diffused_new[diffused_new==0] = np.nan


    if write_flag:
        infile = C2_folder+'/C11.bin'
        """Write files to disk"""
        if os.path.exists(C2_folder+'/C11.bin'):
            infile = C2_folder+'/C11.bin'

        ofileps = C2_folder+'/Ps_iSOmega.bin'
        write_bin(ofileps,surface_new,infile)
        
        ofilepd = C2_folder+'/Pd_iSOmega.bin'
        write_bin(ofilepd,double_bounce_new,infile)
        
        ofilepv = C2_folder+'/Pv_iSOmega.bin'
        write_bin(ofilepv,diffused_new,infile)
        

    return ofileps, ofilepd, ofilepv
