import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel,conv2d,eig22, time_it

@time_it
def misomega(infolder, outname=None, chi_in=45, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    input_filepaths = [
        os.path.join(infolder, "C11.bin"), 
        os.path.join(infolder, "C12_real.bin"),
        os.path.join(infolder, "C12_imag.bin"),
        os.path.join(infolder, "C22.bin")
    ]

    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "Ps_miSOmega.tif"))
        output_filepaths.append(os.path.join(infolder, "Pd_miSOmega.tif"))
        output_filepaths.append(os.path.join(infolder, "Pv_miSOmega.tif"))
        

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, 
                write_flag=write_flag,processing_func=process_chunk_misomega,block_size=(512, 512), 
                max_workers=max_workers,  num_outputs=3, chi_in=chi_in,psi_in=psi_in)

def process_chunk_misomega(chunks, window_size,input_filepaths,chi_in,psi_in):

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
    s0 = np.abs(c11_T1 + c22_T1)
    s1 = np.abs(c11_T1 - c22_T1)
    s2 = np.abs(c12_T1 + c21_T1)
    s3 = np.where(chi_in >= 0, 1j * (c12_T1 - c21_T1), -1j * (c12_T1 - c21_T1))
    s3 = np.real(s3)

    ## Stokes child parameters
    SC = ((s0)-(s3))/2;
    OC = ((s0)+(s3))/2;
    
    CPR = np.divide(SC,OC)  ##SC/OC
    # CPR = SC/OC
    ##scattered fields    
    dop= np.sqrt(np.power(s1,2) + np.power(s2,2) + np.power(s3,2))/(s0)
    Psi = 0.5*((180/np.pi)*np.arctan2(s2,s1))
    DOCP = (-s3)/(dop*s0);
    Chi = 0.5*((180/np.pi)*np.arcsin(DOCP))
    ##---------------------------------
    
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
    

    return surface_new, double_bounce_new, diffused_new