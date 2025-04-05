import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel,conv2d,eig22, time_it
from polsartools.cprvicpp import process_chunk_cprvicpp
from polsartools.testcprvi import sum_filt

import polsartools
import traceback
import pickle
@time_it
def cprvi(infolder, outname=None, chi_in=45, psi_in=0,window_size=1,write_flag=True,max_workers=None):
    input_filepaths = [
        os.path.join(infolder, "C11.bin"), 
        os.path.join(infolder, "C12_real.bin"),
        os.path.join(infolder, "C12_imag.bin"),
        os.path.join(infolder, "C22.bin")
    ]

    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "cprvi.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, 
                write_flag=write_flag,processing_func=process_chunk_cprvi,block_size=(512, 512), 
                max_workers=max_workers,  num_outputs=1, chi_in=chi_in,psi_in=psi_in)

def process_chunk_cprvi(chunks, window_size,input_filepaths,chi_in,psi_in):

    chunk_arrays = [np.array(ch) for ch in chunks]  
    vi_c_raw = process_chunk_cprvicpp( chunk_arrays, window_size, input_filepaths, chi_in, psi_in )
    # vi_c_raw = sum_filt(chunk_arrays, window_size, input_filepaths, chi_in, psi_in )
    
    return np.array(vi_c_raw, copy=True)  
    
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    ncols,nrows = np.shape(c11_T1)

    # if window_size>1:
    #     c11_T1 = conv2d(np.real(c11_T1),kernel)+1j*conv2d(np.imag(c11_T1),kernel)
    #     c12_T1 = conv2d(np.real(c12_T1),kernel)+1j*conv2d(np.imag(c12_T1),kernel)
    #     c21_T1 = conv2d(np.real(c21_T1),kernel)+1j*conv2d(np.imag(c21_T1),kernel)
    #     c22_T1 = conv2d(np.real(c22_T1),kernel)+1j*conv2d(np.imag(c22_T1),kernel)

    # # Compute Stokes parameters
    # s0 = c11_T1 + c22_T1
    # s1 = c11_T1 - c22_T1
    # s2 = c12_T1 + c21_T1
    # s3 = np.where(chi_in >= 0, 1j * (c12_T1 - c21_T1), -1j * (c12_T1 - c21_T1))

    #     # Compute K_T matrix elements
    # K_T = np.zeros((ncols, nrows, 4, 4), dtype=complex)
    # K_T[:, :, 0, 0] = s0
    # K_T[:, :, 0, 2] = s2
    # K_T[:, :, 1, 3] = s1
    # K_T[:, :, 3, 3] = s3

    # # GD_VI - VV-VH/VV-HH (Depolarization calculation)
    # SC = (s0 - s3) / 2
    # OC = (s0 + s3) / 2

    # min_sc_oc = np.real(np.minimum(SC, OC))
    # max_sc_oc = np.real(np.maximum(SC, OC))

    # # Handle edge cases to avoid division by zero or very small values
    # max_sc_oc = np.where(max_sc_oc == 0, 1e-10, max_sc_oc)  # Avoid division by zero
    
    # fp22 = min_sc_oc / max_sc_oc

    # K_depol = np.zeros((ncols, nrows, 4, 4))
    # K_depol[:, :, 0, 0] = 1

    # # GD_DEPOL Calculation
    # num = np.einsum('ijkl,ijlk->ij', K_T, K_depol)
    # den1 = np.sqrt(np.abs(np.einsum('ijkl,ijlk->ij', K_T, K_T)))
    # den2 = np.sqrt(np.abs(np.einsum('ijkl,ijlk->ij', K_depol, K_depol)))
    # den = den1 * den2

    # # Handle division by zero in depolarization calculation
    # den = np.where(den == 0, 1e-10, den)  # Avoid division by zero

    # GD_t1_depol = np.real(2 * np.arccos(np.clip(num / den, -1, 1)) * 180 / np.pi / 180)
    # l_lambda = (3 / 2) * GD_t1_depol

    # # Final calculation of vi_c
    # vi_c = np.real((1 - l_lambda) * np.power(fp22, 2 * l_lambda))

    # return vi_c
    
    fp22 = np.zeros((ncols,nrows))
    l_lambda = np.zeros((ncols,nrows))

    wsi=wsj=window_size

    inci=int(np.fix(wsi/2)) # Up & down movement margin from the central row
    incj=int(np.fix(wsj/2)) # Left & right movement from the central column
    # % Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

    starti=int(np.fix(wsi/2)) # Starting row for window processing
    startj=int(np.fix(wsj/2)) # Starting column for window processing

    stopi= int(nrows-inci)-1 # Stop row for window processing
    stopj= int(ncols-incj)-1 # Stop column for window processing

    for ii in np.arange(startj,stopj+1):

        for jj in np.arange(starti,stopi+1):
            
            C11c = np.nanmean(c11_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            C12c = np.nanmean(c12_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            C21c = np.nanmean(c21_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            C22c = np.nanmean(c22_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample

            C0 = np.array([[C11c,C12c], [C21c, C22c]]);
            
            # %% GD_VI -- VV-VH/VV-HH
            if np.isnan(np.real(C0)).any() or np.isinf(np.real(C0)).any() or np.isneginf(np.real(C0)).any():
                C0 = np.array([[0,0],[0,0]])

            # Stokes Parameter
            s0 = C11c + C22c;
            s1 = C11c - C22c;
            s2 = (C12c + C21c);

            if (chi_in >= 0):
                s3 = (1j*(C12c - C21c)); # The sign is according to RC or LC sign !!
            if (chi_in < 0):
                s3 = -(1j*(C12c - C21c)); # The sign is according to RC or LC sign !!
            
            k11 = s0
            k12 = 0
            k13 = s2
            k14 = 0
            k21 = k12 
            k22 = 0
            k23 = 0
            k24 = s1
            k31 = k13
            k32 = k23
            k33 = 0 
            k34 = 0;
            k41 = k14; 
            k42 = k24; 
            k43 = k34; 
            k44 = s3;

            K_T = 0.5*np.array([[k11,k12,k13,k14], [k21,k22,k23,k24], 
                [k31, k32, k33, k34], [k41,k42,k43,k44]])       

            # Stokes vector child products
            SC = ((s0)-(s3))/2;
            OC = ((s0)+(s3))/2;
            
            min_sc_oc = min(SC,OC);
            max_sc_oc = max(SC,OC);                
            
            K_depol = np.array([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]);
            
            # GD_DEPOL
            
            num1 = np.matmul(K_T.T,K_depol);
            num = np.trace(num1);
            den1 = np.sqrt(abs(np.trace(np.matmul(K_T.T,K_T))));
            den2 = np.sqrt(abs(np.trace(np.matmul(K_depol.T,K_depol))));
            den = den1*den2;
            
            temp_aa = np.real(2*np.arccos(num/den)*180/np.pi);
            GD_t1_depol = np.real(temp_aa/180);
            
                                                           
            l_lambda[ii,jj] = (3/2)*GD_t1_depol;
            
            #GD_VI -- RH-RV/LH-LV
            
            fp22[ii,jj] = (min_sc_oc/max_sc_oc);


    vi_c = np.real((1 - l_lambda)*np.power(fp22, 2*l_lambda))

    return vi_c