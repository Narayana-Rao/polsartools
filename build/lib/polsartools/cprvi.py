from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d, load_C2

def cprvi(C2_folder,chi_in=45,window_size=1,write_flag=None):

    C2_stack = load_C2(C2_folder)

    nrows  = np.shape(C2_stack)[1]
    ncols = np.shape(C2_stack)[0]

    C11 = C2_stack[:,:,0]
    C12 = C2_stack[:,:,1]
    C21 = C2_stack[:,:,2]
    C22 = C2_stack[:,:,3]

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
            
            C11c = np.nanmean(C11[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            C12c = np.nanmean(C12[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            C21c = np.nanmean(C21[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            C22c = np.nanmean(C22[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample

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

    if write_flag:
        infile = os.path.join(C2_folder,'C11.bin')
        """Write files to disk"""
        if os.path.exists(os.path.join(C2_folder,'C11.bin')):
            infile = os.path.join(C2_folder,'C11.bin')

        ofile = os.path.join(C2_folder,'CpRVI.bin')
        write_bin(ofile,vi_c,infile)
                    

    return vi_c
