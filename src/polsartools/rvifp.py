from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d


def grvi(T3_folder,window_size=1,write_flag=None):

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
                
    nrows  = np.shape(T11)[1]
    ncols = np.shape(T11)[0]
    # nrows  = 100
    # ncols = 100
   
    temp_rvi = np.zeros((ncols,nrows))
                
    # %% for window processing
    wsi=wsj=ws
    
    inci=int(np.fix(wsi/2)) # Up & down movement margin from the central row
    incj=int(np.fix(wsj/2)) # Left & right movement from the central column
    # % Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999
    
    starti=int(np.fix(wsi/2)) # Starting row for window processing
    startj=int(np.fix(wsj/2)) # Starting column for window processing
    
    stopi= int(nrows-inci)-1 # Stop row for window processing
    stopj= int(ncols-incj)-1 # Stop column for window processing
            
    for ii in np.arange(startj,stopj+1):

        # self.progress.emit(str(ii)+'/'+str(nrows))
        self.pBar.emit(int((ii/ncols)*90))
        for jj in np.arange(starti,stopi+1):
    
            t11s = np.nanmean(t11_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            t12s = np.nanmean(t12_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            t13s = np.nanmean(t13_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            
            t21s = np.nanmean(t21_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            t22s = np.nanmean(t22_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            t23s = np.nanmean(t23_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            
            t31s = np.nanmean(t31_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            t32s = np.nanmean(t32_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            t33s = np.nanmean(t33_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
    
            T_T1 = np.array([[t11s, t12s, t13s], [t21s, t22s, t23s], [t31s, t32s, t33s]])
    
            # %% RVI
            if np.isnan(np.real(T_T1)).any() or np.isinf(np.real(T_T1)).any() or np.isneginf(np.real(T_T1)).any():
                T_T1 = np.array([[0,0],[0,0]])
                temp_rvi[ii,jj] = 0
                # self.progress.emit(str('invalid Value encountered!!'))
                continue
                
            e_v = -np.sort(-np.linalg.eigvals(T_T1)); # sorting in descending order
            e_v1 = e_v[0]; e_v2 = e_v[1]; e_v3 = e_v[2];

            # self.progress.emit(str('Eigen val Done'))
            
            p1 = e_v1/(e_v1 + e_v2 + e_v3);
            p2 = e_v2/(e_v1 + e_v2 + e_v3);
            p3 = e_v3/(e_v1 + e_v2 + e_v3);
            
            p1=0 if p1<0 else p1
            p2=0 if p2<0 else p2
            p3=0 if p3<0 else p3
                
            p1=1 if p1>1 else p1
            p2=1 if p2>1 else p2
            p3=1 if p3>1 else p3
            
            
            
            temp_rvi[ii,jj] = np.real((4*p3)/(p1 + p2 + p3));
    
    # %% RVI scaled (0 - 1)   
    rvi = temp_rvi;   
    idx = np.argwhere(rvi>1)

    rvi[idx] = (3/4)*rvi[idx];
    rvi[~idx] = rvi[~idx];
    rvi[rvi==0] = np.NaN


    if write_flag:
        infile = T3_folder+'/T11.bin'
        """Write files to disk"""
        if os.path.exists(T3_folder+'/T11.bin'):
            infile = iFolder+'/T11.bin'
        elif os.path.exists(T3_folder+'/C11.bin'):
            infile = T3_folder+'/C11.bin'

        ofilegrvi = T3_folder+'/RVIFP.bin'
        write_bin(ofilegrvi,vi,infile)
                    

    return rvi
