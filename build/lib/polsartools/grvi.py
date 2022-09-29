from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d


def grvi(T3_folder,window_size=1,write_flag=None):

    T11 = read_bin(os.path.join(T3_folder,"T11.bin"))
    T22 = read_bin(os.path.join(T3_folder,"T22.bin"))
    T33 = read_bin(os.path.join(T3_folder,"T33.bin"))

    T12_i = read_bin(os.path.join(T3_folder,'T12_imag.bin'))
    T12_r = read_bin(os.path.join(T3_folder,'T12_real.bin'))
    T13_i = read_bin(os.path.join(T3_folder,'T13_imag.bin'))
    T13_r = read_bin(os.path.join(T3_folder,'T13_real.bin'))
    T23_i = read_bin(os.path.join(T3_folder,'T23_imag.bin'))
    T23_r = read_bin(os.path.join(T3_folder,'T23_real.bin'))

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

                
    nrows  = np.shape(T11)[1]
    ncols = np.shape(T11)[0]
    # nrows  = 100
    # ncols = 100
    del T11, T22, T33,T12, T13, T23, T12_i,T12_r, T13_i, T13_r, T23_i, T23_r

    span = np.zeros((ncols,nrows))
    rho_13_hhvv = np.zeros((ncols,nrows))
    # temp_rvi = np.zeros((ncols,nrows))
    fp22 = np.zeros((ncols,nrows))
    GD_t1_t = np.zeros((ncols,nrows))
    GD_t1_d = np.zeros((ncols,nrows))
    GD_t1_rv = np.zeros((ncols,nrows))
    GD_t1_nd = np.zeros((ncols,nrows))
    GD_t1_c = np.zeros((ncols,nrows))
    GD_t1_lh = np.zeros((ncols,nrows))
    GD_t1_rh = np.zeros((ncols,nrows))
    beta = np.zeros((ncols,nrows))
    beta_1 = np.zeros((ncols,nrows))
    f = np.zeros((ncols,nrows))
    a = np.zeros((ncols,nrows))
    b = np.zeros((ncols,nrows))
    temp_gamma = np.zeros((ncols,nrows))
    t_d = np.zeros((46,1))
    t_nd = np.zeros((46,1))
    t_t = np.zeros((46,1))
    t_c = np.zeros((46,1))
    theta_map = np.zeros((ncols,nrows))
    
    D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
    # %% for window processing
    wsi=wsj=window_size
    
    inci=int(np.fix(wsi/2)) # Up & down movement margin from the central row
    incj=int(np.fix(wsj/2)) # Left & right movement from the central column
    # % Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999
    
    starti=int(np.fix(wsi/2)) # Starting row for window processing
    startj=int(np.fix(wsj/2)) # Starting column for window processing
    
    stopi= int(nrows-inci)-1 # Stop row for window processing
    stopj= int(ncols-incj)-1 # Stop column for window processing
            # %% Elementary targets
    
    M_d = np.array([[1,0,0,0], [ 0,1,0,0], [ 0,0,-1,0], [ 0,0,0,1]])
    M_nd = np.array([[0.625,0.375,0,0], [ 0.375,0.625,0,0], [ 0,0,-0.5,0], [ 0,0,0,0.5]])
    M_t = np.array([[1,0,0,0], [ 0,1,0,0], [ 0,0,1,0], [ 0, 0,0,-1]])
    M_c = np.array([[0.625,0.375,0,0], [ 0.375,0.625,0,0], [0,0,0.5,0], [ 0,0,0,-0.5]])
    M_lh = np.array([[1,0,0,-1], [ 0,0,0,0], [ 0,0,0,0], [ -1,0,0,1]])
    M_rh = np.array([[1,0,0,1], [ 0,0,0,0], [ 0,0,0,0], [ 1,0,0,1]])
    
    for ii in np.arange(startj,stopj+1):

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
    
            #Coherency matrix
            C_T1 = np.matmul(np.matmul((D.T),T_T1),D);
            
            span[ii,jj] = np.real(t11s + t22s + t33s)
            temp_span = span[ii,jj]
            # self.progress.emit(str('span Done'))
            Temp_T1 = T_T1
            
            t11 = Temp_T1[0,0]; t12 = Temp_T1[0,1];        t13 = Temp_T1[0,2]
            t21 = np.conj(t12); t22 = Temp_T1[1,1];        t23 = Temp_T1[1,2]
            t31 = np.conj(t13); t32 = np.conj(t23);        t33 = Temp_T1[2,2]
            
            # %% Ratio of VV/HH (Used in Yamagichi Volume)
            
            hh = 0.5*(t11 + t22 + 2*np.real(t12)); #HH
            hv = t33; #HV
            vv = 0.5*(t11 + t22 - 2*np.real(t12)); #VV
            vol_con = 10*np.log10(vv/hh);
    
    
            # %% Kennaugh Matrix
            
            m11 = t11+t22+t33; m12 = t12+t21; m13 = t13+t31; m14 = -1j*(t23 - t32);
            m21 = t12+t21; m22 = t11+t22-t33; m23 = t23+t32; m24 = -1j*(t13-t31);
            m31 = t13+t31; m32 = t23+t32; m33 = t11-t22+t33; m34 = 1j*(t12-t21);
            m41 = -1j*(t23-t32); m42 = -1j*(t13-t31); m43 = 1j*(t12-t21); m44 = -t11+t22+t33;
            
            M_T = 0.5*np.array([[m11, m12, m13, m14], [m21, m22, m23, m24], [m31, m32, m33, m34], [m41, m42, m43, m44]]);
            
            
            M_T_theta = M_T;
            
            # %% GVSM
            
            t011 = M_T_theta[0,0] + M_T_theta[1,1] + M_T_theta[2,2] - M_T_theta[3,3];
            t012 = M_T_theta[0,1] - 1j*M_T_theta[2,3];
            t013 = M_T_theta[0,2] + 1j*M_T_theta[1,3];
            t021 = np.conj(t012);
            t022 = M_T_theta[0,0] + M_T_theta[1,1] - M_T_theta[2,2] + M_T_theta[3,3];
            t023 = M_T_theta[1,2] +1j*M_T_theta[0,3];
            t031 = np.conj(t013);
            t032 = np.conj(t023);
            t033 = M_T_theta[0,0] - M_T_theta[1,1] + M_T_theta[2,2] + M_T_theta[3,3];
            
            # %% T to C
            
            T0 = np.array([[t011/2, t012, t013], [t021, t022/2, t023], [t031, t032, t033/2]]);
            C0 = np.matmul(np.matmul((D.T),T0),D);
    
            # %% Gamma/Rho
            
            gamma = np.real(C0[0,0]/C0[2,2]); rho = 1/3;
            temp_gamma[ii,jj]= np.real(gamma); #% variable to save
            
            # %% Covariance matrix
            
            c11 = gamma; c12 = 0; c13 = rho*np.sqrt(gamma);
            c21 = 0; c22 = 0.5*(1 + gamma) - rho*np.sqrt(gamma); c23 = 0;
            c31 = np.conj(rho)*np.sqrt(gamma); c32 = 0; c33 = 1;
            
            R = (3/2)*(1 + gamma) - rho*np.sqrt(gamma);
            C1 = (1/R)*np.array([[c11, c12, c13], [c21, c22, c23], [c31, c32, c33]]);
            # self.progress.emit(str('gamma and R Done'))
            # %% Coherency matrix
            
            T1 = np.matmul(np.matmul(D,C1),(D.T));
            
            t11 = T1[0,0]; t12 = T1[0,1]; t13 = T1[0,2];
            t21 = T1[1,0]; t22 = T1[1,1]; t23 = T1[1,2];
            t31 = T1[2,0]; t32 = T1[2,1]; t33 = T1[2,2];
            
            m11 = t11+t22+t33; m12 = t12+t21; m13 = t13+t31; m14 = -1j*(t23 - t32);
            m21 = t12+t21; m22 = t11+t22-t33; m23 = t23+t32; m24 = -1j*(t13-t31);
            m31 = t13+t31; m32 = t23+t32; m33 = t11-t22+t33; m34 = 1j*(t12-t21);
            m41 = -1j*(t23-t32); m42 = -1j*(t13-t31); m43 = 1j*(t12-t21); m44 = -t11+t22+t33;
            
            # %% Generalized Random Volume (Antropov et al.)
            
            M_rv = np.real(np.array([[m11, m12, m13, m14], [m21, m22, m23, m24], [m31, m32, m33, m34], [m41, m42, m43, m44]]));
            
            f[ii,jj] = 1;
    
            # %% GD Volume
            
            num1 = np.matmul(((M_T_theta).T),M_rv); #% volume
            num = np.trace(num1);
            den1 = np.sqrt(abs(np.trace(np.matmul(((M_T_theta).T),M_T_theta))));
            den2 = np.sqrt(abs(np.trace(np.matmul(((M_rv).T),M_rv))));
            den = den1*den2;
            temp_aa = np.real(2*np.arccos(num/den)*180/np.pi);
            GD_t1_rv[ii,jj] = np.real(temp_aa/180);
            # self.progress.emit(str('GD volume Done'))
            # %% GD ALL
            
            num1 = np.matmul(((M_T_theta).T),M_c); #% cylinder
            num = np.trace(num1);
            den1 = np.sqrt(abs(np.trace(np.matmul(((M_T_theta).T),M_T_theta))));
            den2 = np.sqrt(abs(np.trace(np.matmul(((M_c).T),M_c))));
            den = den1*den2;
            temp_aa = np.real(2*np.arccos(num/den)*180/np.pi);
            GD_t1_c[ii,jj] = np.real(temp_aa/180);
            # self.progress.emit(str('GD cylider Done'))
            
            num1 = np.matmul(((M_T_theta).T),M_t); #% trihedral
            num = np.trace(num1);
            den1 = np.sqrt(abs(np.trace(np.matmul(((M_T_theta).T),M_T_theta))));
            den2 = np.sqrt(abs(np.trace(np.matmul(((M_t).T),M_t))));
            den = den1*den2;
            temp_aa = 2*np.arccos(num/den)*180/np.pi;
            GD_t1_t[ii,jj] = np.real(temp_aa/180);
            # self.progress.emit(str('GD trihedral Done'))
            
            num1 = np.matmul(((M_T_theta).T),M_d); #% dihedral
            num = np.trace(num1);
            den1 = np.sqrt(abs(np.trace(np.matmul(((M_T_theta).T),M_T_theta))));
            den2 = np.sqrt(abs(np.trace(np.matmul(((M_d).T),M_d))));
            den = den1*den2;
            temp_aa = 2*np.arccos(num/den)*180/np.pi;
            GD_t1_d[ii,jj] = np.real(temp_aa/180);
            # self.progress.emit(str('GD dihedral Done'))
            
            num1 = np.matmul(((M_T_theta).T),M_nd); #% n-dihedral
            num = np.trace(num1);
            den1 = np.sqrt(abs(np.trace(np.matmul(((M_T_theta).T),M_T_theta))));
            den2 = np.sqrt(abs(np.trace(np.matmul(((M_nd).T),M_nd))));
            den = den1*den2;
            temp_aa = 2*np.arccos(num/den)*180/np.pi;
            GD_t1_nd[ii,jj] = np.real(temp_aa/180);
            # self.progress.emit(str('GD n-dihedral Done'))
            
            # %% VI
            
            t_t = GD_t1_t[ii,jj];
            t_d = GD_t1_d[ii,jj];
            t_c = GD_t1_c[ii,jj];
            t_nd = GD_t1_nd[ii,jj];
            
            a[ii,jj] = np.nanmax([t_t, t_d, t_c, t_nd]);
            b[ii,jj] = np.nanmin([t_t, t_d, t_c, t_nd]);
            beta[ii,jj] = (b[ii,jj]/a[ii,jj])**2;
            
    # %% GRVI
    f[f==0]=np.NaN
    vi = np.power(beta, GD_t1_rv)*(1 - (1/f)*GD_t1_rv);
    
    x1 = np.power(beta, GD_t1_rv);
    x2 = (1 - GD_t1_rv);
    
    f =  np.nan_to_num(f)
    idx1 = np.argwhere(GD_t1_rv>f)
    vi[idx1] = 0;
    vi[~idx1] = vi[~idx1];
    

    if write_flag:
        infile = os.path.join(T3_folder,'T11.bin')
        """Write files to disk"""
        if os.path.exists(os.path.join(T3_folder,'T11.bin')):
            infile = os.path.join(iFolder,'T11.bin')
        elif os.path.exists(os.path.join(T3_folder,'C11.bin')):
            infile = os.path.join(T3_folder,'C11.bin')

        ofilegrvi = os.path.join(T3_folder,'GRVI.bin')
        write_bin(ofilegrvi,vi,infile)
                    

    return vi
