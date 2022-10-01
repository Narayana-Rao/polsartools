from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d,load_C2,eig22


def dpdesc(cpol,xpol,window_size=1,write_flag=None):

    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

    c11 = conv2d(read_bin(cpol),kernel)
    c22 = conv2d(read_bin(xpol),kernel)

    q = c12/c11
    q[q>=1]=1
    mc = (1-q)/(1+q)
    p1 = 1/(1+q)
    p2 = q/(1+q)
    Hc = -1*(p1*np.log2(p1)+p2*np.log2(p2))
    thetac = np.arctan(((1-q)**2)/(1-q+q**2)) * (180/np.pi)

    if write_flag:
        infile = copol
        ofile = os.path.join(os.path.dirname(copol),'mc.bin')
        write_bin(ofile,mc,infile)
        ofile = os.path.join(os.path.dirname(copol),'Hc.bin')
        write_bin(ofile,Hc,infile)
        ofile = os.path.join(os.path.dirname(copol),'Thetac.bin')
        write_bin(ofile,thetac,infile)
                    
    return np.real(mc),np.real(Hc),np.real(thetac)
