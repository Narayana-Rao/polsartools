from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d,load_C2,eig22


def dprvi(cpol,xpol,window_size=1,write_flag=None):

    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

    c11 = conv2d(read_bin(cpol),kernel)
    c22 = conv2d(read_bin(xpol),kernel)

    q = c12/c11
    q[q>=1]=1
    DpRVIc = (q*(q+3)/((q+1)**2))


    if write_flag:
        infile = copol
        ofile = os.path.join(os.path.dirname(copol),'DpRVIc.bin')
        write_bin(ofile,DpRVIc,infile)
                    
    return np.real(DpRVIc)
