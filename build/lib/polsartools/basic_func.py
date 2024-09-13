from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

def conv2d(a, f):
    filt = np.zeros(a.shape)
    wspad = int(f.shape[0]/2)
    s = f.shape + tuple(np.subtract(a.shape, f.shape) + 1)
    strd = np.lib.stride_tricks.as_strided
    subM = strd(a, shape = s, strides = a.strides * 2)
    filt_data = np.einsum('ij,ijkl->kl', f, subM)
    filt[wspad:wspad+filt_data.shape[0],wspad:wspad+filt_data.shape[1]] = filt_data
    return filt

def read_bin(file):
    ds = gdal.Open(file)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr

def write_bin(file,wdata,refData):        
    ds = gdal.Open(refData)
    [cols, rows] = wdata.shape   
    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform(ds.GetGeoTransform())
    outdata.SetProjection(ds.GetProjection())   
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    outdata.FlushCache() 

def load_C2(folder):

    C11 = read_bin(os.path.join(folder,"C11.bin"))
    C22 = read_bin(os.path.join(folder,"C22.bin"))

    C12_i = read_bin(os.path.join(folder,'C12_imag.bin'))
    C12_r = read_bin(os.path.join(folder,'C12_real.bin'))

    C12 = C12_r + 1j*C12_i

    return np.dstack((C11,C12,np.conj(C12),C22))

def load_C3(folder):
    
    C11 = read_bin(os.path.join(folder,"C11.bin"))
    C22 = read_bin(os.path.join(folder,"C22.bin"))
    C33 = read_bin(os.path.join(folder,"C33.bin"))

    C12_i = read_bin(os.path.join(folder,'C12_imag.bin'))
    C12_r = read_bin(os.path.join(folder,'C12_real.bin'))
    C13_i = read_bin(os.path.join(folder,'C13_imag.bin'))
    C13_r = read_bin(os.path.join(folder,'C13_real.bin'))
    C23_i = read_bin(os.path.join(folder,'C23_imag.bin'))
    C23_r = read_bin(os.path.join(folder,'C23_real.bin'))
        
    C12 = C12_r + 1j*C12_i
    C13 = C13_r + 1j*C13_i
    C23 = C23_r + 1j*C23_i
    
    return np.dstack((C11,C12,C13,np.conj(C12),C22,C23,np.conj(C13),np.conj(C23),C33))


def load_T3(folder):
    
    T11 = read_bin(os.path.join(folder,"T11.bin"))
    T22 = read_bin(os.path.join(folder,"T22.bin"))
    T33 = read_bin(os.path.join(folder,"T33.bin"))

    T12_i = read_bin(os.path.join(folder,'T12_imag.bin'))
    T12_r = read_bin(os.path.join(folder,'T12_real.bin'))
    T13_i = read_bin(os.path.join(folder,'T13_imag.bin'))
    T13_r = read_bin(os.path.join(folder,'T13_real.bin'))
    T23_i = read_bin(os.path.join(folder,'T23_imag.bin'))
    T23_r = read_bin(os.path.join(folder,'T23_real.bin'))
        
    T12 = T12_r + 1j*T12_i
    T13 = T13_r + 1j*T13_i
    T23 = T23_r + 1j*T23_i
    
    return np.dstack((T11,T12,T13,np.conj(T12),T22,T23,np.conj(T13),np.conj(T23),T33))

def load_T2(folder):

    T11 = read_bin(os.path.join(folder,"T11.bin"))
    T22 = read_bin(os.path.join(folder,"T22.bin"))

    T12_i = read_bin(os.path.join(folder,'T12_imag.bin'))
    T12_r = read_bin(os.path.join(folder,'T12_real.bin'))
    T12 = T12_r + 1j*T12_i

    return np.dstack((T11,T12,np.conj(T12),T22))

def eig22(c2):
    c11 = c2[:,:,0].flatten()
    c12 = c2[:,:,1].flatten()
    c21 = c2[:,:,2].flatten()
    c22 = c2[:,:,3].flatten()
    trace = -(c11+c22)
    det = c11*c22-c12*c21
    # const= 1
    sqdiscr = np.sqrt(trace*trace - 4*det);
    lambda1 = -(trace + sqdiscr)*0.5;
    lambda2 = -(trace - sqdiscr)*0.5;
    
    return lambda1,lambda2