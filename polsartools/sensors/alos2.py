
import numpy as np
from osgeo import gdal
import os, glob
import matplotlib.pyplot as plt
from skimage.util.shape import view_as_blocks
from polsartools.utils.utils import time_it
def read_bin(file):
    ds = gdal.Open(file)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr
def mlook(data,az,rg): 
    temp = data[0:data.shape[0]-data.shape[0]%az,0:data.shape[1]-data.shape[1]%rg]
    blocks = view_as_blocks(temp, block_shape=(az, rg))
    flatten = blocks.reshape(blocks.shape[0], blocks.shape[1], -1)
    mean = np.nanmean(flatten, axis=2)
    return mean

def read_a2(file):
    
    fp = open(file,mode='rb')
    fp.seek(232)
    ch = int(fp.read(4))
    # print(ch)
    
    fp.seek(236)
    nline = int(fp.read(8))
    # print(nline)
    fp.seek(248)
    npixel = int(fp.read(8))
    # print(npixel)
    
    nrec = 544 + npixel*8
    # print(nrec)
    fp.seek(720)
    data = np.frombuffer(fp.read(int(nrec * nline)), dtype='>f4')
    data = np.array(data).reshape(-1,int(nrec/4)) 
    # print(np.shape(data))
    
    data = data[:,int(544/4):int(nrec/4)] 
    slc = data[:,::2] + 1j*data[:,1::2]
    # print(np.shape(slc))
    del data
    
    return slc

def write_a2_bin(file,wdata):
    
    # ds = gdal.Open(refData)
    [cols, rows] = wdata.shape

    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    # outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    # outdata.SetProjection(ds.GetProjection())##sets same projection as input
    
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    # outdata.GetRasterBand(1).SetNoDataValue(0)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!!

@time_it    
def alos2_fbd_l11(inFolder,azlks=3,rglks=2,calfac_dB=-83):
    
    hh_file = list(glob.glob(os.path.join(inFolder,'IMG-HH-*-FBDR1.1__A')) + \
        glob.glob(os.path.join(inFolder, 'IMG-HH-*-FBDR1.1__D')))[0]

    hv_file = list(glob.glob(os.path.join(inFolder,'IMG-HV-*-FBDR1.1__A')) + \
        glob.glob(os.path.join(inFolder, 'IMG-HV-*-FBDR1.1__D')))[0]

    calfac_linear = np.sqrt(10 ** ((calfac_dB - 32) / 10))
    
    S11 = read_a2(hh_file).astype(np.complex64)*calfac_linear 
    S12 = read_a2(hv_file).astype(np.complex64)*calfac_linear 
    
    
    C11 = mlook(np.abs(S11)**2,azlks,rglks).astype(np.float32)
    C22 = mlook(np.abs(S12)**2,azlks,rglks).astype(np.float32)
    
    C12 = mlook(S11*np.conjugate(S12),azlks,rglks).astype(np.complex64)
    rows,cols = C11.shape
    
    inFolder = os.path.dirname(hh_file)   
    C2Folder = os.path.join(inFolder,os.path.basename(hh_file).split('-HH-')[1].split('1.1')[0],'C2')

    del S11,S12
    
    os.makedirs(C2Folder,exist_ok=True)
    
    write_a2_bin( os.path.join(C2Folder,'C11.bin'), C11)
    print(f"Saved file {C2Folder}/C11.bin")
    del C11  
    write_a2_bin( os.path.join(C2Folder,'C22.bin'), C22)
    print(f"Saved file {C2Folder}/C22.bin")
    del C22
     
    write_a2_bin( os.path.join(C2Folder,'C12_real.bin'), np.real(C12))
    print(f"Saved file {C2Folder}/C12_real.bin")
    write_a2_bin( os.path.join(C2Folder,'C12_imag.bin'), np.imag(C12))
    print(f"Saved file {C2Folder}/C12_imag.bin")
    del C12
    
    
    file = open(C2Folder +'/config.txt',"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\npp1'%(rows,cols))
    file.close()  
