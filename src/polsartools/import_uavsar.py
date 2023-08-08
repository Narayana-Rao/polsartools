
import numpy as np
import glob
from osgeo import gdal,osr
import os
import sys

def write_bin_uav(file,wdata,lat,lon,dx,dy):
    
    [cols, rows] = wdata.shape

    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform([lon,dx,0,lat,0,dy])
    outdata.SetProjection("EPSG:4326")
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    # outdata.GetRasterBand(1).SetNoDataValue(np.NaN)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!!

def uavsar_grd_c3(inFolder):
    annFile = open(glob.glob(inFolder+'/*.ann')[0], 'r')
    for line in annFile:
        if "grd_mag.set_rows" in line:
            rows = int(line.split('=')[1].split(';')[0])
        if "grd_mag.set_cols" in line:
            cols = int(line.split('=')[1].split(';')[0])
        if "grd_mag.row_addr" in line:
            lat = float(line.split('=')[1].split(';')[0])
        if "grd_mag.col_addr" in line:
            lon = float(line.split('=')[1].split(';')[0])
        if "grd_mag.row_mult" in line and "(deg/pixel)" in line:
            dy = float(line.split('=')[1].split(';')[0])
        if "grd_mag.col_mult" in line and "(deg/pixel)" in line:
            dx = float(line.split('=')[1].split(';')[0])

    outFolder = inFolder+'/C3'
    if not os.path.isdir(outFolder):
        os.mkdir(outFolder)

    hhhh = np.fromfile(glob.glob(inFolder+'/*HHHH*.grd')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C11.bin',hhhh,lat,lon,dx,dy)
    del hhhh
    vvvv = np.fromfile(glob.glob(inFolder+'/*VVVV*.grd')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C33.bin',vvvv,lat,lon,dx,dy)
    del vvvv
    hvhv = np.fromfile(glob.glob(inFolder+'/*HVHV*.grd')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C22.bin',2*hvhv,lat,lon,dx,dy)
    del hvhv
    hhhv = np.fromfile(glob.glob(inFolder+'/*HHHV*.grd')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C12_real.bin',np.real(np.sqrt(2)*hhhv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C12_imag.bin',np.imag(np.sqrt(2)*hhhv),lat,lon,dx,dy)
    del hhhv
    hhvv = np.fromfile(glob.glob(inFolder+'/*HHVV*.grd')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C13_real.bin',np.real(hhvv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C13_imag.bin',np.imag(hhvv),lat,lon,dx,dy)
    del hhvv
    hvvv = np.fromfile(glob.glob(inFolder+'/*HVVV*.grd')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C23_real.bin',np.real(np.sqrt(2)*hvvv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C23_imag.bin',np.imag(np.sqrt(2)*hvvv),lat,lon,dx,dy)
    del hvvv

    file = open(outFolder +'/config.txt',"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()   

def uavsar_mlc_c3(inFolder):
    annFile = open(glob.glob(inFolder+'/*.ann')[0], 'r')
    for line in annFile:
        if "mlc_mag.set_rows" in line:
            rows = int(line.split('=')[1].split(';')[0])
        if "mlc_mag.set_cols" in line:
            cols = int(line.split('=')[1].split(';')[0])
        if "set_plon" in line  and "(deg)" in line:
            lon = float(line.split('=')[1].split(';')[0])
        if "set_plat" in line  and "(deg)" in line:
            lat = float(line.split('=')[1].split(';')[0])
        if "mlc_mag.row_mult" in line and "(m/pixel) " in line:
            dy = float(line.split('=')[1].split(';')[0])*0.00001
        if "mlc_mag.col_mult" in line and "(m/pixel) " in line:
            dx = float(line.split('=')[1].split(';')[0])*0.00001

    outFolder = inFolder+'/C3'
    if not os.path.isdir(outFolder):
        os.mkdir(outFolder)

        
    hhhh = np.fromfile(glob.glob(inFolder+'/*HHHH*.mlc')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C11.bin',hhhh,lat,lon,dx,dy)

    del hhhh
    vvvv = np.fromfile(glob.glob(inFolder+'/*VVVV*.mlc')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C33.bin',vvvv,lat,lon,dx,dy)
    del vvvv
    hvhv = np.fromfile(glob.glob(inFolder+'/*HVHV*.mlc')[0], dtype='<f',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C22.bin',2*hvhv,lat,lon,dx,dy)
    del hvhv
    hhhv = np.fromfile(glob.glob(inFolder+'/*HHHV*.mlc')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C12_real.bin',np.real(np.sqrt(2)*hhhv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C12_imag.bin',np.imag(np.sqrt(2)*hhhv),lat,lon,dx,dy)
    del hhhv
    hhvv = np.fromfile(glob.glob(inFolder+'/*HHVV*.mlc')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C13_real.bin',np.real(hhvv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C13_imag.bin',np.imag(hhvv),lat,lon,dx,dy)
    del hhvv
    hvvv = np.fromfile(glob.glob(inFolder+'/*HVVV*.mlc')[0], dtype='<F',).reshape(rows,cols)
    write_bin_uav(outFolder+'/C23_real.bin',np.real(np.sqrt(2)*hvvv),lat,lon,dx,dy)
    write_bin_uav(outFolder+'/C23_imag.bin',np.imag(np.sqrt(2)*hvvv),lat,lon,dx,dy)
    del hvvv

    file = open(outFolder +'/config.txt',"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()   
