import numpy as np
from osgeo import gdal
import os,h5py
import xml.etree.ElementTree as ET
from polsartools.utils.utils import time_it
from polsartools.utils.io_utils import mlook, write_T3, write_C3
from polsartools.utils.geo_utils import geocode_grid, intp_grid, update_vrt, write_latlon


def write_s2_bin(file,wdata):
    [cols, rows] = wdata.shape
    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_CFloat32)
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    outdata.FlushCache()



def asar_geo_T3(inFile,azlks,rglks,cc_linear):
    print("Considering S12 = S21")
           
    try:
        h5File = h5py.File(inFile,"r")
    except:
        raise('Invalid .h5 file !!')
    
    
    h5File = h5py.File(inFile,"r")
    if '/science/LSAR' in h5File:
        freq_band = 'L'            
        print("Detected L-band data ")
    elif '/science/SSAR' in h5File:
        freq_band = 'S'
        print("Detected S-band data")
    else:
        print("Neither LSAR nor SSAR data found in the file.")
        h5File.close()
        return
    
    S11 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HH'])/cc_linear
    S12 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HV'])/cc_linear
    S21 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VH'])/cc_linear
    S22 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VV'])/cc_linear
    
    # xCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinateSpacing'])
    # yCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinateSpacing'])
    # xCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinates'])
    # yCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinates'])
    # projection = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/radarGrid/projection'])

    sceneCenterAlongTrackSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/sceneCenterAlongTrackSpacing' ])
    sceneCenterGroundRangeSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/sceneCenterGroundRangeSpacing'])
    slantRange = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/slantRange'])
    slantRangeSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/slantRangeSpacing'])
    coordinateX = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/geolocationGrid/coordinateX' ])[1,:,:] # at 0 height
    coordinateY = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/geolocationGrid/coordinateY'])[1,:,:] # at 0 height
    epsg = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/geolocationGrid/epsg'])

    multi_look_factor = slantRangeSpacing/sceneCenterAlongTrackSpacing


    h5File.close()


    mlrows,mlcols = S11.shape[0]//azlks,S11.shape[1]//rglks


    inFolder = os.path.dirname(inFile)   
    if not inFolder:
        inFolder = "."
    T3Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T3')
    # os.makedirs(T3Folder,exist_ok=True)

    if not os.path.isdir(T3Folder):
        print("T3 folder does not exist. \nCreating folder {}".format(T3Folder))
        # os.mkdir(T3Folder)
        os.makedirs(T3Folder,exist_ok=True)
    
    os.chdir(T3Folder)
    ######################################################
    lat_interp = intp_grid(coordinateY,(mlrows,mlcols))
    write_latlon('lat.bin', lat_interp)
    long_interp = intp_grid(coordinateX,(mlrows,mlcols))
    write_latlon('lon.bin', long_interp)

    # print('lat lon saved')
    gdal.Translate('lat.vrt', 'lat.bin')
    gdal.Translate('lon.vrt', 'lon.bin')

    minX, maxX = np.nanmin(long_interp),np.nanmax(long_interp)
    minY,maxY = np.nanmin(lat_interp),np.nanmax(lat_interp)

    x_res_m = sceneCenterAlongTrackSpacing*azlks  # in meters
    y_res_m = slantRangeSpacing*rglks
    mean_lat = np.nanmean(lat_interp) 
    x_res_deg = x_res_m / (111320 * np.cos(np.deg2rad(mean_lat)))
    y_res_deg = y_res_m / 111320.0

    bbox = [minX, minY, maxX, maxY]

    #######################################################

    # 3x1 Pauli Coherency Matrix
    Kp = (1/np.sqrt(2))*np.array([S11+S22, S11-S22, S12+S21])

    del S11,S12,S21,S22


    def write_element(data,outname):
        write_latlon('data.bin', data)
        gdal.Translate('data.vrt', 'data.bin')
        update_vrt('data.vrt')
        input_vrt = "data.vrt"
        geocode_grid(input_vrt, outname, 'bin', x_res_deg, y_res_deg, bbox,resampleAlg="near")
        os.remove('data.vrt')
        os.remove('data.bin')
        os.remove('data.hdr')
        
    # 3x3 Pauli Coherency Matrix elements
    T11 = mlook(np.abs(Kp[0])**2,azlks,rglks).astype(np.float32)
    write_element(T11,'T11.bin')
    print(f"Saved file {os.path.join(T3Folder,'T11.bin')}")
    del T11
    T22 = mlook(np.abs(Kp[1])**2,azlks,rglks).astype(np.float32)
    write_element(T22,'T22.bin')
    print(f"Saved file {os.path.join(T3Folder,'T22.bin')}")
    del T22
    T33 = mlook(np.abs(Kp[2])**2,azlks,rglks).astype(np.float32)
    write_element(T33,'T33.bin')
    print(f"Saved file {os.path.join(T3Folder,'T33.bin')}")
    del T33
    T12 = mlook(Kp[0]*np.conj(Kp[1]),azlks,rglks).astype(np.complex64)
    write_element(np.real(T12),'T12_real.bin') 
    print(f"Saved file {os.path.join(T3Folder,'T12_real.bin')}")
    write_element(np.imag(T12),'T12_imag.bin')
    print(f"Saved file {os.path.join(T3Folder,'T12_imag.bin')}")
    del T12
    
    T13 = mlook(Kp[0]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)
    write_element(np.real(T13),'T13_real.bin')
    print(f"Saved file {os.path.join(T3Folder,'T13_real.bin')}")
    write_element(np.imag(T13),'T13_imag.bin')
    print(f"Saved file {os.path.join(T3Folder,'T13_imag.bin')}")
    del T13
    
    T23 = mlook(Kp[1]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)
    write_element(np.real(T23),'T23_real.bin')
    print(f"Saved file {os.path.join(T3Folder,'T23_real.bin')}")
    write_element(np.imag(T23),'T23_imag.bin')
    print(f"Saved file {os.path.join(T3Folder,'T23_imag.bin')}")

    del Kp,T23
    
    os.remove('lat.bin')
    os.remove('lat.vrt')
    os.remove('lat.hdr')
    os.remove('lon.bin')
    os.remove('lon.vrt')
    os.remove('lon.hdr')
    
def asar_T3(inFile,azlks,rglks,cc_linear):
    print("Considering S12 = S21")
           
    try:
        h5File = h5py.File(inFile,"r")
    except:
        raise('Invalid .h5 file !!')
    
    
    h5File = h5py.File(inFile,"r")
    if '/science/LSAR' in h5File:
        freq_band = 'L'            
        print("Detected L-band data ")
    elif '/science/SSAR' in h5File:
        freq_band = 'S'
        print("Detected S-band data")
    else:
        print("Neither LSAR nor SSAR data found in the file.")
        h5File.close()
        return
    
    S11 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HH'])/cc_linear
    S12 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HV'])/cc_linear
    S21 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VH'])/cc_linear
    S22 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VV'])/cc_linear

    h5File.close()

    # 3x1 Pauli Coherency Matrix
    Kp = (1/np.sqrt(2))*np.array([S11+S22, S11-S22, S12+S21])

    del S11,S12,S21,S22

    # 3x3 Pauli Coherency Matrix elements
    T11 = mlook(np.abs(Kp[0])**2,azlks,rglks).astype(np.float32)
    T22 = mlook(np.abs(Kp[1])**2,azlks,rglks).astype(np.float32)
    T33 = mlook(np.abs(Kp[2])**2,azlks,rglks).astype(np.float32)

    T12 = mlook(Kp[0]*np.conj(Kp[1]),azlks,rglks).astype(np.complex64)
    T13 = mlook(Kp[0]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)
    T23 = mlook(Kp[1]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)

    del Kp
    
    
    inFolder = os.path.dirname(inFile)   
    if not inFolder:
        inFolder = "."
    T3Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T3')
    os.makedirs(T3Folder,exist_ok=True)

    if not os.path.isdir(T3Folder):
        print("T3 folder does not exist. \nCreating folder {}".format(T3Folder))
        os.mkdir(T3Folder)
        
    # write_T3(np.dstack([T11,T12,T13,np.conjugate(T12),T22,T23,np.conjugate(T13),np.conjugate(T23),T33]),T3Folder)
    write_T3([np.real(T11),np.real(T12),np.imag(T12),np.real(T13),np.imag(T13),
                np.real(T22),np.real(T23),np.imag(T23),
                np.real(T33)],T3Folder)

def asar_geo_C3(inFile,azlks,rglks,cc_linear):
    print("Considering S12 = S21")
    try:
        h5File = h5py.File(inFile,"r")
    except:
        raise('Invalid .h5 file !!')
    
    
    h5File = h5py.File(inFile,"r")
    if '/science/LSAR' in h5File:
        freq_band = 'L'            
        print("Detected L-band data ")
    elif '/science/SSAR' in h5File:
        freq_band = 'S'
        print(" Detected S-band data")
    else:
        print("Neither LSAR nor SSAR data found in the file.")
        h5File.close()
        return
    
    S11 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HH'])/cc_linear
    S12 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HV'])/cc_linear
    S21 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VH'])/cc_linear
    S22 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VV'])/cc_linear
    
    sceneCenterAlongTrackSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/sceneCenterAlongTrackSpacing' ])
    sceneCenterGroundRangeSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/sceneCenterGroundRangeSpacing'])
    slantRange = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/slantRange'])
    slantRangeSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/slantRangeSpacing'])
    coordinateX = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/geolocationGrid/coordinateX' ])[1,:,:] # at 0 height
    coordinateY = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/geolocationGrid/coordinateY'])[1,:,:] # at 0 height
    epsg = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/geolocationGrid/epsg'])

    multi_look_factor = slantRangeSpacing/sceneCenterAlongTrackSpacing


    h5File.close()


    mlrows,mlcols = S11.shape[0]//azlks,S11.shape[1]//rglks


    inFolder = os.path.dirname(inFile)   
    if not inFolder:
        inFolder = "."
    C3Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C3')
    
    
    if not os.path.isdir(C3Folder):
        print("C3 folder does not exist. \nCreating folder {}".format(C3Folder))
        os.makedirs(C3Folder,exist_ok=True)
    
    os.chdir(C3Folder)
    ######################################################
    lat_interp = intp_grid(coordinateY,(mlrows,mlcols))
    write_latlon('lat.bin', lat_interp)
    long_interp = intp_grid(coordinateX,(mlrows,mlcols))
    write_latlon('lon.bin', long_interp)

    # print('lat lon saved')
    gdal.Translate('lat.vrt', 'lat.bin')
    gdal.Translate('lon.vrt', 'lon.bin')

    minX, maxX = np.nanmin(long_interp),np.nanmax(long_interp)
    minY,maxY = np.nanmin(lat_interp),np.nanmax(lat_interp)

    x_res_m = sceneCenterAlongTrackSpacing*azlks  # in meters
    y_res_m = slantRangeSpacing*rglks
    mean_lat = np.nanmean(lat_interp) 
    x_res_deg = x_res_m / (111320 * np.cos(np.deg2rad(mean_lat)))
    y_res_deg = y_res_m / 111320.0

    bbox = [minX, minY, maxX, maxY]

    #######################################################


    # Kl- 3-D Lexicographic feature vector
    Kl = np.array([S11, np.sqrt(2)*0.5*(S12+S21), S22])
    del S11, S12, S21, S22
    def write_element(data,outname):
        write_latlon('data.bin', data)
        gdal.Translate('data.vrt', 'data.bin')
        update_vrt('data.vrt')
        input_vrt = "data.vrt"
        geocode_grid(input_vrt, outname, 'bin', x_res_deg, y_res_deg, bbox,resampleAlg="near")
        os.remove('data.vrt')
        os.remove('data.bin')
        os.remove('data.hdr')
        
        
    # 3x3 COVARIANCE Matrix elements

    C11 = mlook(np.abs(Kl[0])**2,azlks,rglks).astype(np.float32)
    write_element(C11,'C11.bin')
    print(f"Saved file {os.path.join(C3Folder,'C11.bin')}")
    del C11
    C22 = mlook(np.abs(Kl[1])**2,azlks,rglks).astype(np.float32)
    write_element(C22,'C22.bin')
    print(f"Saved file {os.path.join(C3Folder,'C22.bin')}")
    del C22
    C33 = mlook(np.abs(Kl[2])**2,azlks,rglks).astype(np.float32)
    write_element(C33,'C33.bin')
    print(f"Saved file {os.path.join(C3Folder,'C33.bin')}")
    del C33 

    C12 = mlook(Kl[0]*np.conj(Kl[1]),azlks,rglks).astype(np.complex64)
    write_element(np.real(C12),'C12_real.bin')
    print(f"Saved file {os.path.join(C3Folder,'C12_real.bin')}")
    write_element(np.imag(C12),'C12_imag.bin')
    print(f"Saved file {os.path.join(C3Folder,'C12_imag.bin')}")
    del C12
    C13 = mlook(Kl[0]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
    write_element(np.real(C13),'C13_real.bin')
    print(f"Saved file {os.path.join(C3Folder,'C13_real.bin')}")
    write_element(np.imag(C13),'C13_imag.bin')
    print(f"Saved file {os.path.join(C3Folder,'C13_imag.bin')}")
    del C13
    C23 = mlook(Kl[1]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
    write_element(np.real(C23),'C23_real.bin')
    print(f"Saved file {os.path.join(C3Folder,'C23_real.bin')}")
    write_element(np.imag(C23),'C23_imag.bin')
    print(f"Saved file {os.path.join(C3Folder,'C23_imag.bin')}")
    del C23, Kl
    os.remove('lat.bin')
    os.remove('lat.vrt')
    os.remove('lat.hdr')
    os.remove('lon.bin')
    os.remove('lon.vrt')
    os.remove('lon.hdr')


def asar_C3(inFile,azlks,rglks,cc_linear):    
    
    print("Considering S12 = S21")
    try:
        h5File = h5py.File(inFile,"r")
    except:
        raise('Invalid .h5 file !!')
    
    
    h5File = h5py.File(inFile,"r")
    if '/science/LSAR' in h5File:
        freq_band = 'L'            
        print("Detected L-band data ")
    elif '/science/SSAR' in h5File:
        freq_band = 'S'
        print(" Detected S-band data")
    else:
        print("Neither LSAR nor SSAR data found in the file.")
        h5File.close()
        return
    
    S11 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HH'])/cc_linear
    S12 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HV'])/cc_linear
    S21 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VH'])/cc_linear
    S22 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VV'])/cc_linear
    
    # xCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinateSpacing'])
    # yCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinateSpacing'])
    # xCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinates'])
    # yCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinates'])
    # projection = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/radarGrid/projection'])

    h5File.close()


    # Kl- 3-D Lexicographic feature vector
    Kl = np.array([S11, np.sqrt(2)*0.5*(S12+S21), S22])
    del S11, S12, S21, S22

    # 3x3 COVARIANCE Matrix elements

    C11 = mlook(np.abs(Kl[0])**2,azlks,rglks).astype(np.float32)
    C22 = mlook(np.abs(Kl[1])**2,azlks,rglks).astype(np.float32)
    C33 = mlook(np.abs(Kl[2])**2,azlks,rglks).astype(np.float32)

    C12 = mlook(Kl[0]*np.conj(Kl[1]),azlks,rglks).astype(np.complex64)
    C13 = mlook(Kl[0]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
    C23 = mlook(Kl[1]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)


    inFolder = os.path.dirname(inFile)   
    if not inFolder:
        inFolder = "."
    C3Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C3')
    os.makedirs(C3Folder,exist_ok=True)
    
    if not os.path.isdir(C3Folder):
        print("C3 folder does not exist. \nCreating folder {}".format(C3Folder))
        os.mkdir(C3Folder)
    
    # write_C3(np.dstack([C11,C12,C13,np.conjugate(C12),C22,C23,np.conjugate(C13),np.conjugate(C23),C33]),C3Folder)
    write_C3([np.real(C11),np.real(C12),np.imag(C12),np.real(C13),np.imag(C13),
                np.real(C22),np.real(C23),np.imag(C23),
                np.real(C33)],C3Folder)
@time_it
def isro_asar(inFile,matrix='T3',azlks=8,rglks=6,geocode_flag=False,calibration_constant = 42):
    
    cc_linear = np.sqrt(10**(calibration_constant/10))

    if matrix == 'S2':


        print("Considering S12 = S21")
            
        try:
            h5File = h5py.File(inFile,"r")
        except:
            raise('Invalid .h5 file !!')
        
        h5File = h5py.File(inFile,"r")
        if '/science/LSAR' in h5File:
            freq_band = 'L'            
            print("Detected L-band data ")
        elif '/science/SSAR' in h5File:
            freq_band = 'S'
            print(" Detected S-band data")
        else:
            print("Neither LSAR nor SSAR data found in the file.")
            h5File.close()
            return
        
        S11 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HH'])/cc_linear
        S12 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HV'])/cc_linear
        S21 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VH'])/cc_linear
        S22 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VV'])/cc_linear

        h5File.close()
        
        
        rows,cols,_ = S11.shape
        inFolder = os.path.dirname(inFile)   
        if not inFolder:
            inFolder = "."
        out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'S2')
        os.makedirs(out_dir,exist_ok=True)
       
        
        out_file = os.path.join(out_dir,'s11.bin')
        write_s2_bin(out_file,S11)
        print("Saved file "+out_file)
        
        out_file = os.path.join(out_dir,'s12.bin')
        write_s2_bin(out_file,(S12+S21)/2)
        print("Saved file "+out_file)
        
        out_file = os.path.join(out_dir,'s21.bin')
        write_s2_bin(out_file,(S12+S21)/2)
        print("Saved file "+out_file)
        
        out_file = os.path.join(out_dir,'s22.bin')
        write_s2_bin(out_file,S22)
        print("Saved file "+out_file)
        
        file = open(out_dir +'/config.txt',"w+")
        file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
        file.close() 
        
    elif matrix == 'T3':
        
        if geocode_flag:
            asar_geo_T3(inFile,azlks,rglks,cc_linear)
        else:
            asar_T3(inFile,azlks,rglks,cc_linear)
        
        
        
    elif matrix == 'C3':
        if geocode_flag:
            asar_geo_C3(inFile,azlks,rglks,cc_linear)
        else:
            asar_C3(inFile,azlks,rglks,cc_linear)
        

    else:
        raise ValueError('Invalid matrix type. Valid types are "S2", "T3" and "C3"')