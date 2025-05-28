import numpy as np
from osgeo import gdal
import os
import xml.etree.ElementTree as ET
from polsartools.utils.utils import time_it
from polsartools.utils.io_utils import write_T3, write_C3,mlook
import h5py

def write_s2_bin(file,wdata):
    [cols, rows] = wdata.shape
    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_CFloat32)
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    outdata.FlushCache()

@time_it
def isro_asar(inFile,matrix='T3',azlks=20,rglks=10,calibration_constant = 42):
    
    
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
        
        # xCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinateSpacing'])
        # yCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinateSpacing'])
        # xCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinates'])
        # yCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinates'])
        # projection = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/radarGrid/projection'])

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
        
        
    elif matrix == 'C3':
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
        
        
    else:
        raise ValueError('Invalid matrix type. Valid types are "S2", "T3" and "C3"')