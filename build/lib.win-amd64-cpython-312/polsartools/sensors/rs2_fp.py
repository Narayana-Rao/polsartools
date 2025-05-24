import numpy as np
from osgeo import gdal
import os
import xml.etree.ElementTree as ET
from polsartools.utils.utils import time_it
from polsartools.utils.io_utils import write_T3, write_C3,mlook

def read_rs2_tif(file):
    ds = gdal.Open(file)
    band1 = ds.GetRasterBand(1).ReadAsArray()
    band2 = ds.GetRasterBand(2).ReadAsArray()
    ds=None
    return np.dstack((band1,band2))

def write_s2_bin(file,wdata):
    [cols, rows] = wdata.shape
    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_CFloat32)
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    outdata.FlushCache()

@time_it
def rs2_fp(inFolder,matrix='S2',type='sigma0',azlks=8,rglks=2):
    
    if type == 'sigma0':
        tree = ET.parse(os.path.join(inFolder,"lutSigma.xml"))
        root = tree.getroot()
        lut = root.find('gains').text.strip()
        lut = np.fromstring(lut, sep=' ')
    elif type == 'gamma0':
        tree = ET.parse(os.path.join(inFolder,"lutGamma.xml"))
        root = tree.getroot()
        lut = root.find('gains').text.strip()
        lut = np.fromstring(lut, sep=' ')
    elif type=='beta0':
        tree = ET.parse(os.path.join(inFolder,"lutBeta.xml"))
        root = tree.getroot()
        lut = root.find('gains').text.strip()
        lut = np.fromstring(lut, sep=' ')
    else:
        raise ValueError(f'Unknown type {type} \n Available types: sigma0,gamma0,beta0')

    if matrix == 'S2':

        out_dir = os.path.join(inFolder,"S2")
        os.makedirs(out_dir,exist_ok=True)

        print("Considering S12 = S21")
        inFile = os.path.join(inFolder,"imagery_HH.tif")
        data = read_rs2_tif(inFile)
        out_file = os.path.join(out_dir,'s11.bin')
        write_s2_bin(out_file,data[:,:,0]/lut+1j*(data[:,:,1]/lut))
        print("Saved file "+out_file)
        
        rows,cols,_ = data.shape
        
        inFile = os.path.join(inFolder,"imagery_HV.tif")
        data_xy = read_rs2_tif(inFile)

        inFile = os.path.join(inFolder,"imagery_VH.tif")
        data_yx = read_rs2_tif(inFile)

        data = (data_xy+data_yx)*0.5
        del data_xy,data_yx

        out_file = os.path.join(out_dir,'s12.bin')
        
        write_s2_bin(out_file,data[:,:,0]/lut+1j*(data[:,:,1]/lut))
        print("Saved file "+out_file)
        out_file = os.path.join(out_dir,'s21.bin')
        write_s2_bin(out_file,data[:,:,0]/lut+1j*(data[:,:,1]/lut))
        print("Saved file "+out_file)

        inFile = os.path.join(inFolder,"imagery_VV.tif")
        data = read_rs2_tif(inFile)
        out_file = os.path.join(out_dir,'s22.bin')
        write_s2_bin(out_file,data[:,:,0]/lut+1j*(data[:,:,1]/lut))
        print("Saved file "+out_file)
        
        file = open(out_dir +'/config.txt',"w+")
        file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
        file.close() 
        
    elif matrix == 'T3':
        print("Considering S12 = S21")
        # Kp- 3-D Pauli feature vector
        # Kp = (1/np.sqrt(2))*np.array([S2[0,0]+S2[1,1], S2[0,0]-S2[1,1], S2[1,0]])
        # Kp = (1/np.sqrt(2))*np.array([S2[0,0]+S2[1,1], S2[0,0]-S2[1,1], S2[0,1]])
        
        inFile = os.path.join(inFolder,"imagery_HH.tif")
        data = read_rs2_tif(inFile)
        s11 = data[:,:,0]/lut+1j*(data[:,:,1]/lut)

        inFile = os.path.join(inFolder,"imagery_HV.tif")
        data_xy = read_rs2_tif(inFile)

        inFile = os.path.join(inFolder,"imagery_VH.tif")
        data_yx = read_rs2_tif(inFile)
        # Symmetry assumption
        data = (data_xy+data_yx)*0.5
        del data_xy,data_yx
        s12 = data[:,:,0]/lut+1j*(data[:,:,1]/lut)

        inFile = os.path.join(inFolder,"imagery_VV.tif")
        data = read_rs2_tif(inFile)
        s22 = data[:,:,0]/lut+1j*(data[:,:,1]/lut)

        Kp = (1/np.sqrt(2))*np.array([s11+s22, s11-s22, 2*s12])

        del s11,s12,s22

        # 3x3 Pauli Coherency Matrix elements
        T11 = mlook(np.abs(Kp[0])**2,azlks,rglks)
        T22 = mlook(np.abs(Kp[1])**2,azlks,rglks)
        T33 = mlook(np.abs(Kp[2])**2,azlks,rglks)

        T12 = mlook(Kp[0]*np.conj(Kp[1]),azlks,rglks)
        T13 = mlook(Kp[0]*np.conj(Kp[2]),azlks,rglks)
        T23 = mlook(Kp[1]*np.conj(Kp[2]),azlks,rglks)

        del Kp
        T3Folder = os.path.join(inFolder,'T3')

        if not os.path.isdir(T3Folder):
            print("T3 folder does not exist. \nCreating folder {}".format(T3Folder))
            os.mkdir(T3Folder)
            
        write_T3(np.dstack([T11,T12,T13,np.conjugate(T12),T22,T23,np.conjugate(T13),np.conjugate(T23),T33]),T3Folder)
        
        
    elif matrix == 'C3':
        print("Considering S12 = S21")
        inFile = os.path.join(inFolder,"imagery_HH.tif")
        data = read_rs2_tif(inFile)
        s11 = data[:,:,0]/lut+1j*(data[:,:,1]/lut)

        inFile = os.path.join(inFolder,"imagery_HV.tif")
        data_xy = read_rs2_tif(inFile)

        inFile = os.path.join(inFolder,"imagery_VH.tif")
        data_yx = read_rs2_tif(inFile)
        # Symmetry assumption
        data = (data_xy+data_yx)*0.5
        del data_xy,data_yx
        s12 = data[:,:,0]/lut+1j*(data[:,:,1]/lut)

        inFile = os.path.join(inFolder,"imagery_VV.tif")
        data = read_rs2_tif(inFile)
        s22 = data[:,:,0]/lut+1j*(data[:,:,1]/lut)

        # Kl- 3-D Lexicographic feature vector
        Kl = np.array([s11, np.sqrt(2)*s12, s22])
        del s11,s12,s22

        # 3x3 COVARIANCE Matrix elements

        C11 = mlook(np.abs(Kl[0])**2,azlks,rglks)
        C22 = mlook(np.abs(Kl[1])**2,azlks,rglks)
        C33 = mlook(np.abs(Kl[2])**2,azlks,rglks)

        C12 = mlook(Kl[0]*np.conj(Kl[1]),azlks,rglks)
        C13 = mlook(Kl[0]*np.conj(Kl[2]),azlks,rglks)
        C23 = mlook(Kl[1]*np.conj(Kl[2]),azlks,rglks)

        C3Folder = os.path.join(inFolder,'C3')

        if not os.path.isdir(C3Folder):
            print("C3 folder does not exist. \nCreating folder {}".format(C3Folder))
            os.mkdir(C3Folder)
        
        write_C3(np.dstack([C11,C12,C13,np.conjugate(C12),C22,C23,np.conjugate(C13),np.conjugate(C23),C33]),C3Folder)
        
        
    else:
        raise ValueError('Invalid matrix type. Valid types are "S2", "T3" and "C3"')