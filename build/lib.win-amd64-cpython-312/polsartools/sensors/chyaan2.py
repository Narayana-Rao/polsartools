# -*- coding: utf-8 -*-
"""
Created on Sun May 18 17:40:11 2023

@author: nbhogapurapu

mismatch in corss-pol T33 etc between MIDAS and this have to cross-verify with polsarpro
other elements are fine 
calibration fine


"""


import glob,shutil,os
import numpy as np
from osgeo import gdal 

from skimage.util.shape import view_as_blocks
from polsartools.utils.utils import time_it
def read_bin(file):
    ds = gdal.Open(file,gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr
def mlook(data,az,rg):
    
    temp = data[0:data.shape[0]-data.shape[0]%az,0:data.shape[1]-data.shape[1]%rg]

    blocks = view_as_blocks(temp, block_shape=(az, rg))
    # collapse the last two dimensions in one
    flatten = blocks.reshape(blocks.shape[0], blocks.shape[1], -1)
    # resampling the image by taking either the `mean`,
    # the `max` or the `median` value of each blocks.
    mean = np.nanmean(flatten, axis=2)
    return mean
 
def write_T3(T3_stack,folder):
    
    out_file = folder +'/T11.bin'
    write_bin(out_file,np.real(T3_stack[:,:,0]))
    
    out_file = folder +'/T12_real.bin'
    write_bin(out_file,np.real(T3_stack[:,:,1]))
    out_file = folder +'/T12_imag.bin'
    write_bin(out_file,np.imag(T3_stack[:,:,1]))
    
    out_file = folder +'/T13_real.bin'
    write_bin(out_file,np.real(T3_stack[:,:,2]))
    out_file = folder +'/T13_imag.bin'
    write_bin(out_file,np.imag(T3_stack[:,:,2]))
    
    out_file = folder +'/T22.bin'
    write_bin(out_file,np.real(T3_stack[:,:,4]))
    
    out_file = folder +'/T23_real.bin'
    write_bin(out_file,np.real(T3_stack[:,:,5]))
    out_file = folder +'/T23_imag.bin'
    write_bin(out_file,np.imag(T3_stack[:,:,5]))
    
    out_file = folder +'/T33.bin'
    write_bin(out_file,np.real(T3_stack[:,:,8]))
    
    rows, cols = np.shape(T3_stack[:,:,8])
    file = folder +'/config.txt'
    file = open(file,"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()   

def write_C3(C3_stack,folder):
    
    out_file = folder +'/C11.bin'
    write_bin(out_file,np.real(C3_stack[:,:,0]))
    
    out_file = folder +'/C12_real.bin'
    write_bin(out_file,np.real(C3_stack[:,:,1]))
    out_file = folder +'/C12_imag.bin'
    write_bin(out_file,np.imag(C3_stack[:,:,1]))
    
    out_file = folder +'/C13_real.bin'
    write_bin(out_file,np.real(C3_stack[:,:,2]))
    out_file = folder +'/C13_imag.bin'
    write_bin(out_file,np.imag(C3_stack[:,:,2]))
    
    
    
    out_file = folder +'/C22.bin'
    write_bin(out_file,np.real(C3_stack[:,:,4]))
    
    out_file = folder +'/C23_real.bin'
    write_bin(out_file,np.real(C3_stack[:,:,5]))
    out_file = folder +'/C23_imag.bin'
    write_bin(out_file,np.imag(C3_stack[:,:,5]))
    
    out_file = folder +'/C33.bin'
    write_bin(out_file,np.real(C3_stack[:,:,8]))
    
    
    rows, cols = np.shape(C3_stack[:,:,8])
    file = folder +'/config.txt'
    file = open(file,"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()




def write_bin_s2(file,wdata,refData):
    
    # ds = gdal.Open(refData)
    [cols, rows] = wdata.shape

    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    # outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    # outdata.SetProjection(ds.GetProjection())##sets same projection as input
    
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    # outdata.GetRasterBand(1).SetNoDataValue(np.NaN)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!!   
    

def write_bin(file,wdata):
    
    # ds = gdal.Open(refData)
    [cols, rows] = wdata.shape

    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    # outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    # outdata.SetProjection(ds.GetProjection())##sets same projection as input
    
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    # outdata.GetRasterBand(1).SetNoDataValue(np.NaN)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!! 

@time_it
def chyaan2_fp(inFolder,matrix='T3'):
    
    S2Folder = os.path.join(inFolder,'S2')

    if not os.path.isdir(S2Folder):
        os.mkdir(S2Folder)

    ds = gdal.Open(glob.glob(inFolder+'/data/calibrated/*/*sli*_hh_*.tif')[0])
    cols = ds.RasterXSize  
    rows = ds.RasterYSize 

    # %%

    # Define band mappings and suffixes
    polarizations = ['hh', 'hv', 'vh', 'vv']
    el = ['s11','s12','s21','s22']
    suffixes = ['Real', 'Imag']


    # Function to process a band
    def convert_band(input_file, band_number, output_file):
        gdal.Translate(output_file, gdal.Open(input_file), format='ENVI', bandList=[band_number], outputType=gdal.GDT_Float32)

    # Loop over polarizations and element identifiers
    for pol, el_id in zip(polarizations, el):
        input_file = glob.glob(os.path.join(inFolder, f'data/calibrated/*/*sli*_{pol}_*.tif'))[0]
        for i, suffix in enumerate(suffixes, 1):  # Band 1 for 'Real', Band 2 for 'Imag'
            output_file = os.path.join(S2Folder, f"{el_id}_{suffix}.bin")
            convert_band(input_file, i, output_file)


    #%%
    xmlFile = glob.glob(inFolder+'/data/calibrated/*/*sli*.xml')[0]
    fxml = open(xmlFile, 'r')
    for line in fxml:
        if "output_line_spacing" in line:
            # print("output_line_spacing: ", line.split('>')[1].split('<')[0])
            ols = float(line.split('>')[1].split('<')[0])
        if "output_pixel_spacing" in line:
            # print("output_pixel_spacing: ", line.split('>')[1].split('<')[0])
            ops = float(line.split('>')[1].split('<')[0])
        if "isda:incidence_angle" in line:
            # print("incidence_angle: ", line.split('>')[1].split('<')[0])
            inc= float( line.split('>')[1].split('<')[0])
        if "isda:calibration_constant" in line:
            cc = float( line.split('>')[1].split('<')[0])
        if "isda:pulse_bandwidth" in line:
            bw = float( line.split('>')[1].split('<')[0])/1000000
    fxml.close() 
    gRange = ops/np.sin(inc*np.pi/180)
    mlf = int(np.round(gRange/ols,0))

    lines = ['output_line_spacing '+ str(ols)+'\n',
            'output_pixel_spacing '+ str(ops)+'\n',
            'ground_range '+ str(gRange)+'\n',
            'mlook_factor '+ str(mlf)+'\n',
            'incidence_angle '+ str(inc)+'\n',
            'calibration_constant '+str(cc)+'\n',
            'pulse_bandwidth '+str(bw)+'\n',
            'lines '+ str(rows)+'\n',
            'samples '+str(cols)+'\n'
            
            ]

    mlFile = S2Folder+'/ml_Info.txt'
    with open(mlFile, 'w+') as f:
        f.writelines(lines)
    f.close()


    calFactor = 1/np.sqrt(10**(cc/10))

    for file in glob.glob(S2Folder+"/*.bin"):
        write_bin_s2(file,read_bin(file)*calFactor,file)

    S2Folder = S2Folder

    HHi_ds  = gdal.Open(S2Folder+'/s11_Imag.bin')
    HHr_ds  = gdal.Open(S2Folder+'/s11_Real.bin')
    HVr_ds  = gdal.Open(S2Folder+'/s12_Real.bin')
    HVi_ds  = gdal.Open(S2Folder+'/s12_Imag.bin')
    VHr_ds  = gdal.Open(S2Folder+'/s21_Real.bin')
    VHi_ds  = gdal.Open(S2Folder+'/s21_Imag.bin')
    VVr_ds  = gdal.Open(S2Folder+'/s22_Real.bin')
    VVi_ds  = gdal.Open(S2Folder+'/s22_Imag.bin')


    S2 = np.array([
                [HHr_ds.GetRasterBand(1).ReadAsArray() + 1j*(HHi_ds.GetRasterBand(1).ReadAsArray()), HVr_ds.GetRasterBand(1).ReadAsArray()+ 1j*(HVi_ds.GetRasterBand(1).ReadAsArray()) ],
                [VHr_ds.GetRasterBand(1).ReadAsArray() + 1j*(VHi_ds.GetRasterBand(1).ReadAsArray()), VVr_ds.GetRasterBand(1).ReadAsArray()+ 1j*(VVi_ds.GetRasterBand(1).ReadAsArray()) ]
                ])
    HHi_ds = None
    HHr_ds = None
    HVr_ds = None
    HVi_ds = None
    VHr_ds = None
    VHi_ds = None
    VVr_ds = None
    VVi_ds = None
    az = mlf
    rg = 1    
    
    if matrix == 'T3':
        # Kp- 3-D Pauli feature vector
        Kp = (1/np.sqrt(2))*np.array([S2[0,0]+S2[1,1], S2[0,0]-S2[1,1], S2[1,0]])
        # Kp = (1/np.sqrt(2))*np.array([S2[0,0]+S2[1,1], S2[0,0]-S2[1,1], S2[0,1]])
        # Kp = (1/np.sqrt(2))*np.array([S2[0,0]+S2[1,1], S2[0,0]-S2[1,1], (S2[1,0]+S2[0,1])/2])

        del S2


        # 3x3 Pauli Coherency Matrix elements

        T11 = mlook(np.abs(Kp[0])**2,az,rg)
        T22 = mlook(np.abs(Kp[1])**2,az,rg)
        T33 = mlook(np.abs(Kp[2])**2,az,rg)

        T12 = mlook(Kp[0]*np.conj(Kp[1]),az,rg)
        T13 = mlook(Kp[0]*np.conj(Kp[2]),az,rg)
        T23 = mlook(Kp[1]*np.conj(Kp[2]),az,rg)

        T3Folder = os.path.join(inFolder,'T3')

        if not os.path.isdir(T3Folder):
            os.mkdir(T3Folder)
            
        write_T3(np.dstack([T11,T12,T13,np.conjugate(T12),T22,T23,np.conjugate(T13),np.conjugate(T23),T33]),T3Folder)
    elif matrix == 'C3':
        # Kl- 3-D Lexicographic feature vector
        Kl = np.array([S2[0,0], np.sqrt(2)*S2[1,0], S2[1,1]])
        del S2


        # 3x3 COVARIANCE Matrix elements

        C11 = mlook(np.abs(Kl[0])**2,az,rg)
        C22 = mlook(np.abs(Kl[1])**2,az,rg)
        C33 = mlook(np.abs(Kl[2])**2,az,rg)

        C12 = mlook(Kl[0]*np.conj(Kl[1]),az,rg)
        C13 = mlook(Kl[0]*np.conj(Kl[2]),az,rg)
        C23 = mlook(Kl[1]*np.conj(Kl[2]),az,rg)

        C3Folder = os.path.join(inFolder,'C3')

        if not os.path.isdir(C3Folder):
            os.mkdir(C3Folder)
        
        write_C3(np.dstack([C11,C12,C13,np.conjugate(C12),C22,C23,np.conjugate(C13),np.conjugate(C23),C33]),C3Folder)
        
    else:
        raise ValueError('matrix must be either T3 or C3')
        

    shutil.rmtree(S2Folder)
    
