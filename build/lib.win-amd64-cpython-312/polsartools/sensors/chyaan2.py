import glob,os
import numpy as np
from osgeo import gdal 

from polsartools.utils.utils import time_it
from polsartools.utils.io_utils import mlook, write_T3, write_C3


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


def read_bin(file):
    ds = gdal.Open(file,gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr


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
def chyaan2_fp(inFolder,matrix='T3',azlks=None,rglks=None):
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
    # multi-llok factor 
    mlf = int(np.round(gRange/ols,0))

    ds = gdal.Open(glob.glob(inFolder+'/data/calibrated/*/*sli*_hh_*.tif')[0])
    cols = ds.RasterXSize  
    rows = ds.RasterYSize 


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

    calFactor = 1/np.sqrt(10**(cc/10))

    if matrix == 'S2':

        out_dir = os.path.join(inFolder,"S2")
        os.makedirs(out_dir,exist_ok=True)

        print("Considering S12 = S21")

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_hh_*.tif'))[0]
        data = read_rs2_tif(inFile)
        out_file = os.path.join(out_dir,'s11.bin')
        write_s2_bin(out_file,data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor))
        print("Saved file "+out_file)
        
        rows, cols, _ = data.shape

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_hv_*.tif'))[0]
        data_xy = read_rs2_tif(inFile)


        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_vh_*.tif'))[0]
        data_yx = read_rs2_tif(inFile)

        data = (data_xy+data_yx)*0.5
        del data_xy,data_yx

        out_file = os.path.join(out_dir,'s12.bin')
        
        write_s2_bin(out_file,data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor))
        print("Saved file "+out_file)
        out_file = os.path.join(out_dir,'s21.bin')
        write_s2_bin(out_file,data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor))
        print("Saved file "+out_file)

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_vv_*.tif'))[0]
        data = read_rs2_tif(inFile)
        out_file = os.path.join(out_dir,'s22.bin')
        write_s2_bin(out_file,data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor))
        print("Saved file "+out_file)
        
        file = open(out_dir +'/config.txt',"w+")
        file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
        file.close() 
        
    elif matrix == 'T3':
        print("Considering S12 = S21")
        
        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_hh_*.tif'))[0]
        data = read_rs2_tif(inFile)
        s11 = data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor)

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_hv_*.tif'))[0]
        data_xy = read_rs2_tif(inFile)

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_vh_*.tif'))[0]
        data_yx = read_rs2_tif(inFile)
        
        # Symmetry assumption
        data = (data_xy+data_yx)*0.5
        del data_xy,data_yx
        s12 = data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor)

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_vv_*.tif'))[0]
        data = read_rs2_tif(inFile)
        s22 = data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor)
        
        # Kp- 3-D Pauli feature vector
        Kp = (1/np.sqrt(2))*np.array([s11+s22, s11-s22, 2*s12])

        del s11,s12,s22


        if azlks == None and rglks == None:
            azlks = mlf
            rglks = 1
                
        print(f'Using multi-look factor: azlks = {azlks}, rglks = {rglks}')

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

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_hh_*.tif'))[0]
        data = read_rs2_tif(inFile)
        s11 = data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor)

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_hv_*.tif'))[0]
        data_xy = read_rs2_tif(inFile)

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_vh_*.tif'))[0]
        data_yx = read_rs2_tif(inFile)
        
        # Symmetry assumption
        data = (data_xy+data_yx)*0.5
        del data_xy,data_yx
        s12 = data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor)

        inFile = glob.glob(os.path.join(inFolder, 'data/calibrated/*/*sli*_vv_*.tif'))[0]
        data = read_rs2_tif(inFile)
        s22 = data[:,:,0]*calFactor+1j*(data[:,:,1]*calFactor)

        # Kl- 3-D Lexicographic feature vector
        Kl = np.array([s11, np.sqrt(2)*s12, s22])
        del s11,s12,s22


        if azlks == None and rglks == None:
            azlks = mlf
            rglks = 1
                
        print(f'Using multi-look factor: azlks = {azlks}, rglks = {rglks}')

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

