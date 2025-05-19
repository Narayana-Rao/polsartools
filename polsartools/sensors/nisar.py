
import numpy as np
from osgeo import gdal,osr
import h5py,os,tempfile
from skimage.util.shape import view_as_blocks
from polsartools.utils.utils import time_it
#%%
def read_bin(file):
    ds = gdal.Open(file)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr
def write_rslc_bin(file,wdata):
    
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
def mlook(data,az,rg):
    temp = data[0:data.shape[0]-data.shape[0]%az,0:data.shape[1]-data.shape[1]%rg]
    blocks = view_as_blocks(temp, block_shape=(az, rg))
    flatten = blocks.reshape(blocks.shape[0], blocks.shape[1], -1)
    mean = np.nanmean(flatten, axis=2)
    return mean

def gslc_gtif(data,x_coords,y_coords,projection_epsg,x_res,y_res,output_file,gdal_driver='GTiff'):

    driver = gdal.GetDriverByName(gdal_driver)
    dataset = driver.Create(output_file, len(x_coords), len(y_coords), 1, gdal.GDT_Float32)

    x_min = x_coords[0]
    y_max = y_coords[0]

    geotransform = (x_min, x_res, 0, y_max, 0, y_res)
    dataset.SetGeoTransform(geotransform)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(projection_epsg)
    dataset.SetProjection(srs.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(data)
    dataset.FlushCache()
    dataset = None

def rst_mlook(input_raster, az, rg, output_raster,projection_epsg,
              gdal_driver='GTiff'):
    ds = gdal.Open(input_raster)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    data[data==0]=np.nan

    result = mlook(data, az, rg)
    
    geo_transform = ds.GetGeoTransform()
    projection = ds.GetProjection()

    new_x_size, new_y_size = result.shape[1], result.shape[0]
    new_geo_transform = list(geo_transform)
    new_geo_transform[1] = geo_transform[1] * (data.shape[1] / new_x_size)  
    new_geo_transform[5] = geo_transform[5] * (data.shape[0] / new_y_size)  
    driver = gdal.GetDriverByName(gdal_driver)
    out_ds = driver.Create(output_raster, new_x_size, new_y_size, 1, gdal.GDT_Float32)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(projection_epsg)
    out_ds.SetProjection(srs.ExportToWkt())
    
    out_ds.SetGeoTransform(new_geo_transform)
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(result)
    out_band.FlushCache()
    out_ds = None
    ds = None
    
    return np.shape(result)[0],np.shape(result)[1]

@time_it  
def nisar_gslc(inFile,azlks=20,rglks=10):
    
    inFolder = os.path.dirname(inFile)   
    C2Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2')

    try:
        h5File = h5py.File(inFile,"r")
        xCoordinateSpacing = np.array(h5File['/science/LSAR/GSLC/grids/frequencyA/xCoordinateSpacing'])
        yCoordinateSpacing = np.array(h5File['/science/LSAR/GSLC/grids/frequencyA/yCoordinateSpacing'])
        xCoordinates = np.array(h5File['/science/LSAR/GSLC/grids/frequencyA/xCoordinates'])
        yCoordinates = np.array(h5File['/science/LSAR/GSLC/grids/frequencyA/yCoordinates'])
        projection = np.array(h5File['/science/LSAR/GSLC/metadata/radarGrid/projection'])
    except:
        raise('Invalid GSLC file !!')

    S11 = np.array(h5File['/science/LSAR/GSLC/grids/frequencyA/HH'])
    S12 = np.array(h5File['/science/LSAR/GSLC/grids/frequencyA/HH'])
    
    C11 = np.abs(S11)**2
    C22 = np.abs(S12)**2
    C12 = S11*np.conjugate(S12)
    
    del S11,S12
    
    # tempFile = 'temp.tif'    
    with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tempFile:
        tempFilePath = tempFile.name  # Get the file path

    os.makedirs(C2Folder,exist_ok=True)
    
    gslc_gtif(C11,xCoordinates,yCoordinates,int(projection),
              xCoordinateSpacing,yCoordinateSpacing,tempFilePath,gdal_driver='GTiff')
    del C11
    rst_mlook(tempFilePath, azlks, rglks, os.path.join(C2Folder,'C11.bin'),
              int(projection),gdal_driver='ENVI')
    print(f"Saved file {C2Folder}/C11.bin")
    
    gslc_gtif(C22,xCoordinates,yCoordinates,int(projection),
              xCoordinateSpacing,yCoordinateSpacing,tempFilePath,gdal_driver='GTiff')
    del C22
    rst_mlook(tempFilePath, azlks, rglks, os.path.join(C2Folder,'C22.bin'),
              int(projection),gdal_driver='ENVI')
    print(f"Saved file {C2Folder}/C22.bin")
    
    gslc_gtif(np.real(C12),xCoordinates,yCoordinates,int(projection),
              xCoordinateSpacing,yCoordinateSpacing,tempFilePath,gdal_driver='GTiff')
    
    rst_mlook(tempFilePath, azlks, rglks, os.path.join(C2Folder,'C12_real.bin'),
              int(projection),gdal_driver='ENVI')
    
    print(f"Saved file {C2Folder}/C12_real.bin")
    gslc_gtif(np.imag(C12),xCoordinates,yCoordinates,int(projection),
              xCoordinateSpacing,yCoordinateSpacing,tempFilePath,gdal_driver='GTiff')
    
    rows,cols = rst_mlook(tempFilePath, azlks, rglks, os.path.join(C2Folder,'C12_imag.bin'),
              int(projection),gdal_driver='ENVI')
    
    print(f"Saved file {C2Folder}/C12_imag.bin")
    
    
    file = open(C2Folder +'/config.txt',"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\npp1'%(rows,cols))
    file.close()  
    
    os.remove(tempFilePath)
    h5File.close()
    
@time_it  
def nisar_rslc(inFile,azlks=20,rglks=10):
    
    inFolder = os.path.dirname(inFile)   
    C2Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2')

    try:
        h5File = h5py.File(inFile,"r")
        # sceneCenterAlongTrackSpacing = np.array(h5File['/science/LSAR/RSLC/swaths/frequencyA/sceneCenterAlongTrackSpacing' ])
        # sceneCenterGroundRangeSpacing = np.array(h5File['/science/LSAR/RSLC/swaths/frequencyA/sceneCenterGroundRangeSpacing'])
        # slantRange = np.array(h5File['/science/LSAR/RSLC/swaths/frequencyA/slantRange'])
        # slantRangeSpacing = np.array(h5File['/science/LSAR/RSLC/swaths/frequencyA/slantRangeSpacing'])


        # alongTrackUnitVectorX = np.array(h5File['/science/LSAR/RSLC/metadata/geolocationGrid/alongTrackUnitVectorX' ])
        # alongTrackUnitVectorY = np.array(h5File['/science/LSAR/RSLC/metadata/geolocationGrid/alongTrackUnitVectorY' ])
        # coordinateX = np.array(h5File['/science/LSAR/RSLC/metadata/geolocationGrid/coordinateX' ])
        # coordinateY = np.array(h5File['/science/LSAR/RSLC/metadata/geolocationGrid/coordinateY'])
        # elevationAngle = np.array(h5File['/science/LSAR/RSLC/metadata/geolocationGrid/elevationAngle'])
        # epsg = np.array(h5File['/science/LSAR/RSLC/metadata/geolocationGrid/epsg'])
        
    except:
        raise('Invalid RSLC file !!')

    S11 = np.array(h5File['/science/LSAR/RSLC/swaths/frequencyA/HH'])
    S12 = np.array(h5File['/science/LSAR/RSLC/swaths/frequencyA/HV'])
    
    C11 = np.abs(S11)**2
    C22 = np.abs(S12)**2
    C12 = S11*np.conjugate(S12)
    
    del S11,S12

    os.makedirs(C2Folder,exist_ok=True)
    C11 = mlook(C11,azlks,rglks)
    rows,cols = C11.shape
    write_rslc_bin( os.path.join(C2Folder,'C11.bin'), C11)
    print(f"Saved file {C2Folder}/C11.bin")
    del C11
    
    C22 = mlook(C22,azlks,rglks)   
    write_rslc_bin( os.path.join(C2Folder,'C22.bin'), C22)
    print(f"Saved file {C2Folder}/C22.bin")
    del C22
    
    C12 = mlook(C12,azlks,rglks)   
    write_rslc_bin( os.path.join(C2Folder,'C12_real.bin'), np.real(C12))
    print(f"Saved file {C2Folder}/C12_real.bin")
    write_rslc_bin( os.path.join(C2Folder,'C12_imag.bin'), np.imag(C12))
    print(f"Saved file {C2Folder}/C12_imag.bin")
    del C12
    
    
    file = open(C2Folder +'/config.txt',"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\npp1'%(rows,cols))
    file.close()  
    
    h5File.close()
    
    
#%%



# inFile = "NISAR_L2_PR_GSLC_002_030_A_019_2800_SHNA_A_20081127T061000_20081127T061014_D00404_N_F_J_001.h5"
# # inFile = "NISAR_L1_PR_RSLC_001_001_A_001_2000_SHNA_A_20100612T062402_20100612T062416_T00888_M_F_J_888.h5"

# nisar_gslc(inFile,azlks=20,rglks=10,C2Folder = 'C2')

