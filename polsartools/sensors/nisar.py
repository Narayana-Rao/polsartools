
import numpy as np
from osgeo import gdal,osr
import h5py,os,tempfile
from skimage.util.shape import view_as_blocks
from polsartools.utils.utils import time_it
from polsartools.utils.io_utils import mlook

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
def nisar_gslc(inFile,azlks=22,rglks=10):
    """
    Extracts the C2 matrix elements (C11, C22, and C12) from a NISAR GSLC HDF5 file 
    and saves them into respective binary files.

    Example:
    --------
    >>> nisar_gslc("path_to_file.h5", azlks=30, rglks=15)
    This will extract the C2 matrix elements from the specified NISAR GSLC file 
    and save them in the 'C2' folder.
    
    Parameters:
    -----------
    inFile : str
        The path to the NISAR GSLC HDF5 file containing the radar data.

    azlks : int, optional (default=20)
        The number of azimuth looks for multi-looking. 

    rglks : int, optional (default=10)
        The number of range looks for multi-looking. 

    Returns:
    --------
    None
        The function does not return any value. Instead, it creates a folder 
        named `C2` (if not already present) and saves the following binary files:
        
        - `C11.bin`: Contains the C11 matrix elements.
        - `C22.bin`: Contains the C22 matrix elements.
        - `C12_real.bin`: Contains the real part of the C12 matrix.
        - `C12_imag.bin`: Contains the imaginary part of the C12 matrix.
        - `config.txt`: A text file containing grid dimensions and polarimetric configuration.

    Raises:
    -------
    Exception
        If the GSLC HDF5 file is invalid or cannot be read.


    """
    inFolder = os.path.dirname(inFile)   
    if not inFolder:
        inFolder = "."
    C2Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2')

    try:
        h5File = h5py.File(inFile,"r")
    except:
        raise('Invalid .h5 file !!')
    
    
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
    
    
    xCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/xCoordinateSpacing'])
    yCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/yCoordinateSpacing'])
    xCoordinates = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/xCoordinates'])
    yCoordinates = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/yCoordinates'])
    projection = np.array(h5File[f'/science/{freq_band}SAR/GSLC/metadata/radarGrid/projection'])
    S11 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/HH'])
    S12 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/HV'])
    listOfPolarizations = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/listOfPolarizations']).astype(str)
    len(listOfPolarizations)
    h5File.close()

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
def nisar_rslc(inFile,azlks=22,rglks=10):
    """
    Extracts the C2 matrix elements (C11, C22, and C12) from a NISAR RSLC HDF5 file 
    and saves them into respective binary files.
    
    Example:
    --------
    >>> nisar_rslc("path_to_file.h5", azlks=30, rglks=15)
    This will extract the C2 matrix elements from the specified NISAR RSLC file 
    and save them in the 'C2' folder.
    
    Parameters:
    -----------
    inFile : str
        The path to the NISAR RSLC HDF5 file containing the radar data.

    azlks : int, optional (default=20)
        The number of azimuth looks for multi-looking. 

    rglks : int, optional (default=10)
        The number of range looks for multi-looking. 

    Returns:
    --------
    None
        The function does not return any value. Instead, it creates a folder 
        named `C2` (if not already present) and saves the following binary files:
        
        - `C11.bin`: Contains the C11 matrix elements.
        - `C22.bin`: Contains the C22 matrix elements.
        - `C12_real.bin`: Contains the real part of the C12 matrix.
        - `C12_imag.bin`: Contains the imaginary part of the C12 matrix.
        - `config.txt`: A text file containing grid dimensions and polarimetric configuration.

    Raises:
    -------
    Exception
        If the RSLC HDF5 file is invalid or cannot be read.


    """
    
    inFolder = os.path.dirname(inFile)
    if not inFolder:
        inFolder = "."   
    C2Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2')
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
    
    
    
    S11 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HH'])
    S12 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HV'])
    
    # xCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinateSpacing'])
    # yCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinateSpacing'])
    # xCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinates'])
    # yCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinates'])
    # projection = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/radarGrid/projection'])
    listOfPolarizations = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/listOfPolarizations']).astype(str)
    len(listOfPolarizations)
    
    
    h5File.close()
    
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
    
