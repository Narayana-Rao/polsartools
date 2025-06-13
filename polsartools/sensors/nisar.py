
import numpy as np
from osgeo import gdal,osr
import h5py,os,tempfile
from skimage.util.shape import view_as_blocks
from polsartools.utils.utils import time_it
from polsartools.utils.io_utils import mlook,write_s2_bin_ref, write_s2_ct_ref
#%%
def read_bin(file):
    ds = gdal.Open(file)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr

def nisar_gslc_s2(inFile, S11, S12, S21, S22, xCoordinates, yCoordinates, projection, xCoordinateSpacing, yCoordinateSpacing):
    # Extract folder path
    inFolder = os.path.dirname(inFile) or "."
    S2Folder = os.path.join(inFolder, os.path.basename(inFile).split('.h5')[0], 'S2')
    os.makedirs(S2Folder, exist_ok=True)

    rows, cols = S11.shape  # Determine matrix dimensions

    # Define the matrices to be saved
    matrices = {
        's11.bin': S11,
        's12.bin': (S12 + S21) * 0.5,
        's21.bin': (S12 + S21) * 0.5,
        's22.bin': S22
    }

    for file_name, matrix in matrices.items():
        out_file = os.path.join(S2Folder, file_name)
        gslc_gtif(matrix, xCoordinates, yCoordinates, int(projection),
                  xCoordinateSpacing, yCoordinateSpacing, out_file,
                  gdal_driver='ENVI', dtype=gdal.GDT_CFloat64)
        print(f"Saved file {out_file}")

    # Write configuration file
    config_file_path = os.path.join(S2Folder, 'config.txt')
    with open(config_file_path, "w") as file:
        file.write(f'Nrow\n{rows}\n---------\nNcol\n{cols}\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull')

def nisar_gslc_ct(matrix_type, matrixFolder, K, xCoordinates, yCoordinates, projection, xCoordinateSpacing, yCoordinateSpacing, azlks, rglks):
    if matrix_type not in ["C3", "T3", "C4", "T4"]:
        raise ValueError("Invalid matrix type. Choose 'C3', 'T3', 'C4', or 'T4'.")
    
    prefix = matrix_type[0]  # Extract first character (C or T)
    
    # # Extract folder path
    # inFolder = os.path.dirname(inFile) or "."
    # matrixFolder = os.path.join(inFolder, os.path.basename(inFile).split('.h5')[0], matrix_type)
    # os.makedirs(matrixFolder, exist_ok=True)

    # Compute matrix elements
    matrices = {
        f'{prefix}11': np.abs(K[0])**2,
        f'{prefix}22': np.abs(K[1])**2,
        f'{prefix}33': np.abs(K[2])**2,
        f'{prefix}12_real': np.real(K[0] * np.conj(K[1])),
        f'{prefix}12_imag': np.imag(K[0] * np.conj(K[1])),
        f'{prefix}13_real': np.real(K[0] * np.conj(K[2])),
        f'{prefix}13_imag': np.imag(K[0] * np.conj(K[2])),
        f'{prefix}23_real': np.real(K[1] * np.conj(K[2])),
        f'{prefix}23_imag': np.imag(K[1] * np.conj(K[2])),
    }

    # Extend for 4-element case (C4/T4)
    if len(K) == 4:
        matrices.update({
            f'{prefix}44': np.abs(K[3])**2,
            f'{prefix}14_real': np.real(K[0] * np.conj(K[3])),
            f'{prefix}14_imag': np.imag(K[0] * np.conj(K[3])),
            f'{prefix}24_real': np.real(K[1] * np.conj(K[3])),
            f'{prefix}24_imag': np.imag(K[1] * np.conj(K[3])),
            f'{prefix}34_real': np.real(K[2] * np.conj(K[3])),
            f'{prefix}34_imag': np.imag(K[2] * np.conj(K[3])),
        })

    with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tempFile:
        tempFilePath = tempFile.name

    rows, cols = None, None
    for key, matrix in matrices.items():
        gslc_gtif(matrix.astype(np.float32 if 'real' not in key and 'imag' not in key else np.complex64),
                  xCoordinates, yCoordinates, int(projection),
                  xCoordinateSpacing, yCoordinateSpacing, tempFilePath, gdal_driver='GTiff')

        rows, cols = mlook_geo(tempFilePath, azlks, rglks, os.path.join(matrixFolder, f"{key}.bin"),
                               int(projection), gdal_driver='ENVI')

        print(f"Saved file {matrixFolder}/{key}.bin")

    # Write configuration file
    config_file_path = os.path.join(matrixFolder, 'config.txt')
    with open(config_file_path, "w") as file:
        file.write(f'Nrow\n{rows}\n---------\nNcol\n{cols}\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull')

    os.remove(tempFilePath)


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

def gslc_gtif(data,x_coords,y_coords,projection_epsg,x_res,y_res,output_file,gdal_driver='GTiff',dtype=gdal.GDT_Float32):

    driver = gdal.GetDriverByName(gdal_driver)
    dataset = driver.Create(output_file, len(x_coords), len(y_coords), 1, dtype)

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

def mlook_geo(input_raster, az, rg, output_raster,projection_epsg,
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
def nisar_gslc(inFile,azlks=22,rglks=10,matrixType='C3'):
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
        # print("Detected L-band data ")

    elif '/science/SSAR' in h5File:
        freq_band = 'S'
        # print(" Detected S-band data")

    else:
        print("Neither LSAR nor SSAR data found in the file.")
        h5File.close()
        return
    
    
    xCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/xCoordinateSpacing'])
    yCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/yCoordinateSpacing'])
    xCoordinates = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/xCoordinates'])
    yCoordinates = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/yCoordinates'])
    projection = np.array(h5File[f'/science/{freq_band}SAR/GSLC/metadata/radarGrid/projection'])

    
    listOfPolarizations = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/listOfPolarizations']).astype(str)
    nchannels = len(listOfPolarizations)
    print(f"Detected {freq_band}-band {listOfPolarizations} ")
    
    
    # S11 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/HH'])
    # S12 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/HV'])
    
    
    if nchannels==2:
        print("Extracting C2 matrix elements...")
        if 'HH' in listOfPolarizations and 'HV' in listOfPolarizations:
            S11 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/HH'])
            S12 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/HV'])
        elif 'VV' in listOfPolarizations and 'VH' in listOfPolarizations:
            S11 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/VV'])
            S12 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/VH'])
        elif 'HH' in listOfPolarizations and 'VV' in listOfPolarizations:
            S11 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/HH'])
            S12 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/VV'])
        else:
            print("No HH, HV, VV, or VH polarizations found in the file.")
            h5File.close()
            return

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
        mlook_geo(tempFilePath, azlks, rglks, os.path.join(C2Folder,'C11.bin'),
                int(projection),gdal_driver='ENVI')
        print(f"Saved file {C2Folder}/C11.bin")
        
        gslc_gtif(C22,xCoordinates,yCoordinates,int(projection),
                xCoordinateSpacing,yCoordinateSpacing,tempFilePath,gdal_driver='GTiff')
        del C22
        mlook_geo(tempFilePath, azlks, rglks, os.path.join(C2Folder,'C22.bin'),
                int(projection),gdal_driver='ENVI')
        print(f"Saved file {C2Folder}/C22.bin")
        
        gslc_gtif(np.real(C12),xCoordinates,yCoordinates,int(projection),
                xCoordinateSpacing,yCoordinateSpacing,tempFilePath,gdal_driver='GTiff')
        
        mlook_geo(tempFilePath, azlks, rglks, os.path.join(C2Folder,'C12_real.bin'),
                int(projection),gdal_driver='ENVI')
        
        print(f"Saved file {C2Folder}/C12_real.bin")
        gslc_gtif(np.imag(C12),xCoordinates,yCoordinates,int(projection),
                xCoordinateSpacing,yCoordinateSpacing,tempFilePath,gdal_driver='GTiff')
        
        rows,cols = mlook_geo(tempFilePath, azlks, rglks, os.path.join(C2Folder,'C12_imag.bin'),
                int(projection),gdal_driver='ENVI')
        
        print(f"Saved file {C2Folder}/C12_imag.bin")
        
        
        file = open(C2Folder +'/config.txt',"w+")
        file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\npp1'%(rows,cols))
        file.close()  
        
        os.remove(tempFilePath)
        h5File.close()
    elif nchannels==4:
        # print("Detected full polarimetric data")
        S11 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/HH'])
        S12 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/HV'])        
        S21 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/VH'])
        S22 = np.array(h5File[f'/science/{freq_band}SAR/GSLC/grids/frequencyA/VV'])
        
        h5File.close()
        if matrixType=='S2':
            print("Extracting S2 matrix elements...")
            nisar_gslc_s2(inFile, S11, S12, S21, S22,
                          xCoordinates, yCoordinates, projection, 
                          xCoordinateSpacing, yCoordinateSpacing)
            
        elif matrixType=='C3':
            
            print("Extracting C3 matrix elements...")
            
            inFolder = os.path.dirname(inFile)   
            if not inFolder:
                inFolder = "."
            C3Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C3')
            os.makedirs(C3Folder,exist_ok=True)
            
            # 3x3 covariance Matrix
            nisar_gslc_ct("C3",C3Folder, np.array([S11, np.sqrt(2)*0.5*(S12+S21), S22]), 
                        xCoordinates, yCoordinates, 
                            projection, xCoordinateSpacing, yCoordinateSpacing, 
                            azlks, rglks)
            
            # os.remove(tempFilePath)
        elif matrixType=='T3':
            print("Extracting T3 matrix elements...")
            # 3x1 Pauli Coherency Matrix
            # Kp = (1/np.sqrt(2))*np.array([S11+S22, S11-S22, S12+S21])
            inFolder = os.path.dirname(inFile)   
            if not inFolder:
                inFolder = "."
            T3Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T3')
            os.makedirs(T3Folder,exist_ok=True)
            
            nisar_gslc_ct("T3",T3Folder, (1/np.sqrt(2))*np.array([S11+S22, S11-S22, S12+S21]), 
            xCoordinates, yCoordinates, 
                projection, xCoordinateSpacing, yCoordinateSpacing, 
                azlks, rglks)
            
        elif matrixType=='C4':
            print("Extracting C4 matrix elements...")
            
            inFolder = os.path.dirname(inFile)   
            if not inFolder:
                inFolder = "."
            C4Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C4')
            os.makedirs(C4Folder,exist_ok=True)
            
            # 4x4 covariance Matrix
            nisar_gslc_ct("C4",C4Folder, np.array([S11, S12, S21, S22]), 
            xCoordinates, yCoordinates, 
                projection, xCoordinateSpacing, yCoordinateSpacing, 
                azlks, rglks)
            
        elif matrixType=='T4':
            print("Extracting T4 matrix elements...")
            
            inFolder = os.path.dirname(inFile)   
            if not inFolder:
                inFolder = "."
            T4Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T4')
            os.makedirs(T4Folder,exist_ok=True)
            
            # 4x4 covariance Matrix
            nisar_gslc_ct("T4",T4Folder, (1/np.sqrt(2))*np.array([S11+S22, S11-S22, S12+S21, 1j*(S12-S21)]), 
            xCoordinates, yCoordinates, 
                projection, xCoordinateSpacing, yCoordinateSpacing, 
                azlks, rglks)
           
            
        else:
            raise ValueError(f"Unknown matrix type {matrixType} \n Valid options are 'S2', 'T3', 'C3'")
    else:
        print("No polarimetric data found in the file.")
        h5File.close()
        return
        
@time_it  
def nisar_rslc(inFile,azlks=22,rglks=10, matrixType='C3'):
    """
    Extracts the C2 (for dual-pol), S2/C3/T3 (for full-pol) matrix elements from a NISAR RSLC HDF5 file 
    and saves them into respective binary files.
    
    Example:
    --------
    >>> nisar_rslc("path_to_file.h5", azlks=30, rglks=15)
    This will extract the C2 (for dual-pol), S2/C3/T3 (for full-pol) matrix elements from the specified NISAR RSLC file 
    and save them in the respective folders.
    
    Parameters:
    -----------
    inFile : str
        The path to the NISAR RSLC HDF5 file containing the radar data.

    azlks : int, optional (default=22)
        The number of azimuth looks for multi-looking. 

    rglks : int, optional (default=10)
        The number of range looks for multi-looking. 
        
    matrixType : str, optional (default = C2 for dual-pol, C3 for full-pol)
        Type of matrix to extract. Valid options are 'C2', 'S2', 'C3', and 'T3'.
         

    Returns:
    --------
    None
        The function does not return any value. Instead, it creates a folder 
        named `C2/S2/C3/T3` (if not already present) and saves the following binary files:


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
        # print("Detected L-band data ")
    elif '/science/SSAR' in h5File:
        freq_band = 'S' 
        # print(" Detected S-band data")
    else:
        print("Neither LSAR nor SSAR data found in the file.")
        h5File.close()
        return
    
    # xCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinateSpacing'])
    # yCoordinateSpacing = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinateSpacing'])
    # xCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/xCoordinates'])
    # yCoordinates = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/yCoordinates'])
    # projection = np.array(h5File[f'/science/{freq_band}SAR/RSLC/metadata/radarGrid/projection'])
    
    if f'/science/{freq_band}SAR/RSLC' in h5File:
        base_path = f'/science/{freq_band}SAR/RSLC/swaths/frequencyA'
    elif f'/science/{freq_band}SAR/SLC' in h5File:
        base_path = f'/science/{freq_band}SAR/SLC/swaths/frequencyA'
    else:
        print("Neither RSLC nor SLC found in the HDF5 structure.")
        h5File.close()
        return
    
    
    listOfPolarizations = np.array(h5File[f'{base_path}/listOfPolarizations']).astype(str)
    nchannels = len(listOfPolarizations)

    print(f"Detected {freq_band}-band {listOfPolarizations} ")
    
    if nchannels==2:
        print("Extracting C2 matrix elements...")
        if 'HH' in listOfPolarizations and 'HV' in listOfPolarizations:
            S11 = np.array(h5File[f'{base_path}/HH'])
            S12 = np.array(h5File[f'{base_path}/HV'])
        elif 'VV' in listOfPolarizations and 'VH' in listOfPolarizations:
            S11 = np.array(h5File[f'{base_path}/VV'])
            S12 = np.array(h5File[f'{base_path}/VH'])
        elif 'HH' in listOfPolarizations and 'VV' in listOfPolarizations:
            S11 = np.array(h5File[f'{base_path}/HH'])
            S12 = np.array(h5File[f'{base_path}/VV'])
        else:
            print("No HH, HV, VV, or VH polarizations found in the file.")
            h5File.close()
            return
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
    elif nchannels==4:
        
        S11 = h5File[f'{base_path}/HH']
        S12 = h5File[f'{base_path}/HV']
        S21 = h5File[f'{base_path}/VH']
        S22 = h5File[f'{base_path}/VV']

        if matrixType == 'S2':
            print("Extracting S2 matrix elements...")
            rows, cols = S11.shape
            inFolder = os.path.dirname(inFile) or "."
            out_dir = os.path.join(inFolder, os.path.basename(inFile).split('.h5')[0], 'S2')
            os.makedirs(out_dir, exist_ok=True)

            # Lazy loading: Extract data only at time of writing
            write_s2_bin_ref(os.path.join(out_dir, 's11.bin'), S11[()])
            write_s2_bin_ref(os.path.join(out_dir, 's12.bin'), (S12[()] + S21[()]) / 2)
            write_s2_bin_ref(os.path.join(out_dir, 's21.bin'), (S12[()] + S21[()]) / 2)
            write_s2_bin_ref(os.path.join(out_dir, 's22.bin'), S22[()])

            # Save config file
            with open(os.path.join(out_dir, 'config.txt'), "w") as file:
                file.write(f'Nrow\n{rows}\n---------\nNcol\n{cols}\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull')
            print(f"Saved config file {out_dir}/config.txt")

        elif matrixType in ["C3", "T3", "C4", "T4"]:
            folder_name = matrixType
            print(f"Extracting {matrixType} matrix elements...")

            inFolder = os.path.dirname(inFile) or "."
            matrixFolder = os.path.join(inFolder, os.path.basename(inFile).split('.h5')[0], folder_name)
            os.makedirs(matrixFolder, exist_ok=True)

            # Pass lazy-loaded datasets (do NOT convert them to NumPy arrays)
            if matrixType == "C3":
                write_s2_ct_ref("C3", matrixFolder, [S11, np.sqrt(2) * 0.5 * (S12[()] + S21[()]), S22], azlks, rglks)
            elif matrixType == "T3":
                write_s2_ct_ref("T3", matrixFolder, np.array([
                                                                S11[()] + S22[()], 
                                                                S11[()] - S22[()], 
                                                                S12[()] + S21[()], 
                                                            ]) * (1 / np.sqrt(2)), azlks, rglks)
            elif matrixType == "C4":
                write_s2_ct_ref("C4", matrixFolder, [S11, S12, S21, S22], azlks, rglks)
            elif matrixType == "T4":
                write_s2_ct_ref("T4", matrixFolder, np.array([
                                                                S11[()] + S22[()], 
                                                                S11[()] - S22[()], 
                                                                S12[()] + S21[()], 
                                                                1j * (np.array(S12[()]).astype(np.complex64) - np.array(S21[()]).astype(np.complex64))
                                                            ]) * (1 / np.sqrt(2)), azlks, rglks)
        
        
        # print("Detected full polarimetric data")
        # S11 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HH'])
        # S12 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/HV'])        
        # S21 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VH'])
        # S22 = np.array(h5File[f'/science/{freq_band}SAR/RSLC/swaths/frequencyA/VV'])
        
        # if matrixType=='S2':
        #     print("Extracting S2 matrix elements...")
        #     rows,cols = S11.shape
        #     inFolder = os.path.dirname(inFile)   
        #     if not inFolder:
        #         inFolder = "."
        #     out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'S2')
        #     os.makedirs(out_dir,exist_ok=True)
        
            
        #     out_file = os.path.join(out_dir,'s11.bin')
        #     write_s2_bin(out_file,S11)
        #     print("Saved file "+out_file)
            
        #     out_file = os.path.join(out_dir,'s12.bin')
        #     write_s2_bin(out_file,(S12+S21)/2)
        #     print("Saved file "+out_file)
            
        #     out_file = os.path.join(out_dir,'s21.bin')
        #     write_s2_bin(out_file,(S12+S21)/2)
        #     print("Saved file "+out_file)
            
        #     out_file = os.path.join(out_dir,'s22.bin')
        #     write_s2_bin(out_file,S22)
        #     print("Saved file "+out_file)
            
        #     file = open(out_dir +'/config.txt',"w+")
        #     file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
        #     file.close() 
        
        # elif matrixType=='C3':
        #     print("Extracting C3 matrix elements...")
        #     inFolder = os.path.dirname(inFile)   
        #     if not inFolder:
        #         inFolder = "."
        #     C3Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C3')
        #     os.makedirs(C3Folder,exist_ok=True)
            
        #     write_s2_ct("C3", C3Folder, np.array([S11, np.sqrt(2)*0.5*(S12+S21), S22]), azlks, rglks)

        # elif matrixType=='T3':
            
        #     print("Extracting T3 matrix elements...")
            
        #     inFolder = os.path.dirname(inFile)   
        #     if not inFolder:
        #         inFolder = "."
        #     T3Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T3')
        #     os.makedirs(T3Folder,exist_ok=True)
            
        #     write_s2_ct("T3", T3Folder, (1/np.sqrt(2))*np.array([S11+S22, S11-S22, S12+S21]), azlks, rglks)

        # elif matrixType=='C4':
        #     print("Extracting C4 matrix elements...")
        #     inFolder = os.path.dirname(inFile)   
        #     if not inFolder:
        #         inFolder = "."
        #     C4Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C4')
        #     os.makedirs(C4Folder,exist_ok=True)
            
        #     write_s2_ct("C4", C4Folder, np.array([S11, S12, S21, S22]), azlks, rglks)

        # elif matrixType=='T4':
            
        #     print("Extracting T4 matrix elements...")
            
        #     inFolder = os.path.dirname(inFile)   
        #     if not inFolder:
        #         inFolder = "."
        #     T4Folder = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T4')
        #     os.makedirs(T4Folder,exist_ok=True)
            
        #     write_s2_ct("T4", T4Folder, (1/np.sqrt(2))*np.array([S11+S22, S11-S22, S12+S21, 1j*(S12-S21)]), azlks, rglks)
            
    else:
        print("No polarimetric data found in the file.")
        h5File.close()
        return
        
    h5File.close()
    

    

    
