
import numpy as np
from osgeo import gdal,osr
import os,tempfile
import tables
from skimage.util.shape import view_as_blocks
from polsartools.utils.utils import time_it, mlook_arr
from polsartools.utils.io_utils import write_s2_bin_ref, write_s2_ct_ref
from polsartools.utils.h5_utils import h5_polsar
from netCDF4 import Dataset
#%%

def get_ml_chunk(multiplier, default_size):
    # Rounds up to the next multiple of `multiplier`
    return ((default_size + multiplier - 1) // multiplier) * multiplier
def gslc_meta(inFile):
    band_table = [
        ('/science/LSAR', 'L'),
        ('/science/SSAR', 'S')
    ]

    # Step 1: Identify frequency band and root path
    try:
        with tables.open_file(inFile, mode="r") as h5:
            for path, band in band_table:
                if path in h5:
                    freq_band = band
                    freq_path = path
                    break
            else:
                print("Neither LSAR nor SSAR data found in the file.")
                return None

            # Step 2: Read polarization list
            pol_path = f'{freq_path}/GSLC/grids/frequencyA/listOfPolarizations'
            try:
                listOfPolarizations = np.array(h5.get_node(pol_path).read()).astype(str)
            except tables.NoSuchNodeError:
                print(f"Polarization node missing: {pol_path}")
                listOfPolarizations = None

    except Exception:
        raise RuntimeError("Invalid .h5 file !!")

    # Step 3: Reopen to read raster metadata
    with tables.open_file(inFile, "r") as h5:
        base_grid = f'{freq_path}/GSLC/grids/frequencyA'
        projection_path = f'{freq_path}/GSLC/metadata/radarGrid/projection'

        try:
            projection = np.array(h5.get_node(projection_path).read())
            xSpacing = np.array(h5.get_node(base_grid + '/xCoordinateSpacing').read())
            ySpacing = np.array(h5.get_node(base_grid + '/yCoordinateSpacing').read())
        except tables.NoSuchNodeError as e:
            print(f"⚠️ Missing expected metadata node: {e}")
            return None

    return freq_band,listOfPolarizations, xSpacing, ySpacing, int(projection)

@time_it 
def nisar_gslc(inFile, matrixType='C3', azlks=2, rglks=2, outType='tif',max_workers=None):
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
        
    
    freq_band,listOfPolarizations, xres, yres, projection = gslc_meta(inFile)
    nchannels = len(listOfPolarizations)
    print(f"Detected {freq_band}-band polarization channels: {listOfPolarizations}")
    
    ds = Dataset(inFile, "r")
    xcoords = ds.groups["science"].groups[f"{freq_band}SAR"].groups["GSLC"] \
                   .groups["grids"].groups["frequencyA"] \
                   .variables["xCoordinates"][:]
    ycoords = ds.groups["science"].groups[f"{freq_band}SAR"].groups["GSLC"] \
                   .groups["grids"].groups["frequencyA"] \
                   .variables["yCoordinates"][:]                
    ds=None

    inshape = [len(xcoords),len(ycoords)]
    outshape = [len(xcoords)//rglks,len(ycoords)//azlks]
    
    start_x = min(xcoords)
    start_y = max(ycoords)
    
    # print(projection,start_x,start_y,xres,yres,inshape,outshape)
    # print(min(ycoords),max(ycoords), min(ycoords)+np.abs(yres)*(inshape[1]-1))
    
    inFolder = os.path.dirname(inFile)   
    if not inFolder:
        inFolder = "."
    
    base_path = f'/science/{freq_band}SAR/GSLC/grids/frequencyA'
    
    if nchannels==2:
        out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2')
        temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2','temp')
        print("Extracting C2 matrix elements...")
        if 'HH' in listOfPolarizations and 'HV' in listOfPolarizations:
            h5_polsar(
                        h5_file=inFile,
                        dataset_paths={
                            "HH": f"{base_path}/HH",
                            "HV": f"{base_path}/HV",
                        },
                        output_dir=out_dir,
                        temp_dir=temp_dir,
                        azlks=azlks,
                        rglks=rglks,
                        matrix_type = 'C2HX',
                        apply_multilook=True,
                        chunk_size_x=get_ml_chunk(rglks, 512),
                        chunk_size_y=get_ml_chunk(azlks, 512),
                        max_workers=max_workers,
                        start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
                        outType=outType,
                        dtype = np.float32,
                        inshape=inshape,
                        outshape=outshape
                    )
    

        elif 'VV' in listOfPolarizations and 'VH' in listOfPolarizations:
            h5_polsar(
                        h5_file=inFile,
                        dataset_paths={
                            "VV": f"{base_path}/VV",
                            "VH": f"{base_path}/VH",
                        },
                        output_dir=out_dir,
                        temp_dir=temp_dir,
                        azlks=azlks,
                        rglks=rglks,
                        matrix_type = 'C2VX',
                        apply_multilook=True,
                        chunk_size_x=get_ml_chunk(rglks, 512),
                        chunk_size_y=get_ml_chunk(azlks, 512),
                        max_workers=max_workers,
                        start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
                        outType=outType,
                        dtype = np.float32,
                        inshape=inshape,
                        outshape=outshape
                    )

        elif 'HH' in listOfPolarizations and 'VV' in listOfPolarizations:
            h5_polsar(
                        h5_file=inFile,
                        dataset_paths={
                            "HH": f"{base_path}/HH",
                            "VV": f"{base_path}/VV",
                        },
                        output_dir=out_dir,
                        temp_dir=temp_dir,
                        azlks=azlks,
                        rglks=rglks,
                        matrix_type = 'C2HV',
                        apply_multilook=True,
                        chunk_size_x=get_ml_chunk(rglks, 512),
                        chunk_size_y=get_ml_chunk(azlks, 512),
                        max_workers=max_workers,
                        start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
                        outType=outType,
                        dtype = np.float32,
                        inshape=inshape,
                        outshape=outshape
                    )

        else:
            print("No HH, HV, VV, or VH polarizations found in the file.")

            return
                
    elif nchannels==4:
        if matrixType=='S2':
            print("Extracting S2 matrix elements...")
            out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'S2')
            temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'S2','temp')
            h5_polsar(
                h5_file=inFile,
                dataset_paths={
                    "HH": f"{base_path}/HH",
                    "HV": f"{base_path}/HV",
                    "VH": f"{base_path}/VH",
                    "VV": f"{base_path}/VV",
                },
                output_dir=out_dir,
                temp_dir=temp_dir,
                azlks=azlks,
                rglks=rglks,
                matrix_type = 'S2',
                apply_multilook=False,
                chunk_size_x=get_ml_chunk(rglks, 512),
                chunk_size_y=get_ml_chunk(azlks, 512),
                max_workers=max_workers,
                start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
                outType=outType,
                dtype = np.complex64,
                inshape=inshape,
                outshape=outshape
            )
            
        elif matrixType=='T4':
            print("Extracting T4 matrix elements...")
            out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T4')
            temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T4','temp')
            h5_polsar(
                h5_file=inFile,
                dataset_paths={
                    "HH": f"{base_path}/HH",
                    "HV": f"{base_path}/HV",
                    "VH": f"{base_path}/VH",
                    "VV": f"{base_path}/VV",
                },
                output_dir=out_dir,
                temp_dir=temp_dir,
                azlks=azlks,
                rglks=rglks,
                matrix_type = 'T4',
                apply_multilook=True,
                chunk_size_x=get_ml_chunk(rglks, 512),
                chunk_size_y=get_ml_chunk(azlks, 512),
                max_workers=max_workers,
                start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
                outType=outType,
                dtype = np.float32,
                inshape=inshape,
                outshape=outshape
            )

        elif matrixType=='T3':
            print("Extracting T3 matrix elements...")
            out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T3')
            temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T3','temp')
            h5_polsar(
                h5_file=inFile,
                dataset_paths={
                    "HH": f"{base_path}/HH",
                    "HV": f"{base_path}/HV",
                    "VH": f"{base_path}/VH",
                    "VV": f"{base_path}/VV",
                },
                output_dir=out_dir,
                temp_dir=temp_dir,
                azlks=azlks,
                rglks=rglks,
                matrix_type = 'T3',
                apply_multilook=True,
                chunk_size_x=get_ml_chunk(rglks, 512),
                chunk_size_y=get_ml_chunk(azlks, 512),
                max_workers=max_workers,
                start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=int(projection),
                outType=outType,
                dtype = np.float32,
                inshape=inshape,
                outshape=outshape
            )
        elif matrixType=='C4':
            print("Extracting C4 matrix elements...")
            out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C4')
            temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C4','temp')
            h5_polsar(
                h5_file=inFile,
                dataset_paths={
                    "HH": f"{base_path}/HH",
                    "HV": f"{base_path}/HV",
                    "VH": f"{base_path}/VH",
                    "VV": f"{base_path}/VV",
                },
                output_dir=out_dir,
                temp_dir=temp_dir,
                azlks=azlks,
                rglks=rglks,
                matrix_type = 'C4',
                apply_multilook=True,
                chunk_size_x=get_ml_chunk(rglks, 512),
                chunk_size_y=get_ml_chunk(azlks, 512),
                max_workers=max_workers,
                start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=int(projection),
                outType=outType,
                dtype = np.float32,
                inshape=inshape,
                outshape=outshape            
            )   
        elif matrixType=='C3':
            print("Extracting C3 matrix elements...")
            out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C3')
            temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C3','temp')
            h5_polsar(
                h5_file=inFile,
                dataset_paths={
                    "HH": f"{base_path}/HH",
                    "HV": f"{base_path}/HV",
                    "VH": f"{base_path}/VH",
                    "VV": f"{base_path}/VV",
                },
                output_dir=out_dir,
                temp_dir=temp_dir,
                azlks=azlks,
                rglks=rglks,
                matrix_type = 'C3',
                apply_multilook=True,
                chunk_size_x=get_ml_chunk(rglks, 512),
                chunk_size_y=get_ml_chunk(azlks, 512),
                max_workers=max_workers,
                start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=int(projection),
                outType=outType,            
                dtype = np.float32,
                inshape=inshape,
                outshape=outshape            
            )
        
        elif matrixType=='C2HV':
            print("Extracting C2HV matrix elements...")
            out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2HV')
            temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2HV','temp')
            h5_polsar(
                h5_file=inFile,
                dataset_paths={
                    "HH": f"{base_path}/HH",
                    "VV": f"{base_path}/VV",
                },
                output_dir=out_dir, temp_dir=temp_dir,
                azlks=azlks, rglks=rglks, matrix_type = 'C2HV', apply_multilook=True,
                chunk_size_x=get_ml_chunk(rglks, 512), chunk_size_y=get_ml_chunk(azlks, 512), max_workers=max_workers,
                start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
                outType=outType, dtype = np.float32, inshape=inshape,outshape=outshape            
            )

        elif matrixType=='C2HX':
            print("Extracting C2HX matrix elements...")
            out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2HX')
            temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2HX','temp')
            h5_polsar(
                h5_file=inFile,
                dataset_paths={
                    "HH": f"{base_path}/HH",
                    "HV": f"{base_path}/HV",
                },
                output_dir=out_dir, temp_dir=temp_dir,
                azlks=azlks, rglks=rglks, matrix_type = 'C2HX', apply_multilook=True,
                chunk_size_x=get_ml_chunk(rglks, 512), chunk_size_y=get_ml_chunk(azlks, 512), max_workers=max_workers,
                start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
                outType=outType, dtype = np.float32, inshape=inshape,outshape=outshape            
            )

        elif matrixType=='C2VX':
            print("Extracting C2VX matrix elements...")
            out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2VX')
            temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'C2VX','temp')
            h5_polsar(
                h5_file=inFile,
                dataset_paths={
                    "VV": f"{base_path}/VV",
                    "VH": f"{base_path}/VH",
                },
                output_dir=out_dir, temp_dir=temp_dir,
                azlks=azlks, rglks=rglks, matrix_type = 'C2VX', apply_multilook=True,
                chunk_size_x=get_ml_chunk(rglks, 512), chunk_size_y=get_ml_chunk(azlks, 512), max_workers=max_workers,
                start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
                outType=outType, dtype = np.float32, inshape=inshape,outshape=outshape            
            )

        elif matrixType=='T2HV':
            print("Extracting T2HV matrix elements...")
            out_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T2HV')
            temp_dir = os.path.join(inFolder,os.path.basename(inFile).split('.h5')[0],'T2HV','temp')
            h5_polsar(
                h5_file=inFile,
                dataset_paths={
                    "HH": f"{base_path}/HH",
                    "VV": f"{base_path}/VV",
                },
                output_dir=out_dir, temp_dir=temp_dir,
                azlks=azlks, rglks=rglks, matrix_type = 'T2HV', apply_multilook=True,
                chunk_size_x=get_ml_chunk(rglks, 512), chunk_size_y=get_ml_chunk(azlks, 512), max_workers=max_workers,
                start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
                outType=outType, dtype = np.float32, inshape=inshape,outshape=outshape            
            )
            
        else:
            raise ValueError(f"Unsupported matrix type: {matrixType} please choose from S2, C4, C3, T3, T4, T2HV, C2HV, C2HX, C2VX")



@time_it  
def nisar_gslc_old(inFile,azlks=22,rglks=10,matrixType='C3'):
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
        

def nisar_rslc(inFile,azlks=22,rglks=10, matrixType='C3'):
    return 0 
@time_it  
def nisar_rslc_old(inFile,azlks=22,rglks=10, matrixType='C3'):
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
        C11 = mlook_arr(C11,azlks,rglks)
        rows,cols = C11.shape
        write_rslc_bin( os.path.join(C2Folder,'C11.bin'), C11)
        print(f"Saved file {C2Folder}/C11.bin")
        del C11
        
        C22 = mlook_arr(C22,azlks,rglks)   
        write_rslc_bin( os.path.join(C2Folder,'C22.bin'), C22)
        print(f"Saved file {C2Folder}/C22.bin")
        del C22
        
        C12 = mlook_arr(C12,azlks,rglks)   
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


# def nisar_rslc_tab(inFile, azlks=22, rglks=10, matrixType='C3'):
#     inFolder = os.path.dirname(inFile)
#     if not inFolder:
#         inFolder = "."
    
#     C2Folder = os.path.join(inFolder, os.path.basename(inFile).split('.h5')[0], 'C2')

#     # Open the HDF5 file using PyTables with read-only access
#     with tables.open_file(inFile, "r") as h5File:
#         if '/science/LSAR' in h5File:
#             freq_band = 'L'
#         elif '/science/SSAR' in h5File:
#             freq_band = 'S'
#         else:
#             print("Neither LSAR nor SSAR data found in the file.")
#             return
        
#         # Base path logic (RSLC or SLC)
#         if f'/science/{freq_band}SAR/RSLC' in h5File:
#             base_path = f'/science/{freq_band}SAR/RSLC/swaths/frequencyA'
#         elif f'/science/{freq_band}SAR/SLC' in h5File:
#             base_path = f'/science/{freq_band}SAR/SLC/swaths/frequencyA'
#         else:
#             print("Neither RSLC nor SLC found in the HDF5 structure.")
#             return
        
#         listOfPolarizations = np.array(h5File.getNode(base_path + '/listOfPolarizations')).astype(str)
#         nchannels = len(listOfPolarizations)

#         print(f"Detected {freq_band}-band {listOfPolarizations}")

#         if nchannels == 2:
#             print("Extracting C2 matrix elements...")
#             if 'HH' in listOfPolarizations and 'HV' in listOfPolarizations:
#                 S11 = h5File.getNode(f'{base_path}/HH').read()  # Memory-map the HH dataset
#                 S12 = h5File.getNode(f'{base_path}/HV').read()  # Memory-map the HV dataset
#             elif 'VV' in listOfPolarizations and 'VH' in listOfPolarizations:
#                 S11 = h5File.getNode(f'{base_path}/VV').read()  # Memory-map the VV dataset
#                 S12 = h5File.getNode(f'{base_path}/VH').read()  # Memory-map the VH dataset
#             elif 'HH' in listOfPolarizations and 'VV' in listOfPolarizations:
#                 S11 = h5File.getNode(f'{base_path}/HH').read()  # Memory-map the HH dataset
#                 S12 = h5File.getNode(f'{base_path}/VV').read()  # Memory-map the VV dataset
#             else:
#                 print("No HH, HV, VV, or VH polarizations found in the file.")
#                 return

#             # Memory-mapped arrays
#             # S11 = np.array(S11)  # Convert to NumPy array (mmap)
#             # S12 = np.array(S12)  # Convert to NumPy array (mmap)

#             C11 = np.abs(S11)**2
#             C22 = np.abs(S12)**2
#             C12 = S11 * np.conjugate(S12)

#             del S11, S12

#             os.makedirs(C2Folder, exist_ok=True)

#             C11 = mlook_arr(C11, azlks, rglks)
#             rows, cols = C11.shape
#             write_rslc_bin(os.path.join(C2Folder, 'C11.bin'), C11)
#             print(f"Saved file {C2Folder}/C11.bin")
#             del C11

#             C22 = mlook_arr(C22, azlks, rglks)
#             write_rslc_bin(os.path.join(C2Folder, 'C22.bin'), C22)
#             print(f"Saved file {C2Folder}/C22.bin")
#             del C22

#             C12 = mlook_arr(C12, azlks, rglks)
#             write_rslc_bin(os.path.join(C2Folder, 'C12_real.bin'), np.real(C12))
#             print(f"Saved file {C2Folder}/C12_real.bin")
#             write_rslc_bin(os.path.join(C2Folder, 'C12_imag.bin'), np.imag(C12))
#             print(f"Saved file {C2Folder}/C12_imag.bin")
#             del C12

#             # Save config
#             with open(C2Folder + '/config.txt', "w+") as file:
#                 file.write(f'Nrow\n{rows}\n---------\nNcol\n{cols}\n---------\nPolarCase\nmonostatic\n---------\nPolarType\npp1')

#         elif nchannels == 4:
#             # Memory map the S matrix elements
#             S11 = h5File.getNode(f'{base_path}/HH').read()  # Memory-map the HH dataset
#             S12 = h5File.getNode(f'{base_path}/HV').read()  # Memory-map the HV dataset
#             S21 = h5File.getNode(f'{base_path}/VH').read()  # Memory-map the VH dataset
#             S22 = h5File.getNode(f'{base_path}/VV').read()  # Memory-map the VV dataset

#             if matrixType == 'S2':
#                 print("Extracting S2 matrix elements...")
#                 rows, cols = S11.shape
#                 out_dir = os.path.join(inFolder, os.path.basename(inFile).split('.h5')[0], 'S2')
#                 os.makedirs(out_dir, exist_ok=True)

#                 write_s2_bin_ref(os.path.join(out_dir, 's11.bin'), S11)
#                 write_s2_bin_ref(os.path.join(out_dir, 's12.bin'), (S12 + S21) / 2)
#                 write_s2_bin_ref(os.path.join(out_dir, 's21.bin'), (S12 + S21) / 2)
#                 write_s2_bin_ref(os.path.join(out_dir, 's22.bin'), S22)

#                 # Save config
#                 with open(os.path.join(out_dir, 'config.txt'), "w") as file:
#                     file.write(f'Nrow\n{rows}\n---------\nNcol\n{cols}\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull')

#             elif matrixType in ["C3", "T3", "C4", "T4"]:
#                 folder_name = matrixType
#                 print(f"Extracting {matrixType} matrix elements...")

#                 matrixFolder = os.path.join(inFolder, os.path.basename(inFile).split('.h5')[0], folder_name)
#                 os.makedirs(matrixFolder, exist_ok=True)

#                 # Processing matrix elements with memory mapping
#                 if matrixType == "C3":
#                     write_s2_ct_ref("C3", matrixFolder, [S11, np.sqrt(2) * 0.5 * (S12 + S21), S22], azlks, rglks)
#                 elif matrixType == "T3":
#                     write_s2_ct_ref("T3", matrixFolder, np.array([S11 + S22, S11 - S22, S12 + S21]) * (1 / np.sqrt(2)), azlks, rglks)
#                 elif matrixType == "C4":
#                     write_s2_ct_ref("C4", matrixFolder, [S11, S12, S21, S22], azlks, rglks)
#                 elif matrixType == "T4":
#                     write_s2_ct_ref("T4", matrixFolder, np.array([
#                                                                 S11 + S22, 
#                                                                 S11 - S22, 
#                                                                 S12 + S21, 
#                                                                 1j * (S12 - S21)
#                                                             ]) * (1 / np.sqrt(2)), azlks, rglks)

    

    
