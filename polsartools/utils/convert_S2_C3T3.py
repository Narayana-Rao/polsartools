

from polsartools.utils.utils import time_it
from polsartools.utils.io_utils import write_T3, write_C3,mlook,read_bin
from polsartools.utils.proc_utils import process_chunks_parallel
import numpy as np
import os
from osgeo import gdal


def get_s2_input_filepaths(infolder):
    """
    Searches for s11/s12/s21/s22 files in .bin or .tif format.
    Returns input file paths or raises an exception if files are missing.
    """
    def find_file(base):
        for ext in [".bin", ".tif"]:
            path = os.path.join(infolder, f"{base}{ext}")
            if os.path.isfile(path):
                return path
        return None

    keys = ["s11", "s12", "s21", "s22"]
    input_filepaths = [find_file(k) for k in keys]

    if not all(input_filepaths):
        raise FileNotFoundError("Invalid S2 folder: missing s11/s12/s21/s22 files")
    
    return input_filepaths


def get_output_filepaths(infolder, matrix, outType):
    """
    Returns output filepaths for the specified matrix and output type (bin or tif).
    Also ensures the target directory exists.
    """
    matrix_keys = {
        "T3": ["T11", "T12_real", "T12_imag", "T13_real", "T13_imag",
               "T22", "T23_real", "T23_imag", "T33"],
        "C3": ["C11", "C12_real", "C12_imag", "C13_real", "C13_imag",
               "C22", "C23_real", "C23_imag", "C33"],
        "C4": ["C11", "C12_real", "C12_imag", "C13_real", "C13_imag",
               "C14_real", "C14_imag", "C22", "C23_real", "C23_imag",
               "C24_real", "C24_imag", "C33", "C34_real", "C34_imag", "C44"],
        "T4": ["T11", "T12_real", "T12_imag", "T13_real", "T13_imag",
               "T14_real", "T14_imag", "T22", "T23_real", "T23_imag",
               "T24_real", "T24_imag", "T33", "T34_real", "T34_imag", "T44"]
    }

    if matrix not in matrix_keys:
        raise ValueError(f"Invalid matrix type '{matrix}'")

    ext = ".bin" if outType == "bin" else ".tif"
    outfolder = os.path.join(infolder, matrix)
    os.makedirs(outfolder, exist_ok=True)

    return [os.path.join(outfolder, f"{name}{ext}") for name in matrix_keys[matrix]]

@time_it
def convert_S2_CT(infolder, matrix='T3', azlks=8,rglks=2, cf = 1, 
                  outType="tif", cog_flag=False, cog_overviews = [2, 4, 8, 16], 
                  write_flag=True, max_workers=None,block_size=(512, 512)):
    """
    Converts Full-pol scattering matrix, S2 data into either the multi-looked T3 (coherency) or C3 (covariance) matrix format and saves them in PolSARpro format.

    Example Usage:
    --------------
    >>> convert_S2_CT('/path/to/S2_data', matrix='C3', azlks=10, rglks=5)

    Parameters:
    -----------
    infolder : str
        Path to the input folder scattering matrix data.

    matrix : str, optional, default='T3'
        The matrix type to generate. Options:
        - 'T3' : 3x3 Coherency matrix
        - 'C3' : 3x3 Covariance matrix
        - 'T4' : 4x4 Coherency matrix
        - 'C4' : 4x4 Covariance matrix

    azlks : int, optional, default=8
        Number of azimuth looks for multi-looking (averaging in the azimuth direction).

    rglks : int, optional, default=2
        Number of range looks for multi-looking (averaging in the range direction).

    max_workers : int, optional, default=None
        Number of workers for parallel processing.

    block_size : tuple of (int, int), optional, default=(512, 512)
        Block size for chunk-based parallel processing.

    Returns:
    --------
    None
        The function processes raster data and writes outputs to the specified folder.

    Raises:
    -------
    FileNotFoundError
        If the required Sentinel-2 input files are not found.
    
    Exception
        If an invalid matrix type is provided.

   
    """
    
    window_size=None
    
    input_filepaths =  get_s2_input_filepaths(infolder)
    output_filepaths = get_output_filepaths(infolder, matrix, outType)
  
    """
    GET MULTI-LOOKED RASTER PROPERTIES
       
    """
    
    dataset = gdal.Open(input_filepaths[0], gdal.GA_ReadOnly)
    if dataset is None:
        raise FileNotFoundError(f"Cannot open {input_filepaths[0]}")

    in_cols = dataset.RasterXSize
    in_rows = dataset.RasterYSize
    in_geotransform = dataset.GetGeoTransform()
    in_projection = dataset.GetProjection()

    # Calculate output size after multilooking
    out_x_size = in_cols // rglks
    out_y_size = in_rows // azlks

    # Calculate new geotransform with updated pixel size
    out_geotransform = list(in_geotransform)
    
    if in_geotransform[0]==0.0 and in_geotransform[1]==1.0 and in_geotransform[2]==-0.0 and in_geotransform[3]==0.0 and in_geotransform[4]==-0.0 and in_geotransform[5]==-1.0:
        out_geotransform[1] *= 1 
        out_geotransform[5] *= 1 
    else:
        out_geotransform[1] *= rglks 
        out_geotransform[5] *= azlks
        
    out_geotransform = tuple(out_geotransform)  # back to immutable
    # print(in_geotransform,out_geotransform)
    dataset = None  # Close GDAL dataset
    


    def closest_multiple(value, base):
        return round(value / base) * base

    # Calculate closest multiples of rglks and azlks to 512
    block_x_size = closest_multiple(block_size[0] , rglks)
    block_y_size = closest_multiple(block_size[1] , azlks)

    # print(block_x_size,block_y_size)
    block_size = (block_x_size, block_y_size)
    
    process_chunks_parallel(input_filepaths, list(output_filepaths), 
                             window_size,
                            write_flag,
                            process_chunk_s2ct,
                            *[cf, matrix, azlks, rglks],
                            block_size=block_size, 
                            max_workers=max_workers,  
                            num_outputs=len(output_filepaths),
                            cog_flag=cog_flag,
                            cog_overviews=cog_overviews,
                            out_x_size=out_x_size,
                            out_y_size=out_y_size,
                            out_geotransform=out_geotransform,
                            out_projection=in_projection,
                            azlks=azlks,
                            rglks=rglks
                            )

def process_chunk_s2ct(chunks, *args, **kwargs):
    # print(args[-1],args[-2],args[-3])
    abs_cf = args[-4]
    matrix=args[-3]
    azlks=args[-2]
    rglks=args[-1]
    
    s11 = np.array(chunks[0])*abs_cf
    s12 = np.array(chunks[1])*abs_cf
    s21 = np.array(chunks[2])*abs_cf
    s22 = np.array(chunks[3])*abs_cf
    
    if matrix == 'C4':
        Kl = np.array([s11, s12, s21, s22])

        del s11,s12,s21,s22
        
        C11 = mlook(np.abs(Kl[0])**2,azlks,rglks).astype(np.float32)
        C22 = mlook(np.abs(Kl[1])**2,azlks,rglks).astype(np.float32)
        C33 = mlook(np.abs(Kl[2])**2,azlks,rglks).astype(np.float32)
        C44 = mlook(np.abs(Kl[3])**2,azlks,rglks).astype(np.float32)

        C12 = mlook(Kl[0]*np.conj(Kl[1]),azlks,rglks).astype(np.complex64)
        C13 = mlook(Kl[0]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
        C14 = mlook(Kl[0]*np.conj(Kl[3]),azlks,rglks).astype(np.complex64)
        C23 = mlook(Kl[1]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
        C24 = mlook(Kl[1]*np.conj(Kl[3]),azlks,rglks).astype(np.complex64)
        C34 = mlook(Kl[2]*np.conj(Kl[3]),azlks,rglks).astype(np.complex64)
        
        del Kl
        return np.real(C11),np.real(C12),np.imag(C12),np.real(C13),np.imag(C13),np.real(C14),np.imag(C14), np.real(C22),np.real(C23),np.imag(C23),np.real(C24),np.imag(C24), np.real(C33),np.real(C34),np.imag(C34), np.real(C44)
        
        
    elif matrix == "T4":
        Kp = (1/np.sqrt(2))*np.array([s11+s22, s11-s22, s12+s21, 1j*(s12-s21)])

        del s11,s12,s21,s22

        # 4x4 Pauli Coherency Matrix elements
        T11 = mlook(np.abs(Kp[0])**2,azlks,rglks).astype(np.float32)
        T22 = mlook(np.abs(Kp[1])**2,azlks,rglks).astype(np.float32)
        T33 = mlook(np.abs(Kp[2])**2,azlks,rglks).astype(np.float32)
        T44 = mlook(np.abs(Kp[3])**2,azlks,rglks).astype(np.float32)

        T12 = mlook(Kp[0]*np.conj(Kp[1]),azlks,rglks).astype(np.complex64)
        T13 = mlook(Kp[0]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)
        T14 = mlook(Kp[0]*np.conj(Kp[3]),azlks,rglks).astype(np.complex64)
        T23 = mlook(Kp[1]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)
        T24 = mlook(Kp[1]*np.conj(Kp[3]),azlks,rglks).astype(np.complex64)
        T34 = mlook(Kp[2]*np.conj(Kp[3]),azlks,rglks).astype(np.complex64)  

        del Kp
        return np.real(T11),np.real(T12),np.imag(T12),np.real(T13),np.imag(T13), np.real(T14),np.imag(T14),np.real(T22),np.real(T23),np.imag(T23), np.real(T24),np.imag(T24),np.real(T33),np.real(T34),np.imag(T34),np.real(T44)    


    elif matrix == "T3":
        Kp = (1/np.sqrt(2))*np.array([s11+s22, s11-s22, s12+s21])

        del s11,s12,s21,s22

        # 3x3 Pauli Coherency Matrix elements
        T11 = mlook(np.abs(Kp[0])**2,azlks,rglks).astype(np.float32)
        T22 = mlook(np.abs(Kp[1])**2,azlks,rglks).astype(np.float32)
        T33 = mlook(np.abs(Kp[2])**2,azlks,rglks).astype(np.float32)

        T12 = mlook(Kp[0]*np.conj(Kp[1]),azlks,rglks).astype(np.complex64)
        T13 = mlook(Kp[0]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)
        T23 = mlook(Kp[1]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)

        del Kp
        return np.real(T11),np.real(T12),np.imag(T12),np.real(T13),np.imag(T13),np.real(T22),np.real(T23),np.imag(T23),np.real(T33)

        
    elif matrix=='C3':
        # Kl- 3-D Lexicographic feature vector
        Kl = np.array([s11, np.sqrt(2)*0.5*(s12+s21), s22])
        del s11,s12,s22

        # 3x3 COVARIANCE Matrix elements

        C11 = mlook(np.abs(Kl[0])**2,azlks,rglks).astype(np.float32)
        C22 = mlook(np.abs(Kl[1])**2,azlks,rglks).astype(np.float32)
        C33 = mlook(np.abs(Kl[2])**2,azlks,rglks).astype(np.float32)

        C12 = mlook(Kl[0]*np.conj(Kl[1]),azlks,rglks).astype(np.complex64)
        C13 = mlook(Kl[0]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
        C23 = mlook(Kl[1]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
        return np.real(C11),np.real(C12),np.imag(C12),np.real(C13),np.imag(C13),np.real(C22),np.real(C23),np.imag(C23),np.real(C33)

