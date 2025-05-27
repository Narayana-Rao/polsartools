

from polsartools.utils.utils import time_it
from polsartools.utils.io_utils import write_T3, write_C3,mlook,read_bin
from polsartools.utils.proc_utils import process_chunks_parallel
import numpy as np
import os
from osgeo import gdal



@time_it
def convert_S2_CT(infolder, matrix='T3', azlks=8,rglks=2, window_size=None, 
                  outType="bin", cog_flag=False, cog_overviews = [2, 4, 8, 16], 
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
        - 'T3' : Coherency matrix
        - 'C3' : Covariance matrix

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
    
    if os.path.isfile(os.path.join(infolder,"s11.bin")):
        input_filepaths = [
        os.path.join(infolder, "s11.bin"), 
        os.path.join(infolder, "s12.bin"),
        os.path.join(infolder, "s21.bin"),
        os.path.join(infolder, "s22.bin")
    ]
    else:
        print(f"Invalid S2 folder!!")
    
    output_filepaths = []
    num_outputs = 1
    if matrix == "T3":
        os.makedirs(os.path.join(infolder,"T3"), exist_ok=True)
        output_filepaths.append(os.path.join(infolder,"T3", "T11.bin"))
        output_filepaths.append(os.path.join(infolder,"T3","T12_real.bin"))
        output_filepaths.append(os.path.join(infolder,"T3","T12_imag.bin"))
        output_filepaths.append(os.path.join(infolder,"T3","T13_real.bin"))
        output_filepaths.append(os.path.join(infolder,"T3","T13_imag.bin"))
        output_filepaths.append(os.path.join(infolder,"T3","T22.bin"))
        output_filepaths.append(os.path.join(infolder,"T3","T23_real.bin"))
        output_filepaths.append(os.path.join(infolder,"T3","T23_imag.bin"))
        output_filepaths.append(os.path.join(infolder,"T3","T33.bin"))
    elif matrix == "C3":
        os.makedirs(os.path.join(infolder,"C3"), exist_ok=True)
        output_filepaths.append(os.path.join(infolder,"C3", "C11.bin"))
        output_filepaths.append(os.path.join(infolder,"C3","C12_real.bin"))
        output_filepaths.append(os.path.join(infolder,"C3","C12_imag.bin"))
        output_filepaths.append(os.path.join(infolder,"C3","C13_real.bin"))
        output_filepaths.append(os.path.join(infolder,"C3","C13_imag.bin"))
        output_filepaths.append(os.path.join(infolder,"C3","C22.bin"))
        output_filepaths.append(os.path.join(infolder,"C3","C23_real.bin"))
        output_filepaths.append(os.path.join(infolder,"C3","C23_imag.bin"))
        output_filepaths.append(os.path.join(infolder,"C3","C33.bin"))
    else:
        raise Exception(f"Invalid matrix type!! Available types: ['T3', 'C3']")
    
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
                            *[matrix, azlks, rglks],
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
    matrix=args[-3]
    azlks=args[-2]
    rglks=args[-1]
    
    s11 = np.array(chunks[0])
    s12 = np.array(chunks[1])
    s22 = np.array(chunks[3])
    
    if matrix == "T3":
        Kp = (1/np.sqrt(2))*np.array([s11+s22, s11-s22, 2*s12])

        del s11,s12,s22

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
        Kl = np.array([s11, np.sqrt(2)*s12, s22])
        del s11,s12,s22

        # 3x3 COVARIANCE Matrix elements

        C11 = mlook(np.abs(Kl[0])**2,azlks,rglks).astype(np.float32)
        C22 = mlook(np.abs(Kl[1])**2,azlks,rglks).astype(np.float32)
        C33 = mlook(np.abs(Kl[2])**2,azlks,rglks).astype(np.float32)

        C12 = mlook(Kl[0]*np.conj(Kl[1]),azlks,rglks).astype(np.complex64)
        C13 = mlook(Kl[0]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
        C23 = mlook(Kl[1]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
        return np.real(C11),np.real(C12),np.imag(C12),np.real(C13),np.imag(C13),np.real(C22),np.real(C23),np.imag(C23),np.real(C33)






def convert_S2_CT_old(inFolder, matrix='T3',azlks=8,rglks=2):

    if os.path.isfile(os.path.join(inFolder,"s11.bin")) and os.path.isfile(os.path.join(inFolder,"s12.bin")) and os.path.isfile(os.path.join(inFolder,"s22.bin")):
        s11 = read_bin(os.path.join(inFolder,"s11.bin"))
        s12 = read_bin(os.path.join(inFolder,"s12.bin"))
        s22 = read_bin(os.path.join(inFolder,"s22.bin"))

        if matrix == "T3":
            Kp = (1/np.sqrt(2))*np.array([s11+s22, s11-s22, 2*s12])

            del s11,s12,s22

            # 3x3 Pauli Coherency Matrix elements
            T11 = mlook(np.abs(Kp[0])**2,azlks,rglks).astype(np.float32)
            T22 = mlook(np.abs(Kp[1])**2,azlks,rglks).astype(np.float32)
            T33 = mlook(np.abs(Kp[2])**2,azlks,rglks).astype(np.float32)

            T12 = mlook(Kp[0]*np.conj(Kp[1]),azlks,rglks).astype(np.complex64)
            T13 = mlook(Kp[0]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)
            T23 = mlook(Kp[1]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)

            del Kp
            T3Folder = os.path.join(inFolder,'T3')

            if not os.path.isdir(T3Folder):
                print("T3 folder does not exist. \nCreating folder {}".format(T3Folder))
                os.mkdir(T3Folder)
                
            # write_T3(np.dstack([T11,T12,T13,np.conjugate(T12),T22,T23,np.conjugate(T13),np.conjugate(T23),T33]),T3Folder)
            write_T3([np.real(T11),np.real(T12),np.imag(T12),np.real(T13),np.imag(T13),
                    np.real(T22),np.real(T23),np.imag(T23),
                    np.real(T33)],T3Folder)
            
            
        elif matrix=='C3':
            # Kl- 3-D Lexicographic feature vector
            Kl = np.array([s11, np.sqrt(2)*s12, s22])
            del s11,s12,s22

            # 3x3 COVARIANCE Matrix elements

            C11 = mlook(np.abs(Kl[0])**2,azlks,rglks).astype(np.float32)
            C22 = mlook(np.abs(Kl[1])**2,azlks,rglks).astype(np.float32)
            C33 = mlook(np.abs(Kl[2])**2,azlks,rglks).astype(np.float32)

            C12 = mlook(Kl[0]*np.conj(Kl[1]),azlks,rglks).astype(np.complex64)
            C13 = mlook(Kl[0]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
            C23 = mlook(Kl[1]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)

            C3Folder = os.path.join(inFolder,'C3')

            if not os.path.isdir(C3Folder):
                print("C3 folder does not exist. \nCreating folder {}".format(C3Folder))
                os.mkdir(C3Folder)
            
            # write_C3(np.dstack([C11,C12,C13,np.conjugate(C12),C22,C23,np.conjugate(C13),np.conjugate(C23),C33]),C3Folder)
            write_C3([np.real(C11),np.real(C12),np.imag(C12),np.real(C13),np.imag(C13),
                    np.real(C22),np.real(C23),np.imag(C23),
                    np.real(C33)],C3Folder)
        else:
            print("Matrix type not supported. Supported types are 'T3' and 'C3'")
            
    else:
        print("s11, s12 and s22 files not found in {}".format(inFolder))
        return