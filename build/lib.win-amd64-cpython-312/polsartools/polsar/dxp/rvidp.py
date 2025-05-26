import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it
from .dxp_infiles import dxpc2files
@time_it
def rvidp(infolder,  window_size=1, outType="tif", cog_flag=False, cog_overviews = [2, 4, 8, 16], write_flag=True, max_workers=None,block_size=(512, 512)):
    """
        
        Computes Radar vegetation index from the input dual-polarization (dual-pol) C2 matrix data, and writes
        the output in various formats (GeoTIFF or binary). The computation is performed in parallel for efficiency.

        Example:
        --------
        >>> rvidp("path_to_C2_folder", window_size=5, outType="tif", cog_flag=True)
        This will compute Radar vegetation index from the C2 matrix in the specified folder,
        generating output in Geotiff format with Cloud Optimized GeoTIFF settings enabled.
        
        Parameters:
        -----------
        infolder : str
            Path to the input folder containing C2 matrix data.
        window_size : int, optional
            Size of the processing window (default is 1).
        outType : str, optional
            Output format of the files; can be "tif" (GeoTIFF) or "bin" (binary) (default is "tif").
        cog_flag : bool, optional
            If True, outputs Cloud Optimized GeoTIFF (COG) (default is False).
        cog_overviews : list of int, optional
            List of overview levels to be used for COGs (default is [2, 4, 8, 16]).
        write_flag : bool, optional
            Whether to write the computed output files (default is True).
        max_workers : int, optional
            Number of parallel workers for processing (default is None, which uses one less than the number of available CPU cores).
        block_size : tuple of int, optional
            Size of each processing block (default is (512, 512)), defining the spatial chunk dimensions used in parallel computation.
    
        Returns:
        --------
        None
            The function writes the computed Radar vegetation index to the specified output format.

        Output Files:
        -------------
        - "rvidp.tif" or "rvidp.bin"

    """
    input_filepaths = dxpc2files(infolder)
    output_filepaths = []
    if outType == "bin":
        output_filepaths.append(os.path.join(infolder, "rvidp.bin"))
    else:
        output_filepaths.append(os.path.join(infolder, "rvidp.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_rvidp,block_size=block_size, max_workers=max_workers,  num_outputs=1,
                            cog_flag=cog_flag,
                            cog_overviews=cog_overviews,
                            )


def process_chunk_rvidp(chunks, window_size,*args):
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    if window_size>1:
        c11_T1r = conv2d(np.real(c11_T1),kernel)
        c11_T1i = conv2d(np.imag(c11_T1),kernel)
        c11s = c11_T1r+1j*c11_T1i

        c12_T1r = conv2d(np.real(c12_T1),kernel)
        c12_T1i = conv2d(np.imag(c12_T1),kernel)
        c12s = c12_T1r+1j*c12_T1i


        c21_T1r = conv2d(np.real(c21_T1),kernel)
        c21_T1i = conv2d(np.imag(c21_T1),kernel)
        c21s = c21_T1r+1j*c21_T1i


        c22_T1r = conv2d(np.real(c22_T1),kernel)
        c22_T1i = conv2d(np.imag(c22_T1),kernel)
        c22s = c22_T1r+1j*c22_T1i

        
        c2_det = (c11s*c22s-c12s*c21s)
        c2_trace = c11s+c22s
        rvidp = np.real(4*c22s/c2_trace)
    else:
        c2_det = (c11_T1*c22_T1-c12_T1*c21_T1)
        c2_trace = c11_T1+c22_T1
        rvidp = np.real(4*c22_T1/c2_trace)

    return rvidp