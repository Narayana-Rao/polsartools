import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it,eig22
from .dxp_infiles import dxpc2files
@time_it
def dprvi(infolder,  window_size=1, outType="tif", cog_flag=False, cog_overviews = [2, 4, 8, 16], write_flag=True, max_workers=None,block_size=(512, 512)):
    """
        
        Computes  dual-pol Radar vegetation index from the input dual-polarization (dual-pol) C2 matrix data, and writes
        the output in various formats (GeoTIFF or binary). The computation is performed in parallel for efficiency.

        Example:
        --------
        >>> dprvi("path_to_C2_folder", window_size=5, outType="tif", cog_flag=True)
        This will compute dual-pol Radar vegetation index from the C2 matrix in the specified folder,
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
            The function writes the computed dual-pol  Radar vegetation index to the specified output format.

        Output Files:
        -------------
        - "dprvi.tif" or "dprvi.bin"

    """
    input_filepaths = dxpc2files(infolder)
    output_filepaths = []

    if outType == "bin":
        output_filepaths.append(os.path.join(infolder, "dprvi.bin"))
    else:
        output_filepaths.append(os.path.join(infolder, "dprvi.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_dprvi,block_size=block_size, max_workers=max_workers,  num_outputs=1,
                            cog_flag=cog_flag,
                            cog_overviews=cog_overviews,
                            )
    
def process_chunk_dprvi(chunks, window_size,*args):
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    if window_size>1:
        c11s = conv2d(np.real(c11_T1),kernel)+1j*conv2d(np.imag(c11_T1),kernel)
        c12s = conv2d(np.real(c12_T1),kernel)+1j*conv2d(np.imag(c12_T1),kernel)
        c21s = conv2d(np.real(c21_T1),kernel)+1j*conv2d(np.imag(c21_T1),kernel)
        c22s = conv2d(np.real(c22_T1),kernel)+1j*conv2d(np.imag(c22_T1),kernel)

        c2_det = (c11s*c22s-c12s*c21s)
        c2_trace = c11s+c22s
        # t2_span = t11s*t22s
        m = (np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))
        egv1,egv2 = eig22(np.dstack([c11s,c12s,c21s,c22s]))
        egf = np.vstack([egv1,egv2])
        egfmax = egf.max(axis=0)#.reshape(np.shape(C2_stack[:,:,0]))
        beta = (egfmax/(egv1+egv2)).reshape(np.shape(c11s))
        dprvi = np.real(1-(m*beta))
    else:
        c2_det = (c11_T1*c22_T1-c12_T1*c21_T1)
        c2_trace = c11_T1+c22_T1
        m = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))

        egv1,egv2 = eig22(np.dstack([c11_T1,c12_T1,c21_T1,c22_T1]))
        egf = np.vstack([egv1,egv2])
        egfmax = egf.max(axis=0)#.reshape(np.shape(C2_stack[:,:,0]))
        beta = (egfmax/(egv1+egv2)).reshape(np.shape(c11_T1))
        dprvi = np.real(1-(m*beta))


    return dprvi