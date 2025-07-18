import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it,eig22
from .dxp_infiles import dxpc2files
@time_it
def dprvi(infolder,  window_size=1, outType="tif", cog_flag=False, cog_overviews = [2, 4, 8, 16], write_flag=True, max_workers=None,block_size=(512, 512)):
    """Compute dual-pol Radar Vegetation Index (DpRVI) from C2 matrix data.

    This function processes dual-polarization SAR data to generate the DpRVI, which is useful
    for vegetation monitoring and biomass estimation. The processing is done in parallel
    blocks for improved performance.

    Examples
    --------
    >>> # Basic usage with default parameters
    >>> dprvi("/path/to/c2_data")
    
    >>> # Advanced usage with custom parameters
    >>> dprvi(
    ...     infolder="/path/to/c2_data",
    ...     window_size=3,
    ...     outType="tif",
    ...     cog_flag=True,
    ...     block_size=(1024, 1024)
    ... )
    
    Parameters
    ----------
    infolder : str
        Path to the input folder containing C2 matrix files.
    window_size : int, default=1
        Size of the spatial averaging window. Larger windows reduce speckle noise
        but decrease spatial resolution.
    outType : {'tif', 'bin'}, default='tif'
        Output file format:
        - 'tif': GeoTIFF format with georeferencing information
        - 'bin': Raw binary format
    cog_flag : bool, default=False
        If True, creates a Cloud Optimized GeoTIFF (COG) with internal tiling
        and overviews for efficient web access.
    cog_overviews : list[int], default=[2, 4, 8, 16]
        Overview levels for COG creation. Each number represents the
        decimation factor for that overview level.
    write_flag : bool, default=True
        If True, writes results to disk. If False, only processes data in memory.
    max_workers : int | None, default=None
        Maximum number of parallel processing workers. If None, uses
        CPU count - 1 workers.
    block_size : tuple[int, int], default=(512, 512)
        Size of processing blocks (rows, cols) for parallel computation.
        Larger blocks use more memory but may be more efficient.

    Returns
    -------
    None
        Results are written to disk as either 'dprvi.tif' or 'dprvi.bin'
        in the input folder.

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


    return dprvi.astype(np.float32)