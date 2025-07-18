import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it,eig22
from .dxp_infiles import dxpc2files
@time_it
def dp_desc(cpFile, xpFile,  window_size=1, outType="tif", cog_flag=False, cog_overviews = [2, 4, 8, 16], write_flag=True, max_workers=None,block_size=(512, 512)):
    """Compute dual polarimetric descriptors (mc, Hc, and Thetac) from Dual-pol GRD data.

    This function calculates the mc, Degree of Purity), Hc, Pseudo-Entropy, and Thetac, Pseudo-Scattering Type Parameter 
    using co-polarized (`cpFile`) and cross-polarized (`xpFile`) SAR raster files. 

    Examples
    --------
    >>> # Basic usage with default parameters
    >>> dp_desc("/path/to/copol_file.tif", "/path/to/crosspol_file.tif")

    >>> # Advanced usage with custom parameters
    >>> dp_desc(
    ...     cpFile="/path/to/copol_file.tif",
    ...     xpFile="/path/to/crosspol_file.tif",
    ...     window_size=3,
    ...     outType="tif",
    ...     cog_flag=True,
    ...     block_size=(1024, 1024)
    ... )
    
    Parameters
    ----------
    cpFile : str
        Path to the co-polarized backscatter (linear) SAR raster file.
    xpFile : str
        Path to the cross-polarized backscatter (linear) SAR raster file.
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
        Results are saved as 'mc.tif', 'Hc.tif', and 'Thetac.tif' (or '.bin' if `outType='bin'`).

    """
    input_filepaths = [cpFile,xpFile]
    output_filepaths = []

    if outType == "bin":
        output_filepaths.append(os.path.join(os.path.dirname(cpFile), "mc.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(cpFile), "Hc.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(cpFile), "Thetac.bin"))
    else:
        output_filepaths.append(os.path.join(os.path.dirname(cpFile), "mc.tif"))
        output_filepaths.append(os.path.join(os.path.dirname(cpFile), "Hc.tif"))
        output_filepaths.append(os.path.join(os.path.dirname(cpFile), "Thetac.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_dpdesc,block_size=block_size, max_workers=max_workers,  num_outputs=len(output_filepaths),
                            cog_flag=cog_flag,
                            cog_overviews=cog_overviews,
                            )
    
def process_chunk_dpdesc(chunks, window_size,*args):
    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11 = np.array(chunks[0])
    c22 = np.array(chunks[1])


    if window_size>1:
        c11s = conv2d(c11,kernel)
        c22s = conv2d(c22,kernel)

        q = c22s/c11s
        q[q>=1]=1
        mc = (1-q)/(1+q)
        p1 = 1/(1+q)
        p2 = q/(1+q)
        Hc = -1*(p1*np.log2(p1)+p2*np.log2(p2))
        thetac = np.arctan(((1-q)**2)/(1-q+q**2)) * (180/np.pi)
        
    else:
        q = c22/c11
        q[q>=1]=1
        mc = (1-q)/(1+q)
        p1 = 1/(1+q)
        p2 = q/(1+q)
        Hc = -1*(p1*np.log2(p1)+p2*np.log2(p2))
        thetac = np.arctan(((1-q)**2)/(1-q+q**2)) * (180/np.pi)
        

    return mc.astype(np.float32),Hc.astype(np.float32),thetac.astype(np.float32)