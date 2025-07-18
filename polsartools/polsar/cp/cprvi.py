import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import time_it
from polsartools.cprvicpp import process_chunk_cprvicpp
from .cp_infiles import cpc2files
@time_it
def cprvi(infolder,   chi_in=45, psi_in=0, window_size=1, outType="tif", cog_flag=False, 
          cog_overviews = [2, 4, 8, 16], write_flag=True, 
          max_workers=None,block_size=(512, 512)):

    """Compute compact-pol Radar Vegetation Index (CpRVI) from C2 matrix data.

    This function processes compact-polarimetric SAR data to generate the CP-RVI, which
    is useful for vegetation monitoring and biomass estimation using compact-pol SAR systems.
    The processing is done in parallel blocks for improved performance.

    Examples
    --------
    >>> # Basic usage with default parameters (right circular transmission)
    >>> cprvi("/path/to/cp_data")
    
    >>> # Custom parameters for left circular transmission
    >>> cprvi(
    ...     infolder="/path/to/cp_data",
    ...     chi_in=-45,
    ...     psi_in=0,
    ...     window_size=3,
    ...     outType="tif",
    ...     cog_flag=True
    ... )
    
    Parameters
    ----------
    infolder : str
        Path to the input folder containing compact-pol C2 matrix files.
    chi_in : float, default=45
        Ellipticity angle chi of the transmitted wave in degrees.
        For circular polarization, chi = 45° (right circular) or -45° (left circular).
    psi_in : float, default=0
        Orientation angle psi of the transmitted wave in degrees.
        For circular polarization, typically 0°.
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
        Results are written to disk as either 'cprvi.tif' or 'cprvi.bin'
        in the input folder.

    """
    input_filepaths = cpc2files(infolder)
    output_filepaths = []
    if outType == "bin":
        output_filepaths.append(os.path.join(infolder, "cprvi.bin"))
    else:
        output_filepaths.append(os.path.join(infolder, "cprvi.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), 
                            window_size,
                        write_flag,
                        process_chunk_cprvi,
                        *[chi_in, psi_in],
                        block_size=block_size, 
                        max_workers=max_workers,  
                        num_outputs=len(output_filepaths),
                        cog_flag=cog_flag,
                        cog_overviews=cog_overviews,

                        )
def process_chunk_cprvi(chunks, window_size, *args, **kwargs):
    
    chi_in=args[-2]
    psi_in=args[-1]
    # print(chi_in,psi_in)
    
    chunk_arrays = [np.array(ch) for ch in chunks]  
    # CPP function
    vi_c_raw = process_chunk_cprvicpp( chunk_arrays, window_size, chi_in, psi_in )

    return np.array(vi_c_raw, copy=True).astype(np.float32) 
    
 