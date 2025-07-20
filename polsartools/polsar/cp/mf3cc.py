import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it
from .cp_infiles import cpc2files
@time_it
def mf3cc(infolder,   chi_in=45, psi_in=0, window_size=1, outType="tif", cog_flag=False, 
          cog_overviews = [2, 4, 8, 16], write_flag=True, 
          max_workers=None,block_size=(512, 512),
          progress_callback=None  # for QGIS plugin
          ):
    """Perform Model-Free 3-Component Decomposition for compact-pol SAR data.

    This function implements the model-free three-component decomposition for
    compact-polarimetric SAR data, decomposing the total backscattered power into
    surface (Ps), double-bounce (Pd), and volume (Pv) scattering components, along
    with the scattering-type parameter (Theta_CP).

    Examples
    --------
    >>> # Basic usage with default parameters
    >>> mf3cc("/path/to/cp_data")
    
    >>> # Advanced usage with custom parameters
    >>> mf3cc(
    ...     infolder="/path/to/cp_data",
    ...     chi_in=-45,
    ...     window_size=5,
    ...     outType="tif",
    ...     cog_flag=True,
    ...     block_size=(1024, 1024)
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
        If True, creates Cloud Optimized GeoTIFF (COG) outputs with internal tiling
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
        Writes four output files to disk:
        1. Ps_mf3cc: Surface scattering power component
        2. Pd_mf3cc: Double-bounce scattering power component
        3. Pv_mf3cc: Volume scattering power component
        4. Theta_CP_mf3cc: Scattering-type parameter

    """
    
    input_filepaths = cpc2files(infolder)

    output_filepaths = []
    if outType == "bin":
        output_filepaths.append(os.path.join(infolder, "Ps_mf3cc.bin"))
        output_filepaths.append(os.path.join(infolder, "Pd_mf3cc.bin"))
        output_filepaths.append(os.path.join(infolder, "Pv_mf3cc.bin"))
        output_filepaths.append(os.path.join(infolder, "Theta_CP_mf3cc.bin"))
    else:
        output_filepaths.append(os.path.join(infolder, "Ps_mf3cc.tif"))
        output_filepaths.append(os.path.join(infolder, "Pd_mf3cc.tif"))
        output_filepaths.append(os.path.join(infolder, "Pv_mf3cc.tif"))
        output_filepaths.append(os.path.join(infolder, "Theta_CP_mf3cc.tif"))
        
    process_chunks_parallel(input_filepaths, list(output_filepaths), 
                            window_size,
                        write_flag,
                        process_chunk_mf3cc,
                        *[chi_in, psi_in],
                        block_size=block_size, 
                        max_workers=max_workers,  
                        num_outputs=len(output_filepaths),
                        cog_flag=cog_flag,
                        cog_overviews=cog_overviews,
                        progress_callback=progress_callback
                        )
def process_chunk_mf3cc(chunks, window_size, *args, **kwargs):
    
    chi_in=args[-2]
    psi_in=args[-1]
    # print(chi_in,psi_in):

    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    ncols,nrows = np.shape(c11_T1)

    if window_size>1:
        c11_T1 = conv2d(np.real(c11_T1),kernel)+1j*conv2d(np.imag(c11_T1),kernel)
        c12_T1 = conv2d(np.real(c12_T1),kernel)+1j*conv2d(np.imag(c12_T1),kernel)
        c21_T1 = conv2d(np.real(c21_T1),kernel)+1j*conv2d(np.imag(c21_T1),kernel)
        c22_T1 = conv2d(np.real(c22_T1),kernel)+1j*conv2d(np.imag(c22_T1),kernel)
    
    c2_det = (c11_T1*c22_T1-c12_T1*c21_T1)
    c2_trace = c11_T1+c22_T1
    # t2_span = t11s*t22s
    m1 = np.real(np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))

    # Compute Stokes parameters
    s0 = c11_T1 + c22_T1
    s1 = c11_T1 - c22_T1
    s2 = np.real(c12_T1 + c21_T1)
    s3 = np.where(chi_in >= 0, 1j * (c12_T1 - c21_T1), -1j * (c12_T1 - c21_T1))
    s3 = np.real(s3)

    SC = ((s0)-(s3))/2;
    OC = ((s0)+(s3))/2;

    h = (OC-SC)
    # span = c11s + c22s

    val = ((m1*s0*h))/((SC*OC + (m1**2)*(s0**2)))
    thet = np.real(np.arctan(val))
    theta_CP = np.rad2deg(thet)

    Ps_CP= (((m1*(c2_trace)*(1.0+np.sin(2*thet))/2)))
    Pd_CP= (((m1*(c2_trace)*(1.0-np.sin(2*thet))/2)))
    Pv_CP= (c2_trace*(1.0-m1))
    

    return Ps_CP.astype(np.float32), Pd_CP.astype(np.float32), Pv_CP.astype(np.float32), theta_CP.astype(np.float32)