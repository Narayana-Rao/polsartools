import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it
from .cp_infiles import cpc2files
@time_it
def dopcp(infolder,   chi_in=45, psi_in=0, window_size=1, outType="tif", cog_flag=False, 
          cog_overviews = [2, 4, 8, 16], write_flag=True, 
          max_workers=None,block_size=(512, 512)):
    """Compute Degree of Polarization (DoP) from compact-pol SAR data.

        This function calculates the Degree of Polarization (DoP) from compact-polarimetric
        SAR data, which quantifies the polarized portion of the scattered wave. DoP is a
        useful parameter for characterizing surface properties and scattering mechanisms.

        Examples
        --------
        >>> # Basic usage with default parameters (right circular transmission)
        >>> dopcp("/path/to/cp_data")
        
        >>> # Advanced usage with custom parameters
        >>> dopcp(
        ...     infolder="/path/to/cp_data",
        ...     chi_in=-45,  # left circular transmission
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
            Size of the spatial averaging window. Larger windows provide better
            estimation of DoP but reduce spatial resolution.
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
            Results are written to disk as either 'dopcp.tif' or 'dopcp.bin'
            in the input folder. The output DoP values range from 0 to 1, where:
            - 0 indicates completely unpolarized scattered wave
            - 1 indicates completely polarized scattered wave

        Notes
        -----
        The Degree of Polarization (DoP) is calculated using the Stokes parameters
        derived from the compact-pol coherency matrix. For a partially polarized wave,
        DoP is given by:

        DoP = sqrt(S₁² + S₂² + S₃²) / S₀

        where S₀, S₁, S₂, S₃ are the Stokes parameters.

        Key characteristics:
        - DoP is invariant to the wave polarization basis
        - Values are normalized between 0 and 1
        - Higher values indicate stronger polarized scattering
        - Lower values suggest depolarizing mechanisms

        Common applications:
        - Surface roughness estimation
        - Vegetation density analysis
        - Urban area characterization
        - Sea ice classification


        """
    
    input_filepaths = cpc2files(infolder)

    output_filepaths = []
    if outType == "bin":
        output_filepaths.append(os.path.join(infolder, "dopcp.bin"))
    else:   
        output_filepaths.append(os.path.join(infolder, "dopcp.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), 
                            window_size,
                        write_flag,
                        process_chunk_dopcp,
                        *[chi_in, psi_in],
                        block_size=block_size, 
                        max_workers=max_workers,  
                        num_outputs=len(output_filepaths),
                        cog_flag=cog_flag,
                        cog_overviews=cog_overviews,

                        )

def process_chunk_dopcp(chunks, window_size, *args, **kwargs):
    
    chi_in=args[-2]
    psi_in=args[-1]
    
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
        
    # Compute Stokes parameters
    s0 = c11_T1 + c22_T1
    s1 = c11_T1 - c22_T1
    s2 = c12_T1 + c21_T1
    s3 = np.where(chi_in >= 0, 1j * (c12_T1 - c21_T1), -1j * (c12_T1 - c21_T1))

    dop= np.sqrt(np.power(s1,2) + np.power(s2,2) + np.power(s3,2))/(s0);   

    return np.real(dop).astype(np.float32)