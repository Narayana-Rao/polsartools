import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it
from polsartools.utils.convert_matrices import T3_C3_mat, C3_T3_mat
from .fp_infiles import fp_c3t3files
@time_it

def yam4cfp(infolder,  model="", window_size=1, outType="tif", cog_flag=False, 
                    cog_overviews = [2, 4, 8, 16], write_flag=True, 
                    max_workers=None,block_size=(512, 512)):
    
    """Perform Yamaguchi 4-Component Decomposition for full-pol SAR data.

    This function implements the Yamaguchi 4-component decomposition with three
    different model options: original (Y4O), rotation-corrected (Y4R), and
    extended volume scattering model (Y4S). The decomposition separates the total
    power into surface, double-bounce, volume, and helix scattering components.

    Examples
    --------
    >>> # Original Yamaguchi decomposition
    >>> yam4cfp("/path/to/fullpol_data")
    
    >>> # Rotation-corrected decomposition
    >>> yam4cfp(
    ...     infolder="/path/to/fullpol_data",
    ...     model="y4cr",
    ...     window_size=5,
    ...     outType="tif",
    ...     cog_flag=True
    ... )
    
    >>> # Extended volume model decomposition
    >>> yam4cfp(
    ...     infolder="/path/to/fullpol_data",
    ...     model="y4cs",
    ...     window_size=5
    ... )

    Parameters
    ----------
    infolder : str
        Path to the input folder containing full-pol T3 or C3 matrix files.
    model : {'', 'y4cr', 'y4cs'}, default=''
        Decomposition model to use:
        - '': Original Yamaguchi 4-component (Y4O)
        - 'y4cr': Rotation-corrected Yamaguchi (Y4R)
        - 'y4cs': Extended volume scattering model (Y4S)
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
        Writes four output files to disk for the selected model:
        1. {prefix}_odd: Surface scattering power
        2. {prefix}_dbl: Double-bounce scattering power
        3. {prefix}_vol: Volume scattering power
        4. {prefix}_hlx: Helix scattering power
        where prefix is 'Yam4co', 'Yam4cr', or 'Yam4csr' based on model choice.


    """
    input_filepaths = fp_c3t3files(infolder)

    output_filepaths = []
    if outType == "bin":
        if not model or model=="y4co":
            output_filepaths.append(os.path.join(infolder, "Yam4co_odd.bin"))
            output_filepaths.append(os.path.join(infolder, "Yam4co_dbl.bin"))
            output_filepaths.append(os.path.join(infolder, "Yam4co_vol.bin"))
            output_filepaths.append(os.path.join(infolder, "Yam4co_hlx.bin"))
        elif model=="y4cr":
            output_filepaths.append(os.path.join(infolder, "Yam4cr_odd.bin"))
            output_filepaths.append(os.path.join(infolder, "Yam4cr_dbl.bin"))
            output_filepaths.append(os.path.join(infolder, "Yam4cr_vol.bin"))
            output_filepaths.append(os.path.join(infolder, "Yam4cr_hlx.bin"))
        elif model=="y4cs":
            output_filepaths.append(os.path.join(infolder, "Yam4csr_odd.bin"))
            output_filepaths.append(os.path.join(infolder, "Yam4csr_dbl.bin"))
            output_filepaths.append(os.path.join(infolder, "Yam4csr_vol.bin"))
            output_filepaths.append(os.path.join(infolder, "Yam4csr_hlx.bin"))
        else:
            raise(f"Invalid model!! \n model type argument must be either '' for default or Y4R or S4R")
    
    else:
        if not model or model=="y4co":
            output_filepaths.append(os.path.join(infolder, "Yam4co_odd.tif"))
            output_filepaths.append(os.path.join(infolder, "Yam4co_dbl.tif"))
            output_filepaths.append(os.path.join(infolder, "Yam4co_vol.tif"))
            output_filepaths.append(os.path.join(infolder, "Yam4co_hlx.tif"))
        elif model=="y4cr":
            output_filepaths.append(os.path.join(infolder, "Yam4cr_odd.tif"))
            output_filepaths.append(os.path.join(infolder, "Yam4cr_dbl.tif"))
            output_filepaths.append(os.path.join(infolder, "Yam4cr_vol.tif"))
            output_filepaths.append(os.path.join(infolder, "Yam4cr_hlx.tif"))
        elif model=="y4cs":
            output_filepaths.append(os.path.join(infolder, "Yam4csr_odd.tif"))
            output_filepaths.append(os.path.join(infolder, "Yam4csr_dbl.tif"))
            output_filepaths.append(os.path.join(infolder, "Yam4csr_vol.tif"))
            output_filepaths.append(os.path.join(infolder, "Yam4csr_hlx.tif"))
        else:
            raise(f"Invalid model!! \n model type argument must be either '' for default or Y4R or S4R")
            
    
    process_chunks_parallel(input_filepaths, list(output_filepaths), 
                    window_size, write_flag,
                    process_chunk_yam4cfp,
                    *[model],
                    block_size=block_size, 
                    max_workers=max_workers,  num_outputs=len(output_filepaths),
                    cog_flag=cog_flag,
                    cog_overviews=cog_overviews,
                    )


def unitary_rotation(T3, teta):
    # Direct access to each slice, no transposition or extra memory allocation
    T11 = T3[0]
    T12 = T3[1]
    T13 = T3[2]
    T22 = T3[4]
    T23 = T3[5]
    T33 = T3[8]
    
    # Compute the real and imaginary parts of T12, T13, and T23
    T12_re, T12_im = T12.real, T12.imag
    T13_re, T13_im = T13.real, T13.imag
    T23_re, T23_im = T23.real, T23.imag

    # Apply the unitary rotation directly to the 3D array (in place)
    T3[1] = T12_re * np.cos(teta) + T13_re * np.sin(teta) + 1j * (T12_im * np.cos(teta) + T13_im * np.sin(teta))
    T3[2] = -T12_re * np.sin(teta) + T13_re * np.cos(teta) + 1j * (-T12_im * np.sin(teta) + T13_im * np.cos(teta))
    T3[4] = T22 * np.cos(teta)**2 + 2 * T23_re * np.cos(teta) * np.sin(teta) + T33 * np.sin(teta)**2
    T3[5] = -T22 * np.cos(teta) * np.sin(teta) + T23_re * np.cos(teta)**2 - T23_re * np.sin(teta)**2 + T33 * np.cos(teta) * np.sin(teta) + 1j * (T23_im * np.cos(teta)**2 + T23_im * np.sin(teta)**2)
    T3[8] = T22 * np.sin(teta)**2 + T33 * np.cos(teta)**2 - 2 * T23_re * np.cos(teta) * np.sin(teta)

    return T3

# def process_chunk_yam4cfp(chunks, window_size, input_filepaths, model,*args):
def process_chunk_yam4cfp(chunks, window_size, input_filepaths,  *args, **kwargs):
    model = args[-1]
    # additional_arg1 = args[0] if len(args) > 0 else None
    # additional_arg2 = args[1] if len(args) > 1 else None

    if 'T11' in input_filepaths[0] and 'T22' in input_filepaths[5] and 'T33' in input_filepaths[8]:
        t11_T1 = np.array(chunks[0])
        t12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
        t13_T1 = np.array(chunks[3])+1j*np.array(chunks[4])
        t21_T1 = np.conj(t12_T1)
        t22_T1 = np.array(chunks[5])
        t23_T1 = np.array(chunks[6])+1j*np.array(chunks[7])
        t31_T1 = np.conj(t13_T1)
        t32_T1 = np.conj(t23_T1)
        t33_T1 = np.array(chunks[8])

        T_T1 = np.array([[t11_T1, t12_T1, t13_T1], 
                     [t21_T1, t22_T1, t23_T1], 
                     [t31_T1, t32_T1, t33_T1]])
        C3 = T3_C3_mat(T_T1)
        span = C3[0,0].real+C3[1,1].real+C3[2,2].real
        del C3


    if 'C11' in input_filepaths[0] and 'C22' in input_filepaths[5] and 'C33' in input_filepaths[8]:
        C11 = np.array(chunks[0])
        C12 = np.array(chunks[1])+1j*np.array(chunks[2])
        C13 = np.array(chunks[3])+1j*np.array(chunks[4])
        C21 = np.conj(C12)
        C22 = np.array(chunks[5])
        C23 = np.array(chunks[6])+1j*np.array(chunks[7])
        C31 = np.conj(C13)
        C32 = np.conj(C23)
        C33 = np.array(chunks[8])
        C3 = np.array([[C11, C12, C13], 
                         [C21, C22, C23], 
                         [C31, C32, C33]])
        span = C3[0,0].real+C3[1,1].real+C3[2,2].real
        T_T1 = C3_T3_mat(C3)
        del C3


    if window_size>1:
        kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

        t11f = conv2d(T_T1[0,0,:,:],kernel)
        t12f = conv2d(np.real(T_T1[0,1,:,:]),kernel)+1j*conv2d(np.imag(T_T1[0,1,:,:]),kernel)
        t13f = conv2d(np.real(T_T1[0,2,:,:]),kernel)+1j*conv2d(np.imag(T_T1[0,2,:,:]),kernel)
        
        t21f = np.conj(t12f) 
        t22f = conv2d(T_T1[1,1,:,:],kernel)
        t23f = conv2d(np.real(T_T1[1,2,:,:]),kernel)+1j*conv2d(np.imag(T_T1[1,2,:,:]),kernel)

        t31f = np.conj(t13f) 
        t32f = np.conj(t23f) 
        t33f = conv2d(T_T1[2,2,:,:],kernel)

        T_T1 = np.array([[t11f, t12f, t13f], [t21f, t22f, t23f], [t31f, t32f, t33f]])

    

    _,_,rows,cols = np.shape(T_T1)
    
    SpanMax = np.nanmax(span)
    SpanMin = np.nanmin(span)
    eps = 1e-6
    SpanMin = np.nanmax([SpanMin, eps])

    
    T_T1 = T_T1.reshape(9, rows, cols)
    T_T1 = np.dstack((T_T1[0,:,:],T_T1[1,:,:],T_T1[2,:,:],
                    T_T1[3,:,:],T_T1[4,:,:],T_T1[5,:,:],T_T1[6,:,:],T_T1[7,:,:],T_T1[8,:,:]))


    M_odd = np.zeros((rows,cols))
    M_dbl = np.zeros((rows,cols))
    M_vol =  np.zeros((rows,cols))
    M_hlx =  np.zeros((rows,cols))

    
    # print(model)
    for ii in range(rows):
        for jj in range(cols):
            T3 = np.array((
                    T_T1[ii,jj,0],T_T1[ii,jj,1],T_T1[ii,jj,2],
                    T_T1[ii,jj,3],T_T1[ii,jj,4],T_T1[ii,jj,5],
                    T_T1[ii,jj,6],T_T1[ii,jj,7],T_T1[ii,jj,8],
                        ))        
            if model in ["y4cr", "y4cs"]:
                # print(model)
                teta = 0.5 * np.arctan(2 * T3[5].real / (T3[4].real - T3[8].real))
                T3 = unitary_rotation(T3, teta)

            Pc = 2.0 * np.abs(T3[5].imag)
            HV_type = 1
            
            if model=="y4cs":
                # print(model)
                C1 = T3[0] - T3[4] + (7.0 / 8.0) * T3[8] + (Pc / 16.0)
                if C1 > 0.0:
                    HV_type = 1  # Surface scattering
                else:
                    HV_type = 2  # Double bounce scattering
            
            # Surface scattering
            if HV_type == 1:
                ratio = 10.0 * np.log10((T3[0] + T3[4] - 2.0 * T3[1].real) / (T3[0] + T3[4] + 2.0 * T3[1].real))
                
                if -2.0 < ratio <= 2.0:
                    Pv = 2.0 * (2.0 * T3[8] - Pc)
                else:
                    Pv = (15.0 / 8.0) * (2.0 * T3[8] - Pc)
            # Double bounce scattering
            if HV_type == 2:
                Pv = (15.0 / 16.0) * (2.0 * T3[8] - Pc)

            TP = np.real(T3[0] + T3[4] + T3[8])
            
            if Pv < 0.0:
                # Freeman - Yamaguchi 3-components algorithm
                HHHH = (T3[0] + 2.0 * T3[1].real + T3[4]) / 2.0
                HHVVre = (T3[0] - T3[4]) / 2.0
                HHVVim = -T3[1].imag
                HVHV = T3[8] / 2.0
                VVVV = (T3[0] - 2.0 * T3[1].real + T3[4]) / 2.0
        
                ratio = 10.0 * np.log10(VVVV / HHHH)
                
                if ratio <= -2.0:
                    FV = 15.0 * (HVHV / 4.0)
                    HHHH -= 8.0 * (FV / 15.0)
                    VVVV -= 3.0 * (FV / 15.0)
                    HHVVre -= 2.0 * (FV / 15.0)
                if ratio > 2.0:
                    FV = 15.0 * (HVHV / 4.0)
                    HHHH -= 3.0 * (FV / 15.0)
                    VVVV -= 8.0 * (FV / 15.0)
                    HHVVre -= 2.0 * (FV / 15.0)
                if -2.0 < ratio <= 2.0:
                    FV = 8.0 * (HVHV / 2.0)
                    HHHH -= 3.0 * (FV / 8.0)
                    VVVV -= 3.0 * (FV / 8.0)
                    HHVVre -= 1.0 * (FV / 8.0)
        
                # Case 1: Volume Scatter > Total
                if HHHH <= eps or VVVV <= eps:
                    FD = 0.0
                    FS = 0.0
                    if -2.0 < ratio <= 2.0:
                        FV = (HHHH + 3.0 * (FV / 8.0)) + HVHV + (VVVV + 3.0 * (FV / 8.0))
                    if ratio <= -2.0:
                        FV = (HHHH + 8.0 * (FV / 15.0)) + HVHV + (VVVV + 3.0 * (FV / 15.0))
                    if ratio > 2.0:
                        FV = (HHHH + 3.0 * (FV / 15.0)) + HVHV + (VVVV + 8.0 * (FV / 15.0))
                else:
                    # Data conditioning for non realizable ShhSvv* term
                    rtemp = HHVVre ** 2 + HHVVim ** 2
                    if rtemp > HHHH * VVVV:
                        HHVVre = HHVVre * np.sqrt((HHHH * VVVV) / rtemp)
                        HHVVim = HHVVim * np.sqrt((HHHH * VVVV) / rtemp)
        
                    # Odd Bounce
                    if HHVVre >= 0.0:
                        ALPre = -1.0
                        ALPim = 0.0
                        FD = (HHHH * VVVV - HHVVre ** 2 - HHVVim ** 2) / (HHHH + VVVV + 2.0 * HHVVre)
                        FS = VVVV - FD
                        BETre = (FD + HHVVre) / FS
                        BETim = HHVVim / FS
        
                    # Even Bounce
                    if HHVVre < 0.0:
                        BETre = 1.0
                        BETim = 0.0
                        FS = (HHHH * VVVV - HHVVre ** 2 - HHVVim ** 2) / (HHHH + VVVV - 2.0 * HHVVre)
                        FD = VVVV - FS
                        ALPre = (HHVVre - FS) / FD
                        ALPim = HHVVim / FD
        
                M_odd[ii, jj] = FS * (1 + BETre ** 2 + BETim ** 2)
                M_dbl[ii, jj] = FD * (1 + ALPre ** 2 + ALPim ** 2)
                M_vol[ii, jj] = FV
                M_hlx[ii, jj] = 0.0
        
                # # Apply Span limits
                M_odd[ii, jj] = np.clip(M_odd[ii, jj], SpanMin, SpanMax)
                M_dbl[ii, jj] = np.clip(M_dbl[ii, jj], SpanMin, SpanMax)
                M_vol[ii, jj] = np.clip(M_vol[ii, jj], SpanMin, SpanMax)
            else:
                # Surface scattering (HV_type == 1)
                if HV_type == 1:
                    S = T3[0].real - (Pv / 2.)
                    D = TP - Pv - Pc - S
                    Cre = T3[1].real + T3[2].real
                    Cim = T3[1].imag + T3[2].imag
                
                    if ratio <= -2.:
                        Cre = Cre - (Pv / 6.)
                    if ratio > 2.:
                        Cre = Cre + (Pv / 6.)
                
                    if (Pv + Pc) > TP:
                        Ps = 0.
                        Pd = 0.
                        Pv = TP - Pc
                    else:
                        CO = 2. * T3[0].real + Pc - TP
                        if CO > 0.:
                            Ps = S + (Cre * Cre + Cim * Cim) / S
                            Pd = D - (Cre * Cre + Cim * Cim) / S
                        else:
                            Pd = D + (Cre * Cre + Cim * Cim) / D
                            Ps = S - (Cre * Cre + Cim * Cim) / D
                
                    if Ps < 0.:
                        if Pd < 0.:
                            Ps = 0.
                            Pd = 0.
                            Pv = TP - Pc
                        else:
                            Ps = 0.
                            Pd = TP - Pv - Pc
                    else:
                        if Pd < 0.:
                            Pd = 0.
                            Ps = TP - Pv - Pc
                
                # Double bounce scattering (HV_type == 2)
                elif HV_type == 2:
                    S = T3[0].real
                    D = TP - Pv - Pc - S
                
                    Cre = T3[1].real + T3[2].real
                    Cim = T3[1].imag + T3[2].imag
                    
                    Pd = D + (Cre * Cre + Cim * Cim) / D
                    Ps = S - (Cre * Cre + Cim * Cim) / D
                
                    if Ps < 0.:
                        if Pd < 0.:
                            Ps = 0.
                            Pd = 0.
                            Pv = TP - Pc
                        else:
                            Ps = 0.
                            Pd = TP - Pv - Pc
                    else:
                        if Pd < 0.:
                            Pd = 0.
                            Ps = TP - Pv - Pc
                
                # Ensure values are non-negative and within bounds
                if Ps < 0.:
                    Ps = 0.
                if Pd < 0.:
                    Pd = 0.
                if Pv < 0.:
                    Pv = 0.
                if Pc < 0.:
                    Pc = 0.
                
                if Ps > SpanMax:
                    Ps = SpanMax
                if Pd > SpanMax:
                    Pd = SpanMax
                if Pv > SpanMax:
                    Pv = SpanMax
                if Pc > SpanMax:
                    Pc = SpanMax
                
                # Store the results in the output arrays
                M_odd[ii, jj] = Ps
                M_dbl[ii, jj] = Pd
                M_vol[ii, jj] = Pv
                M_hlx[ii, jj] = Pc
    # print(np.nanmean(M_odd),np.nanmean(M_dbl),np.nanmean(M_vol),np.nanmean(M_hlx))
    return M_odd, M_dbl, M_vol, M_hlx