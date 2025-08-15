import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it
from polsartools.utils.convert_matrices import T3_C3_mat
from .fp_infiles import fp_c3t3files
"""

Praks, J., Koeniguer, E.C. and Hallikainen, M.T., 2009. 
Alternatives to target entropy and alpha angle in SAR polarimetry. 
IEEE Transactions on Geoscience and Remote Sensing, 47(7), pp.2262-2274.

"""

@time_it
def praks_parm_fp(infolder,  window_size=1, outType="tif", cog_flag=False, 
          cog_overviews = [2, 4, 8, 16], write_flag=True, 
          max_workers=None,block_size=(512, 512),
          progress_callback=None,  # for QGIS plugin
          ):
    
    input_filepaths = fp_c3t3files(infolder)

    output_filepaths = []
    if outType == "bin":
        output_filepaths.append(os.path.join(infolder, "FrobeniusNorm.bin"))
        output_filepaths.append(os.path.join(infolder, "ScattPredominance.bin"))
        output_filepaths.append(os.path.join(infolder, "ScatteringDiversity.bin"))
        output_filepaths.append(os.path.join(infolder, "DegreePurity.bin"))
        output_filepaths.append(os.path.join(infolder, "DepolarizationIndex.bin"))
        output_filepaths.append(os.path.join(infolder, "Praks_Alpha.bin"))
        output_filepaths.append(os.path.join(infolder, "Praks_Entropy.bin"))
    else:
        output_filepaths.append(os.path.join(infolder, "FrobeniusNorm.tif"))
        output_filepaths.append(os.path.join(infolder, "ScattPredominance.tif"))
        output_filepaths.append(os.path.join(infolder, "ScatteringDiversity.tif"))
        output_filepaths.append(os.path.join(infolder, "DegreePurity.tif"))
        output_filepaths.append(os.path.join(infolder, "DepolarizationIndex.tif"))
        output_filepaths.append(os.path.join(infolder, "Praks_Alpha.tif"))
        output_filepaths.append(os.path.join(infolder, "Praks_Entropy.tif"))
        
    process_chunks_parallel(input_filepaths, list(output_filepaths), 
                            window_size=window_size, write_flag=write_flag,
                        processing_func=process_chunk_praks,block_size=block_size, 
                        max_workers=max_workers,  num_outputs=len(output_filepaths),
                        cog_flag=cog_flag,
                        cog_overviews=cog_overviews,
                        progress_callback=progress_callback
                        )

def process_chunk_praks(chunks, window_size, input_filepaths, *args):

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

        T3 = np.array([[t11_T1, t12_T1, t13_T1], 
                     [t21_T1, t22_T1, t23_T1], 
                     [t31_T1, t32_T1, t33_T1]])
        T_T1 = T3_C3_mat(T3)


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
        T_T1 = np.array([[C11, C12, C13], 
                         [C21, C22, C23], 
                         [C31, C32, C33]])


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

    # Span (total power)
    span = T_T1[0, 0] + T_T1[1, 1] + T_T1[2, 2]
    span_real = np.real(span) + 1e-12  # Avoid division by zero

    # Normalize M by span
    M_norm = T_T1 / span_real[None, None, :, :]

    # Frobenius norm
    FrobNorm = np.sum(np.abs(M_norm)**2, axis=(0, 1))

    # Scattering power (predominant)
    ScattPred = np.sqrt(FrobNorm)

    # Scattering diversity
    ScattDiv = 1.5 * (1. - FrobNorm)

    # Safe term for purity and depolarization index
    safe_term = np.maximum(FrobNorm - 0.25, 0)

    # Degree of polarization purity
    DegPur = 2.0 * np.sqrt(safe_term)

    # Depolarization index
    DepInd = 1. - 2.0 * np.sqrt(safe_term) / np.sqrt(3.)

    # Alpha angle (in degrees)
    Alpha = np.arccos(np.clip(np.real(M_norm[0, 0]), -1.0, 1.0)) * 180. / np.pi

    # Entropy (from determinant)
    M_adj = M_norm.copy()
    for i in range(3):
        M_adj[i, i] += 0.16

    # Move axes to compute determinant over (n, m)
    M_adj_moved = np.moveaxis(M_adj, [2, 3], [0, 1])  # shape: (n, m, 3, 3)
    det = np.linalg.det(M_adj_moved)

    Entropy = 2.52 + 0.78 * np.log(np.sqrt(np.real(det)**2 + np.imag(det)**2) + 1e-12) / np.log(3.)

    return FrobNorm.astype(np.float32),ScattPred.astype(np.float32),ScattDiv.astype(np.float32),\
        DegPur.astype(np.float32),DepInd.astype(np.float32),Alpha.astype(np.float32),Entropy.astype(np.float32)