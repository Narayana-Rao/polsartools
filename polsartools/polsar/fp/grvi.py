import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel,conv2d,eig22

def grvi(infolder, outname=None, window_size=1,max_workers=None):
    
    input_filepaths = [
    os.path.join(infolder, "T11.bin"),
    os.path.join(infolder, 'T12_imag.bin'),
    os.path.join(infolder, 'T12_real.bin'),
    os.path.join(infolder, 'T13_imag.bin'),
    os.path.join(infolder, 'T13_real.bin'),
    
    os.path.join(infolder, "T22.bin"),
    os.path.join(infolder, 'T23_imag.bin'),
    os.path.join(infolder, 'T23_real.bin'),
    
    os.path.join(infolder, "T33.bin"),

    ]
    if outname is None:
        outname = os.path.join(infolder, "grvi.tif")
    process_chunks_parallel(outname, input_filepaths, processing_func=process_chunk_grvi, window_size=window_size, max_workers=max_workers)

def process_chunk_grvi(chunks, window_size):
    T11 = np.array(chunks[0])
    T12_i = np.array(chunks[1])
    T12_r = np.array(chunks[2])
    T13_i = np.array(chunks[3])
    T13_r = np.array(chunks[4])

    T22 = np.array(chunks[5])
    T23_i = np.array(chunks[6])
    T23_r = np.array(chunks[7])

    T33 = np.array(chunks[8])

    T12 = T12_r + 1j * T12_i
    T13 = T13_r + 1j * T13_i
    T23 = T23_r + 1j * T23_i

    # Create complex coherency matrix
    t11_T1 = T11
    t12_T1 = T12
    t13_T1 = T13
    t21_T1 = np.conj(T12)
    t22_T1 = T22
    t23_T1 = T23
    t31_T1 = np.conj(T13)
    t32_T1 = np.conj(T23)
    t33_T1 = T33

    nrows, ncols = T11.shape
    # Initialize arrays
    span = np.zeros((nrows, ncols))
    GD_t1_rv = np.zeros((nrows, ncols))
    GD_t1_c = np.zeros((nrows, ncols))
    GD_t1_t = np.zeros((nrows, ncols))
    GD_t1_d = np.zeros((nrows, ncols))
    GD_t1_nd = np.zeros((nrows, ncols))
    beta = np.zeros((nrows, ncols))
    a = np.zeros((nrows, ncols))
    b = np.zeros((nrows, ncols))
    
    # Define matrices
    D = (1 / np.sqrt(2)) * np.array([[1, 0, 1], [1, 0, -1], [0, np.sqrt(2), 0]])

    M_d = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]])
    M_nd = np.array([[0.625, 0.375, 0, 0], [0.375, 0.625, 0, 0], [0, 0, -0.5, 0], [0, 0, 0, 0.5]])
    M_t = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
    M_c = np.array([[0.625, 0.375, 0, 0], [0.375, 0.625, 0, 0], [0, 0, 0.5, 0], [0, 0, 0, -0.5]])
    M_lh = np.array([[1, 0, 0, -1], [0, 0, 0, 0], [0, 0, 0, 0], [-1, 0, 0, 1]])
    M_rh = np.array([[1, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 1]])

    # # Vectorize window mean calculations
    # def apply_window_mean(array, size):
    #     padded_array = np.pad(array, pad_width=size, mode='constant', constant_values=np.nan)
    #     shape = (padded_array.shape[0] - size * 2, padded_array.shape[1] - size * 2)
    #     strides = padded_array.strides * 2
    #     windows = np.lib.stride_tricks.as_strided(padded_array, shape=shape + (size * 2, size * 2), strides=strides)
    #     return np.nanmean(windows, axis=(-1, -2))

    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    t11s = conv2d(t11_T1, kernel)
    t12s = conv2d(t12_T1, kernel)
    t13s = conv2d(t13_T1, kernel)
    t21s = conv2d(t21_T1, kernel)
    t22s = conv2d(t22_T1, kernel)
    t23s = conv2d(t23_T1, kernel)
    t31s = conv2d(t31_T1, kernel)
    t32s = conv2d(t32_T1, kernel)
    t33s = conv2d(t33_T1, kernel)

     # Coherency matrix
    D = (1 / np.sqrt(2)) * np.array([[1, 0, 1], [1, 0, -1], [0, np.sqrt(2), 0]])
    T_T1 = np.array([[t11s, t12s, t13s], [t21s, t22s, t23s], [t31s, t32s, t33s]])
    C_T1 = np.matmul(np.matmul(D.T, T_T1), D)

    span = np.real(t11s + t22s + t33s)

    # Kennaugh Matrix
    M_T = np.array([[t11s + t22s + t33s, t12s + t21s, t13s + t31s, -1j * (t23s - t32s)],
                    [t12s + t21s, t11s + t22s - t33s, t23s + t32s, -1j * (t13s - t31s)],
                    [t13s + t31s, t23s + t32s, t11s - t22s + t33s, 1j * (t12s - t21s)],
                    [-1j * (t23s - t32s), -1j * (t13s - t31s), 1j * (t12s - t21s), -t11s + t22s + t33s]])

    # Gamma/Rho
    def compute_gd(M_T_theta, M):
        num1 = np.matmul(M_T_theta.T, M)
        num = np.trace(num1, axis1=1, axis2=2)
        den1 = np.sqrt(np.abs(np.trace(np.matmul(M_T_theta.T, M_T_theta), axis1=1, axis2=2)))
        den2 = np.sqrt(np.abs(np.trace(np.matmul(M.T, M), axis1=1, axis2=2)))
        return np.real(2 * np.arccos(num / (den1 * den2)) * 180 / np.pi) / 180

    GD_t1_rv = compute_gd(M_T, M_rv)
    GD_t1_c = compute_gd(M_T, np.array([[0.625, 0.375, 0, 0], [0.375, 0.625, 0, 0], [0, 0, 0.5, 0], [0, 0, 0, -0.5]]))
    GD_t1_t = compute_gd(M_T, np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]))
    GD_t1_d = compute_gd(M_T, np.array([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]))
    GD_t1_nd = compute_gd(M_T, np.array([[0.625, 0.375, 0, 0], [0.375, 0.625, 0, 0], [0, 0, -0.5, 0], [0, 0, 0, 0.5]]))

    # VI calculation
    beta = np.nan_to_num((np.nanmin([GD_t1_t, GD_t1_d, GD_t1_c, GD_t1_nd], axis=0) /
                          np.nanmax([GD_t1_t, GD_t1_d, GD_t1_c, GD_t1_nd], axis=0))**2)
    vi = np.power(beta, GD_t1_rv) * (1 - (1 / np.nan_to_num(f)) * GD_t1_rv)

    idx1 = GD_t1_rv > np.nan_to_num(f)
    vi[idx1] = 0
    vi[~idx1] = vi[~idx1]

    return np.real(vi)