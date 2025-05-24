import os, tempfile
from osgeo import gdal
import numpy as np
from functools import wraps
import time

def time_it(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        try:
            result = func(*args, **kwargs)
            end_time = time.time()
            duration = int(end_time - start_time)

            if duration < 60:
                print(f"Execution time for {func.__name__}: {duration:.2f} seconds")
            elif duration < 3600:
                minutes, seconds = divmod(duration, 60)
                print(f"Execution time for {func.__name__}: {minutes:02}:{seconds:02} (MM:SS)")
            else:
                hours, remainder = divmod(duration, 3600)
                minutes, seconds = divmod(remainder, 60)
                print(f"Execution time for {func.__name__}: {hours:02}:{minutes:02}:{seconds:02} (HH:MM:SS)")
                
            return result
        except Exception as e:
            raise e
    return wrapper

def conv2d(a, f):
    filt = np.zeros(a.shape)
    wspad = int(f.shape[0]/2)
    s = f.shape + tuple(np.subtract(a.shape, f.shape) + 1)
    strd = np.lib.stride_tricks.as_strided
    subM = strd(a, shape=s, strides=a.strides * 2)
    filt_data = np.einsum('ij,ijkl->kl', f, subM)
    filt[wspad:wspad+filt_data.shape[0], wspad:wspad+filt_data.shape[1]] = filt_data
    return filt

def eig22(c2):
    c11 = c2[:, :, 0].flatten()
    c12 = c2[:, :, 1].flatten()
    c21 = c2[:, :, 2].flatten()
    c22 = c2[:, :, 3].flatten()
    trace = -(c11 + c22)
    det = c11 * c22 - c12 * c21
    sqdiscr = np.sqrt(trace * trace - 4 * det)
    lambda1 = -(trace + sqdiscr) * 0.5
    lambda2 = -(trace - sqdiscr) * 0.5
    return lambda1, lambda2
 
