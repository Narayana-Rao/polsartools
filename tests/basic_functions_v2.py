# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 21:54:20 2024

@author: COE
"""

from osgeo import gdal
import numpy as np
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from scipy.ndimage import uniform_filter
from threading import Lock
import time
file_lock = Lock()
def conv2d(a, f):
    filt = np.zeros(a.shape)
    wspad = int(f.shape[0]/2)
    s = f.shape + tuple(np.subtract(a.shape, f.shape) + 1)
    strd = np.lib.stride_tricks.as_strided
    subM = strd(a, shape = s, strides = a.strides * 2)
    filt_data = np.einsum('ij,ijkl->kl', f, subM)
    filt[wspad:wspad+filt_data.shape[0],wspad:wspad+filt_data.shape[1]] = filt_data
    return filt
def eig22(c2):
    c11 = c2[:,:,0].flatten()
    c12 = c2[:,:,1].flatten()
    c21 = c2[:,:,2].flatten()
    c22 = c2[:,:,3].flatten()
    trace = -(c11+c22)
    det = c11*c22-c12*c21
    # const= 1
    sqdiscr = np.sqrt(trace*trace - 4*det);
    lambda1 = -(trace + sqdiscr)*0.5;
    lambda2 = -(trace - sqdiscr)*0.5;
    
    return lambda1,lambda2

def read_chunk_with_overlap(filepath, x_start, y_start, width, height, window_size):
    """
    Read a chunk from the input raster with an overlap equal to the window size.
    Opens a new dataset for each thread.
    """
    dataset = gdal.Open(filepath, gdal.GA_ReadOnly)
    if dataset is None:
        raise FileNotFoundError(f"Cannot open {filepath}")

    band = dataset.GetRasterBand(1)

    # Ensure we don't read beyond the edges of the raster
    xoff = max(x_start - window_size, 0)
    yoff = max(y_start - window_size, 0)
    xsize = min(width + 2 * window_size, dataset.RasterXSize - xoff)
    ysize = min(height + 2 * window_size, dataset.RasterYSize - yoff)

    chunk = band.ReadAsArray(xoff=xoff, yoff=yoff, win_xsize=xsize, win_ysize=ysize)

    dataset = None  # Close the dataset after reading
    return chunk

def process_dprvi(chunks, window_size):

    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    c11_T1r = conv2d(np.real(c11_T1),kernel)
    c11_T1i = conv2d(np.imag(c11_T1),kernel)
    c11s = c11_T1r+1j*c11_T1i

    c12_T1r = conv2d(np.real(c12_T1),kernel)
    c12_T1i = conv2d(np.imag(c12_T1),kernel)
    c12s = c12_T1r+1j*c12_T1i


    c21_T1r = conv2d(np.real(c21_T1),kernel)
    c21_T1i = conv2d(np.imag(c21_T1),kernel)
    c21s = c21_T1r+1j*c21_T1i


    c22_T1r = conv2d(np.real(c22_T1),kernel)
    c22_T1i = conv2d(np.imag(c22_T1),kernel)
    c22s = c22_T1r+1j*c22_T1i

    c2_det = (c11s*c22s-c12s*c21s)
    c2_trace = c11s+c22s
    # t2_span = t11s*t22s
    m = (np.sqrt(1.0-(4.0*c2_det/np.power(c2_trace,2))))
    egv1,egv2 = eig22(np.dstack([c11s,c12s,c21s,c22s]))
    egf = np.vstack([egv1,egv2])
    egfmax = egf.max(axis=0)#.reshape(np.shape(C2_stack[:,:,0]))
    beta = (egfmax/(egv1+egv2)).reshape(np.shape(c11_T1))
    dprvi = np.real(1-(m*beta))

    return dprvi, c2_trace

def write_chunk_to_temp_file(processed_chunks, x_start, y_start, block_width, block_height, window_size, raster_width, raster_height, num_outputs):
    temp_paths = []
    for i in range(num_outputs):
        temp_fd, temp_path = tempfile.mkstemp(suffix=f'_output_{i}.tif')
        os.close(temp_fd)
        
        driver = gdal.GetDriverByName('GTiff')
        temp_dataset = driver.Create(temp_path, block_width, block_height, 1, gdal.GDT_Float32)
        geotransform = (x_start, 1, 0, y_start, 0, -1)
        temp_dataset.SetGeoTransform(geotransform)
        temp_dataset.SetProjection('')

        temp_band = temp_dataset.GetRasterBand(1)
        
        if x_start == 0 and y_start == 0:
            temp_band.WriteArray(processed_chunks[i][0:-window_size * 2, 0:-window_size * 2])
        if x_start == 0 and y_start != 0:
            temp_band.WriteArray(processed_chunks[i][window_size:-window_size, 0:-window_size * 2])
        if x_start != 0 and y_start == 0:
            temp_band.WriteArray(processed_chunks[i][0:-window_size * 2, window_size:-window_size])
        elif x_start != 0 and y_start != 0:
            temp_band.WriteArray(processed_chunks[i][window_size:-window_size, window_size:-window_size])

        temp_dataset.FlushCache()
        temp_dataset = None

        temp_paths.append(temp_path)

    return temp_paths, x_start, y_start

def merge_temp_files(output_filepaths, temp_files, raster_width, raster_height, geotransform, projection, num_outputs):
    for i in range(num_outputs):
        driver = gdal.GetDriverByName('GTiff')
        output_dataset = driver.Create(output_filepaths[i], raster_width, raster_height, 1, gdal.GDT_Float32)
        output_dataset.SetGeoTransform(geotransform)
        output_dataset.SetProjection(projection)
        output_band = output_dataset.GetRasterBand(1)

        mask = np.zeros((raster_height, raster_width), dtype=np.float32)

        for temp_path_set, x_start, y_start in temp_files:
            temp_dataset = gdal.Open(temp_path_set[i], gdal.GA_ReadOnly)
            temp_band = temp_dataset.GetRasterBand(1)

            temp_chunk = temp_band.ReadAsArray()
            temp_height, temp_width = temp_chunk.shape

            output_band.WriteArray(temp_chunk[:temp_height, :temp_width], xoff=x_start, yoff=y_start)
            mask[y_start:y_start + temp_height, x_start:x_start + temp_width] += temp_chunk[:temp_height, :temp_width]

            temp_dataset = None

        mask[mask == 0] = np.nan
        mask[np.isnan(mask)] = 0

        output_dataset.FlushCache()
        output_dataset = None

def process_and_write_chunk(args, processing_func, num_outputs):
    (input_filepaths, x_start, y_start, read_block_width, read_block_height, window_size, raster_width, raster_height) = args

    chunks = [read_chunk_with_overlap(fp, x_start, y_start, read_block_width, read_block_height, window_size) for fp in input_filepaths]

    processed_chunks = processing_func(chunks, window_size)

    temp_paths, temp_x_start, temp_y_start = write_chunk_to_temp_file(
        processed_chunks,
        x_start, y_start, read_block_width, read_block_height, window_size, raster_width, raster_height, num_outputs
    )

    return temp_paths, temp_x_start, temp_y_start

def process_chunks_parallel(input_filepaths,   output_filepaths, window_size=3,write_flag=True,block_size=(512, 512), max_workers=None, processing_func=process_dprvi, num_outputs=1):
    if len(input_filepaths) not in [2, 4, 9]:
        raise ValueError("This function only supports 2, 4, or 9 input rasters.")

    if max_workers is None:
        max_workers = os.cpu_count()

    input_datasets = [gdal.Open(fp, gdal.GA_ReadOnly) for fp in input_filepaths]
    if any(ds is None for ds in input_datasets):
        raise FileNotFoundError("One or more input files could not be opened.")

    raster_width = input_datasets[0].RasterXSize
    raster_height = input_datasets[0].RasterYSize
    geotransform = input_datasets[0].GetGeoTransform()
    projection = input_datasets[0].GetProjection()

    # num_outputs = 2  # Adjust based on the number of outputs from the processing function

    tasks = []
    merged_arrays = [np.zeros((raster_height, raster_width), dtype=np.float32) for _ in range(num_outputs)] if not write_flag else None
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for y in range(0, raster_height, block_size[1]):
            for x in range(0, raster_width, block_size[0]):
                read_block_width = min(block_size[0], raster_width - x)
                read_block_height = min(block_size[1], raster_height - y)

                args = (input_filepaths, x, y, read_block_width, read_block_height, window_size, raster_width, raster_height)
                tasks.append(executor.submit(process_and_write_chunk, args, processing_func, num_outputs))

        temp_files = []
        for future in as_completed(tasks):
            temp_paths, x_start, y_start = future.result()
            if write_flag:
                temp_files.append((temp_paths, x_start, y_start))
            else:
                for i in range(num_outputs):
                    temp_dataset = gdal.Open(temp_paths[i], gdal.GA_ReadOnly)
                    temp_band = temp_dataset.GetRasterBand(1)
                    temp_chunk = temp_band.ReadAsArray()

                    temp_height, temp_width = temp_chunk.shape
                    merged_arrays[i][y_start:y_start + temp_height, x_start:x_start + temp_width] = temp_chunk
                    temp_dataset = None
                    os.remove(temp_paths[i])

    if write_flag:
        merge_temp_files(output_filepaths, temp_files, raster_width, raster_height, geotransform, projection, num_outputs)
        for temp_path_set, _, _ in temp_files:
            for temp_path in temp_path_set:
                os.remove(temp_path)
    else:
        return merged_arrays




#%%%%%%%%%%%%%%%%%%%%%%
# def dprvi(C2_folder,window_size=1,write_flag=None):

t0 = time.time()

inFolder = r"C:\Users\nbhogapurapu\Desktop\polsartools\tests\winnip_31606_12061_006_120717_L090_CX_07_grd\C2"


# Example usage for 4 input rasters, with a 3x3 window and autodetected cores
input_filepaths = [os.path.join(inFolder,"C11.bin"), 
                    os.path.join(inFolder,"C12_real.bin"), os.path.join(inFolder,"C12_imag.bin"),
                   os.path.join(inFolder,"C22.bin"),
                   # os.path.join(inFolder,"C11.bin"),
                   # os.path.join(inFolder,"C11.bin")
                   ]
# Example usage
# input_filepaths = ["input_raster_1.tif", "input_raster_2.tif"]
# output_filepath1 = "dprvi.tif"
# process_chunks_parallel(input_filepaths, output_filepath1,  block_size=(512, 512), window_size=1, processing_func=process_dprvi)


window_size = 1

output_filepaths = ["output_dprvi.tif", "output_second.tif"]
# process_chunks_parallel(input_filepaths, output_filepaths, block_size=(512, 512), window_size=1, processing_func=process_dprvi, write_flag=True, )
process_chunks_parallel(input_filepaths, output_filepaths, window_size=3,write_flag=True,block_size=(512, 512), max_workers=None, processing_func=process_dprvi, num_outputs=2
                             )
    
print("%0.2f"%(time.time()-t0))
















