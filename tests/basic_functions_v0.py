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

def process_chunk(chunks, window_size):
    combined_chunk = np.mean(np.array(chunks), axis=0)
    filtered_chunk = uniform_filter(combined_chunk, size=window_size, mode='constant')
    return filtered_chunk

def write_chunk_to_temp_file(processed_chunk, x_start, y_start, block_width, block_height, window_size, raster_width, raster_height):
    temp_fd, temp_path = tempfile.mkstemp(suffix='.tif')
    os.close(temp_fd)

    driver = gdal.GetDriverByName('GTiff')
    temp_dataset = driver.Create(temp_path, block_width, block_height, 1, gdal.GDT_Float32)
    geotransform = (x_start, 1, 0, y_start, 0, -1)
    temp_dataset.SetGeoTransform(geotransform)
    temp_dataset.SetProjection('')

    # Write only the non-overlapping part of the chunk
    temp_band = temp_dataset.GetRasterBand(1)
    if x_start==0 and y_start==0:
        temp_band.WriteArray(processed_chunk[0:-window_size*2, 0:-window_size*2])
    if x_start==0 and y_start!=0:
        temp_band.WriteArray(processed_chunk[window_size:-window_size, 0:-window_size*2])
    if x_start!=0 and y_start==0:
        temp_band.WriteArray(processed_chunk[0:-window_size*2, window_size:-window_size])
    elif x_start!=0 and y_start!=0:
        temp_band.WriteArray(processed_chunk[window_size:-window_size, window_size:-window_size])
    # print(processed_chunk.shape,x_start,y_start)
    temp_dataset.FlushCache()
    temp_dataset = None

    return temp_path, x_start, y_start


def merge_temp_files(output_filepath, temp_files, raster_width, raster_height, geotransform, projection):
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_filepath, raster_width, raster_height, 1, gdal.GDT_Float32)
    output_dataset.SetGeoTransform(geotransform)
    output_dataset.SetProjection(projection)
    output_band = output_dataset.GetRasterBand(1)

    # Create a mask to handle overlaps
    mask = np.zeros((raster_height, raster_width), dtype=np.float32)

    for temp_path, x_start, y_start in temp_files:
        temp_dataset = gdal.Open(temp_path, gdal.GA_ReadOnly)
        temp_band = temp_dataset.GetRasterBand(1)

        temp_chunk = temp_band.ReadAsArray()
        temp_height, temp_width = temp_chunk.shape

        if x_start < 0:
            print(x_start)
            temp_width += x_start
            # print(x_start,temp_width)
            x_start = 0
        if y_start < 0:
            temp_height += y_start
            # print(y_start,temp_height)
            y_start = 0

        if x_start + temp_width > raster_width:
            temp_width = raster_width - x_start
        if y_start + temp_height > raster_height:
            temp_height = raster_height - y_start

        # Read and write only the valid part
        output_band.WriteArray(temp_chunk[:temp_height, :temp_width], xoff=x_start, yoff=y_start)
        mask[y_start:y_start + temp_height, x_start:x_start + temp_width] += temp_chunk[:temp_height, :temp_width]

        temp_dataset = None

    # Normalize mask to avoid writing artifacts
    mask[mask == 0] = np.nan
    mask[np.isnan(mask)] = 0

    output_dataset.FlushCache()
    output_dataset = None

def process_and_write_chunk(args):
    (input_filepaths, x_start, y_start, read_block_width, read_block_height, window_size, raster_width, raster_height) = args

    # Read chunks with overlap
    chunks = [read_chunk_with_overlap(fp, x_start, y_start, read_block_width, read_block_height, window_size) for fp in input_filepaths]

    # Process the chunks
    processed_chunk = process_chunk(chunks, window_size)

    # Write the processed chunk to a temporary file
    temp_path, temp_x_start, temp_y_start = write_chunk_to_temp_file(
        processed_chunk,
        x_start, y_start, read_block_width, read_block_height, window_size, raster_width, raster_height
    )

    return temp_path, temp_x_start, temp_y_start

def process_rasters_in_chunks_parallel(output_filepath, input_filepaths, block_size=(512, 512), window_size=3, max_workers=None):
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

    tasks = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for y in range(0, raster_height, block_size[1]):
            for x in range(0, raster_width, block_size[0]):
                read_block_width = min(block_size[0], raster_width - x)
                read_block_height = min(block_size[1], raster_height - y)

                args = (input_filepaths, x, y, read_block_width, read_block_height, window_size, raster_width, raster_height)
                tasks.append(executor.submit(process_and_write_chunk, args))

        temp_files = []
        for future in as_completed(tasks):
            temp_path, x_start, y_start = future.result()
            temp_files.append((temp_path, x_start, y_start))
            # print(temp_path)

    merge_temp_files(output_filepath, temp_files, raster_width, raster_height, geotransform, projection)

    for temp_path, _, _ in temp_files:
        os.remove(temp_path)

def process_whole_raster(input_filepaths, output_filepath, window_size):
    if len(input_filepaths) not in [2, 4, 9]:
        raise ValueError("This function only supports 2, 4, or 9 input rasters.")
    
    # Open all input rasters
    datasets = [gdal.Open(fp, gdal.GA_ReadOnly) for fp in input_filepaths]
    if any(ds is None for ds in datasets):
        raise FileNotFoundError("One or more input files could not be opened.")

    # Read the raster dimensions and geotransform from the first dataset
    raster_width = datasets[0].RasterXSize
    raster_height = datasets[0].RasterYSize
    geotransform = datasets[0].GetGeoTransform()
    projection = datasets[0].GetProjection()

    # Initialize arrays to hold the data from all rasters
    combined_data = np.zeros((raster_height, raster_width, len(input_filepaths)), dtype=np.float32)

    # Read data from each raster
    for i, dataset in enumerate(datasets):
        band = dataset.GetRasterBand(1)
        combined_data[:, :, i] = band.ReadAsArray()

    # Compute the average across rasters
    combined_data_mean = np.mean(combined_data, axis=2)

    # Apply uniform filter
    filtered_raster = uniform_filter(combined_data_mean, size=window_size, mode='constant')

    # Write output raster
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_filepath, raster_width, raster_height, 1, gdal.GDT_Float32)
    output_dataset.SetGeoTransform(geotransform)
    output_dataset.SetProjection(projection)

    output_band = output_dataset.GetRasterBand(1)
    output_band.WriteArray(filtered_raster)

    output_dataset.FlushCache()
    output_dataset = None

    # Close all datasets
    for dataset in datasets:
        dataset = None

def compute_image_difference(img1_path, img2_path):
    dataset1 = gdal.Open(img1_path, gdal.GA_ReadOnly)
    dataset2 = gdal.Open(img2_path, gdal.GA_ReadOnly)
    
    band1 = dataset1.GetRasterBand(1)
    band2 = dataset2.GetRasterBand(1)
    
    array1 = band1.ReadAsArray()
    array2 = band2.ReadAsArray()
    
    difference = np.abs(array1 - array2)
    max_difference = np.max(difference)
    
    print(f"Maximum difference: {max_difference}")
#%%%%%%%%%%%%%%%%%%%%%%

# Example usage to compare outputs:
# compute_image_difference("output_raster_parallel_v9_512.tif", "whole_raster_output.tif")

t0 = time.time()
inFolder = r"C:\Users\nbhogapurapu\Documents\GitHub\PolSAR-tools\sample_data\dual_pol\C2_VVVH"
inFolder = r"C:\Users\nbhogapurapu\Desktop\polsartools\tests\winnip_31606_12061_006_120717_L090_CX_07_grd\C2"


# Example usage for 4 input rasters, with a 3x3 window and autodetected cores
input_filepaths = [os.path.join(inFolder,"C11.bin"), 
                   # os.path.join(inFolder,"C12_real.bin"), os.path.join(inFolder,"C12_imag.bin"),
                   os.path.join(inFolder,"C11.bin"),
                   # os.path.join(inFolder,"C11.bin"),
                   # os.path.join(inFolder,"C11.bin")
                   ]
# Example usage
# input_filepaths = ["input_raster_1.tif", "input_raster_2.tif"]
output_filepath1 = "output_raster_parallel_v92_512.tif"
process_rasters_in_chunks_parallel(output_filepath1, input_filepaths, block_size=(512, 512), window_size=1)
output_filepath2 = "output_raster_parallel_v92_512_val.tif"
process_whole_raster(input_filepaths, output_filepath2, window_size=1)


compute_image_difference(output_filepath1, output_filepath2)

print("%0.2f"%(time.time()-t0))