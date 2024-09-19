import os, tempfile
from osgeo import gdal
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

def process_chunks_parallel(output_filepath, input_filepaths, processing_func, block_size=(512, 512), window_size=1, max_workers=None):
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
                tasks.append(executor.submit(process_and_write_chunk, args, processing_func))

        temp_files = []
        for future in as_completed(tasks):
            temp_path, x_start, y_start = future.result()
            temp_files.append((temp_path, x_start, y_start))

    merge_temp_files(output_filepath, temp_files, raster_width, raster_height, geotransform, projection)

    for temp_path, _, _ in temp_files:
        os.remove(temp_path)

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

def process_chunks_with_function(chunks, window_size, processing_func):
    """
    Process chunks using a specified processing function.
    
    Args:
        chunks (list of numpy.ndarray): List of chunks to process.
        window_size (int): Size of the window for processing.
        processing_func (callable): The processing function to apply.
        
    Returns:
        numpy.ndarray: Processed chunk.
    """
    return processing_func(chunks, window_size)

def process_and_write_chunk(args,processing_func):
    (input_filepaths, x_start, y_start, read_block_width, read_block_height, window_size, raster_width, raster_height) = args

    # Read chunks with overlap
    chunks = [read_chunk_with_overlap(fp, x_start, y_start, read_block_width, read_block_height, window_size) for fp in input_filepaths]

    # Process the chunks using the provided processing function
    processed_chunk = process_chunks_with_function(chunks, window_size, processing_func)

    # Write the processed chunk to a temporary file
    temp_path, temp_x_start, temp_y_start = write_chunk_to_temp_file(
        processed_chunk,
        x_start, y_start, read_block_width, read_block_height, window_size, raster_width, raster_height
    )

    return temp_path, temp_x_start, temp_y_start

