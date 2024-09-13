import os
import numpy as np
import xarray as xr
import rasterio

def read_raster_in_chunks(filepath, chunk_size):
    """
    Read a raster file in chunks using xarray.
    
    Parameters:
        filepath (str): Path to the raster file.
        chunk_size (tuple): Size of the chunks to read (rows, columns).
        
    Returns:
        xarray.DataArray: Data array containing the raster data.
    """
    with xr.open_dataset(filepath, engine='rasterio') as ds:
        # Read the first band in chunks
        data_array = ds.rio.read(1, chunks=chunk_size)
    
    return data_array
    
def process_raster_in_chunks(input_file, output_file, param1, param2, chunk_size):
    """
    Process a raster file in chunks and save the result to a new file using rasterio.
    
    Parameters:
        input_file (str): Path to the input raster file.
        output_file (str): Path to the output raster file.
        param1 (int): Example parameter for processing.
        param2 (int): Example parameter for processing.
        chunk_size (tuple): Size of the chunks to process (rows, columns).
    """
    with rasterio.open(input_file) as src:
        meta = src.meta.copy()
        dtype = str(src.dtypes[0])  # Get data type of the raster
        
        # Read in chunks and process
        for i in range(0, src.height, chunk_size[0]):
            for j in range(0, src.width, chunk_size[1]):
                window = rasterio.windows.Window(j, i, chunk_size[1], chunk_size[0])
                chunk = src.read(1, window=window)
                
                # Example processing: invert values
                processed_chunk = np.flipud(chunk)
                
                # Write processed chunk to output
                with rasterio.open(output_file, 'w', **meta) as dst:
                    dst.write(processed_chunk, 1, window=window)

# def process_raster_in_chunks(input_file, output_file, param1, param2, chunk_size):
#     """
#     Process a raster file in chunks and save the result to a new file.
    
#     Parameters:
#         input_file (str): Path to the input raster file.
#         output_file (str): Path to the output raster file.
#         param1 (int): Example parameter for processing.
#         param2 (int): Example parameter for processing.
#         chunk_size (tuple): Size of the chunks to process (rows, columns).
#     """
#     data_array = read_raster_in_chunks(input_file, chunk_size)
    
#     # Get data type from xarray DataArray
#     dtype = str(data_array.dtype)
    
#     # Example processing: invert values and apply some transformation
#     processed_data = np.flipud(data_array)

#     # Save the processed raster
#     with rasterio.open(input_file) as src:
#         meta = src.meta.copy()
#         meta.update({'count': 1, 'dtype': dtype})

#         with rasterio.open(output_file, 'w', **meta) as dst:
#             dst.write(processed_data, 1)
