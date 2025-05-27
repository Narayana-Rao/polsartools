import os, tempfile
from osgeo import gdal
import numpy as np
from concurrent.futures import ProcessPoolExecutor,as_completed
from tqdm import tqdm

def process_chunks_parallel(
                            input_filepaths,
                            output_filepaths,
                            window_size,
                            write_flag,
                            processing_func,
                            *proc_args,  # NEW
                            block_size=(512, 512),
                            max_workers=None,
                            num_outputs=1,
                            cog_flag=False,
                            cog_overviews=[2, 4, 8, 16],
                            post_proc=None,  # NEW
                            out_x_size=None, out_y_size=None,
                            out_geotransform=None,
                            out_projection=None,
                            **proc_kwargs  # NEW
                        ):


    if len(input_filepaths) not in [2, 4, 9]:
        raise ValueError("This function only supports 2, 4, or 9 input rasters.")

    if max_workers is None:
        max_workers = os.cpu_count()-1  # Use all available CPUs
        # max_workers = 1
    input_datasets = [gdal.Open(fp, gdal.GA_ReadOnly) for fp in input_filepaths]
    if any(ds is None for ds in input_datasets):
        raise FileNotFoundError("One or more input files could not be opened.")

    raster_width = input_datasets[0].RasterXSize
    raster_height = input_datasets[0].RasterYSize
    geotransform = input_datasets[0].GetGeoTransform()
    projection = input_datasets[0].GetProjection()



    # # Override if resized output is requested
    output_width = out_x_size if out_x_size else raster_width
    output_height = out_y_size if out_y_size else raster_height
    output_geotransform = out_geotransform if out_geotransform else geotransform
    output_projection = out_projection if out_projection else projection

    # Ensure block_size does not exceed raster dimensions
    adjusted_block_size_x = min(block_size[0], raster_width)
    adjusted_block_size_y = min(block_size[1], raster_height)

    merged_arrays = [np.zeros((output_height, output_width), dtype=np.float32) for _ in range(num_outputs)] if not write_flag else None

    tasks = []

    """ with progress bar"""
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
    # with ThreadPoolExecutor(max_workers=max_workers) as executor:
        tasks = []
        for y in range(0, raster_height, adjusted_block_size_y):
            for x in range(0, raster_width, adjusted_block_size_x):
                read_block_width = min(adjusted_block_size_x, raster_width - x)
                read_block_height = min(adjusted_block_size_y, raster_height - y)
                # args_ = (input_filepaths, x, y, read_block_width, read_block_height, window_size, raster_width, raster_height, chi_in, psi_in,model)
                args_ = ( input_filepaths, x, y, read_block_width, read_block_height, window_size,  raster_width, raster_height, *proc_args )
                
                # print(f"Submitting task for chunk at ({x}, {y})")
                tasks.append(executor.submit(process_and_write_chunk, args_, processing_func, num_outputs, output_width,  output_height,**proc_kwargs))
        # Initialize tqdm progress bar with the total number of tasks
        with tqdm(total=len(tasks), desc=f"Progress ", unit=" block") as pbar:
            temp_files = []
            for future in as_completed(tasks):
                # try:
                    result = future.result()
                    if result is None:
                        raise ValueError("Block processing returned None")
                    temp_paths, x_start, y_start = result
                    if write_flag:
                        temp_files.append((temp_paths, x_start, y_start))
                    else:
                        for i in range(num_outputs):
                            temp_dataset = gdal.Open(temp_paths[i], gdal.GA_ReadOnly)
                            temp_band = temp_dataset.GetRasterBand(1)
                            temp_chunk = temp_band.ReadAsArray()
                            print(type(temp_chunk))
                            temp_height, temp_width = temp_chunk.shape
                            merged_arrays[i][y_start:y_start + temp_height, x_start:x_start + temp_width] = temp_chunk
                            temp_dataset = None
                            os.remove(temp_paths[i])
                    
                    pbar.update(1)
                # except Exception as e:
                #     print(f"Error in processing task: {e}")
            
    azlks = proc_kwargs.get('azlks', 1)
    rglks = proc_kwargs.get('rglks', 1)
    if write_flag:
        merge_temp_files(output_filepaths, temp_files,  output_width, output_height,     output_geotransform, output_projection, num_outputs,cog_flag,cog_overviews,azlks=azlks, rglks=rglks )
        for temp_path_set, _, _ in temp_files:
            for temp_path in temp_path_set:
                os.remove(temp_path)
                
        if post_proc:
            post_proc(output_filepaths)
    else:
        return merged_arrays

def read_chunk_with_overlap(filepath, x_start, y_start, width, height, window_size):
    dataset = gdal.Open(filepath, gdal.GA_ReadOnly)
    if dataset is None:
        raise FileNotFoundError(f"Cannot open {filepath}")
    
    band = dataset.GetRasterBand(1)
    
    
    if not window_size:  # None or 0
        chunk = band.ReadAsArray(x_start, y_start, width, height)
        dataset = None
        return chunk


    
    # Handle edge cases for window overlap near the borders of the raster
    xoff = max(x_start - window_size, 0)
    yoff = max(y_start - window_size, 0)
    
    xsize = min(width + 2 * window_size, dataset.RasterXSize - xoff)
    ysize = min(height + 2 * window_size, dataset.RasterYSize - yoff)


    # Check if the raster size is less than or equal to  block size read entire raster
    if width >= dataset.RasterXSize and height >= dataset.RasterYSize:
        chunk = band.ReadAsArray(xoff=xoff, yoff=yoff, win_xsize=xsize, win_ysize=ysize)
        dataset = None  
        return chunk


    # Determine if the chunk is near the borders of the dataset
    is_left_border = x_start == 0
    is_right_border = x_start + width >= dataset.RasterXSize
    is_top_border = y_start == 0
    is_bottom_border = y_start + height >= dataset.RasterYSize

    # Handle corner cases
    if is_left_border and is_top_border:
        # Top-left corner: no overlap on top and left
        xoff = x_start
        yoff = y_start
        xsize = min(width + window_size, dataset.RasterXSize - xoff)
        ysize = min(height + window_size, dataset.RasterYSize - yoff)
    elif is_left_border and is_bottom_border:
        # Bottom-left corner: no overlap on bottom and left
        xoff = x_start
        yoff = max(y_start - window_size, 0)
        xsize = min(width + window_size, dataset.RasterXSize - xoff)
        ysize = min(height + window_size, dataset.RasterYSize - yoff)
    elif is_right_border and is_top_border:
        # Top-right corner: no overlap on top and right
        xoff = max(x_start - window_size, 0)
        yoff = y_start
        xsize = min(width + window_size, dataset.RasterXSize - xoff)
        ysize = min(height + window_size, dataset.RasterYSize - yoff)
    elif is_right_border and is_bottom_border:
        # Bottom-right corner: no overlap on bottom and right
        xoff = max(x_start - window_size, 0)
        yoff = max(y_start - window_size, 0)
        xsize = min(width + window_size, dataset.RasterXSize - xoff)
        ysize = min(height + window_size, dataset.RasterYSize - yoff)
    elif is_left_border:
        # Left edge: no overlap on left, adjust only horizontally
        xoff = x_start
        xsize = min(width + window_size, dataset.RasterXSize - xoff)
    elif is_right_border:
        # Right edge: no overlap on right, adjust only horizontally
        xsize = min(width + window_size, dataset.RasterXSize - x_start)
    elif is_top_border:
        # Top edge: no overlap on top, adjust only vertically
        yoff = y_start
        ysize = min(height + window_size, dataset.RasterYSize - yoff)
    elif is_bottom_border:
        # Bottom edge: no overlap on bottom, adjust only vertically
        ysize = min(height + window_size, dataset.RasterYSize - y_start)

    # Read the adjusted chunk with overlap
    chunk = band.ReadAsArray(xoff=xoff, yoff=yoff, win_xsize=xsize, win_ysize=ysize)

    dataset = None  # Close the dataset after reading
    return chunk

def write_chunk_to_temp_file(processed_chunks, x_start, y_start, block_width, block_height, window_size, output_width, output_height, num_outputs):
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


        if not window_size:
            temp_band.WriteArray(processed_chunks[i])
            temp_dataset.FlushCache()
            temp_dataset = None
            temp_paths.append(temp_path)
            continue
        # Check if the block size is equal to the raster size
        # if block_width >= output_width and block_height >= output_height:     
        #     temp_band.WriteArray(processed_chunks[i]) 
        #     temp_dataset.FlushCache()
        #     temp_dataset = None
        #     temp_paths.append(temp_path)
        #     continue  # Skip to the next output
        
        # Determine the boundaries
        
        if window_size:
            is_left_border = x_start == 0
            is_right_border = x_start + block_width >= output_width
            is_top_border = y_start == 0
            is_bottom_border = y_start + block_height >= output_height

            if is_left_border and is_top_border:
                # Top-left corner: reduce only from right and bottom
                temp_band.WriteArray(processed_chunks[i][:-window_size, :][:, :-window_size])
            elif is_left_border and is_bottom_border:
                # Bottom-left corner: reduce only from right and top
                temp_band.WriteArray(processed_chunks[i][window_size:, :][:, :-window_size])
            elif is_right_border and is_top_border:
                # Top-right corner: reduce only from left and bottom
                temp_band.WriteArray(processed_chunks[i][:-window_size, :][:, window_size:])
            elif is_right_border and is_bottom_border:
                # Bottom-right corner: reduce only from left and top
                temp_band.WriteArray(processed_chunks[i][window_size:, :][:, window_size:])
            elif is_left_border:
                # Left edge: reduce from right, top, and bottom
                temp_band.WriteArray(processed_chunks[i][window_size:-window_size,:-window_size ])  # Skip left
            elif is_right_border:
                # Right edge: reduce from left, top, and bottom
                temp_band.WriteArray(processed_chunks[i][ window_size:-window_size,window_size:])  # Skip right
            elif is_top_border:
                # Top edge: reduce from left, right, and bottom
                temp_band.WriteArray(processed_chunks[i][ :-window_size,window_size:-window_size])  # Skip top
            elif is_bottom_border:
                # Bottom edge: reduce from left, right, and top
                temp_band.WriteArray(processed_chunks[i][window_size:,window_size:-window_size])  # Skip bottom
            else:
                # Non-border chunks: apply window size reduction
                temp_band.WriteArray(processed_chunks[i][window_size:-window_size, window_size:-window_size])
        else:
            temp_band.WriteArray(processed_chunks[i])

        temp_dataset.FlushCache()
        temp_dataset = None

        temp_paths.append(temp_path)

    return temp_paths, x_start, y_start

def merge_temp_files(output_filepaths, temp_files, output_width, output_height, output_geotransform, output_projection, num_outputs,cog_flag,cog_overviews,azlks=1, rglks=1):
    for i in range(num_outputs):
        if '.tif' in output_filepaths[0]:
            driver = gdal.GetDriverByName('GTiff')
            options = ['COMPRESS=LZW', 'BIGTIFF=YES', 'TILED=YES']
            if cog_flag:
                options.extend(['COG=YES', 'TFW=YES']) 
            output_dataset = driver.Create(output_filepaths[i], output_width, output_height, 1, gdal.GDT_Float32, options=options)
        elif '.bin' in output_filepaths[0]:
            driver = gdal.GetDriverByName('ENVI')
            output_dataset = driver.Create(output_filepaths[i], output_width, output_height, 1, gdal.GDT_Float32)

        output_dataset.SetGeoTransform(output_geotransform)
        output_dataset.SetProjection(output_projection)
        output_band = output_dataset.GetRasterBand(1)

        for temp_path_set, x_start, y_start in temp_files:
            temp_dataset = gdal.Open(temp_path_set[i], gdal.GA_ReadOnly)
            temp_band = temp_dataset.GetRasterBand(1)

            temp_chunk = temp_band.ReadAsArray()
            temp_height, temp_width = temp_chunk.shape
            
            # Scale the input-based chunk position to the output grid
            scaled_x = x_start 
            scaled_y = y_start 

            # Safety check: trim the chunk if it goes beyond bounds
            chunk_height, chunk_width = temp_chunk.shape
            if scaled_x + chunk_width > output_width:
                temp_chunk = temp_chunk[:, :output_width - scaled_x]
            if scaled_y + chunk_height > output_height:
                temp_chunk = temp_chunk[:output_height - scaled_y, :]
            
            
            if len(temp_files) == 1:
                x_start, y_start = 0, 0
                scaled_x, scaled_y = 0, 0
                
            # output_band.WriteArray(temp_chunk, xoff=x_start, yoff=y_start)
            output_band.WriteArray(temp_chunk, xoff=scaled_x, yoff=scaled_y)
            # print(x_start,y_start)
            temp_dataset = None

        output_dataset.FlushCache()
        output_dataset = None
        
        if '.tif' in output_filepaths[i] and cog_flag:
            # add_overviews(output_filepaths[i])
            gdal.SetConfigOption('COMPRESS_OVERVIEW', 'LZW')
            dataset = gdal.Open(output_filepaths[i], gdal.GA_Update)
            dataset.BuildOverviews("AVERAGE", cog_overviews)
            dataset = None
        
        print(f"Saved file {output_filepaths[i]}")



def process_and_write_chunk(args, processing_func, num_outputs, output_width, output_height, **proc_kwargs):
    # try:
        (input_filepaths, x_start, y_start, read_block_width, read_block_height, window_size, raster_width, raster_height, *proc_args) = args

        chunks = [read_chunk_with_overlap(fp, x_start, y_start, read_block_width, read_block_height, window_size) for fp in input_filepaths]
        # processed_chunks = processing_func(chunks, window_size, input_filepaths, chi_in,psi_in,model)
        processed_chunks = processing_func(
            # chunks, window_size, input_filepaths, *proc_args, **proc_kwargs
            chunks, window_size, input_filepaths, *proc_args, **proc_kwargs
            # chunks, *proc_args, **proc_kwargs
        )

        if num_outputs == 1:
            processed_chunks = [processed_chunks]

    
        azlks = proc_kwargs.get('azlks', 1)
        rglks = proc_kwargs.get('rglks', 1)

        out_x_start = x_start // rglks
        out_y_start = y_start // azlks

        temp_paths, temp_x_start, temp_y_start = write_chunk_to_temp_file(
            processed_chunks, out_x_start, out_y_start, read_block_width//rglks, read_block_height//azlks, window_size,  output_width, output_height,  num_outputs
        )

        return temp_paths, out_x_start, out_y_start
    # except Exception as e:
    #     print(f"Error in processing chunk: {e}")
    #     return None
