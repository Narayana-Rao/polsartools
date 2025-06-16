import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import conv2d,time_it
from polsartools.preprocessing.pre_utils import get_filter_io_paths

from polsartools.preprocessing.rflee_filter import process_chunk_refined_lee
from polsartools.rflee import process_chunk_rfleecpp
@time_it
def boxcar(infolder,  window_size=3, outType="tif", cog_flag=False, cog_overviews = [2, 4, 8, 16], write_flag=True, max_workers=None,block_size=(512, 512)):
    
    input_filepaths, output_filepaths = get_filter_io_paths(infolder, window_size, outType=outType, filter_type="boxcar")

    # Process chunks in parallel
    num_outputs = len(output_filepaths)

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_boxcar,block_size=block_size, max_workers=max_workers,  num_outputs=num_outputs,
                            cog_flag=cog_flag,
                            cog_overviews=cog_overviews,
                            )

def process_chunk_boxcar(chunks, window_size, *args):

    filtered_chunks = []
    for i in range(len(chunks)):
        img = np.array(chunks[i])
        kernel = np.ones((window_size, window_size), np.float32) / (window_size * window_size)
        filtered_chunks.append(conv2d(img, kernel))
    return filtered_chunks


@time_it
def rlee(infolder,  window_size=3, outType="tif", cog_flag=False, cog_overviews = [2, 4, 8, 16], write_flag=True, max_workers=None,block_size=(512, 512)):


    input_filepaths, output_filepaths = get_filter_io_paths(infolder, window_size, outType=outType, filter_type="rlee")
    num_outputs = len(output_filepaths)
    
    ### Python implementation
    # process_chunks_parallel(input_filepaths, output_filepaths, window_size=window_size, write_flag=write_flag,
    #                         processing_func=process_chunk_refined_lee, block_size=(512, 512), max_workers=max_workers,
    #                         num_outputs=num_outputs)
    
    #### Uncomment below to use C++ implementation 

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_rfl,block_size=block_size, max_workers=max_workers,  num_outputs=num_outputs,
                            cog_flag=cog_flag,
                            cog_overviews=cog_overviews,
                            )


def process_chunk_rfl(chunks, window_size, *args):

    # print('before pad',np.shape(chunks[0]))
    
    for i in range(len(chunks)):
        pad_top_left = window_size // 2 
        pad_bottom_right = window_size // 2 +1
        chunks[i] = np.pad(chunks[i], 
                            ((pad_top_left, pad_bottom_right), 
                            (pad_top_left, pad_bottom_right)), 
                            mode='constant', constant_values=0)
    
    
    # print("after pad",np.shape(chunks[0]))
    chunk_arrays = [np.array(ch) for ch in chunks]  
    vi_c_raw = process_chunk_rfleecpp(chunk_arrays, window_size)
    
    proc_chunks=[]
    for chunk in vi_c_raw:
        filt_data = np.array(chunk)
        # filt_data[filt_data == 0] = np.nan
        proc_chunks.append(filt_data)
        # print('mean %0.3f'%np.nanmean(np.array(chunk)),'std %0.3f'%np.nanstd(np.array(chunk)))
        
        
    del vi_c_raw,filt_data,chunk
    
    # print("proc_chunks pad",np.shape(proc_chunks[0]))
    for i in range(len(proc_chunks)):
        # Calculate the padding size
        pad_top_left = window_size // 2 
        pad_bottom_right = window_size // 2 +1
        proc_chunks[i] = proc_chunks[i][pad_top_left:-pad_bottom_right, pad_top_left:-pad_bottom_right]
    
    
    def shift_array(arr,shift):
        # Get the number of rows and columns in the array
        rows, cols = arr.shape

        # Step 1: Move the rightmost 3 columns to the left
        right_columns = arr[:, -shift:]  # Last 3 columns
        remaining_columns = arr[:, :-shift]  # All but last 3 columns

        # Step 2: Move the bottom 3 rows to the top
        bottom_rows = arr[-shift:, :]  # Last 3 rows
        remaining_rows = arr[:-shift, :]  # All but last 3 rows

        # Combine the shifted rows and columns
        shifted_array = np.vstack((bottom_rows, remaining_rows))  # Stack bottom rows to the top
        shifted_array = np.hstack((right_columns, shifted_array[:, :-shift]))  # Stack right columns to the left

        return shifted_array
    
    
    for i in range(len(proc_chunks)):
        # Remove 'window_size' rows from the top and 'window_size' columns from the left
        proc_chunks[i] = shift_array(proc_chunks[i],window_size//2)
    
    
    # print("proc_chunks unpad",np.shape(proc_chunks[0]))

    
    num_chunks = len(proc_chunks) // 2
    out_chunks = []

    for i in range(num_chunks):
        real_part = proc_chunks[2 * i]       # Get the real part from the even indices
        # real_part[real_part == 0] = np.nan
        imag_part = proc_chunks[2 * i + 1]   # Get the imaginary part from the odd indices
        # imag_part[imag_part == 0] = np.nan
        complex_array = real_part + 1j * imag_part  # Create a complex number
        
        out_chunks.append(complex_array)
        # print(np.nanmean(real_part),' ' ,np.nanmean(imag_part))
    filtered_chunks = []

    if len(chunks)==9:
        # print("out_chunks shape:", np.shape(out_chunks))
        filtered_chunks.append(np.real(out_chunks[0]))
        filtered_chunks.append(np.real(out_chunks[1]))
        filtered_chunks.append(np.imag(out_chunks[1]))
        filtered_chunks.append(np.real(out_chunks[2]))
        filtered_chunks.append(np.imag(out_chunks[2]))
        filtered_chunks.append(np.real(out_chunks[4]))
        filtered_chunks.append(np.real(out_chunks[5]))
        filtered_chunks.append(np.imag(out_chunks[5]))
        filtered_chunks.append(np.real(out_chunks[8]))
    if len(chunks)==4:
        filtered_chunks.append(np.real(out_chunks[0]))
        filtered_chunks.append(np.real(out_chunks[1]))
        filtered_chunks.append(np.imag(out_chunks[1]))
        filtered_chunks.append(np.real(out_chunks[3]))
   
    return filtered_chunks
    