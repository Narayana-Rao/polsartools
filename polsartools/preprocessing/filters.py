import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d
from polsartools.preprocessing.pre_utils import get_filter_io_paths

# from polsartools.preprocessing.rflee_filter import process_chunk_refined_lee
from polsartools.rflee import process_chunk_rfleecpp
@time_it
def boxcar(infolder, outname=None, chi_in=0, psi_in=0, window_size=3, write_flag=True, max_workers=None):
    
    input_filepaths, output_filepaths = get_filter_io_paths(infolder, outname, window_size, filter_type="boxcar")

    # Process chunks in parallel
    num_outputs = len(output_filepaths)
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_boxcar, block_size=(512, 512), max_workers=max_workers, 
                            num_outputs=num_outputs)

def process_chunk_boxcar(chunks, window_size, input_filepaths, *args):
    filtered_chunks = []
    for i in range(len(chunks)):
        img = np.array(chunks[i])
        kernel = np.ones((window_size, window_size), np.float32) / (window_size * window_size)
        filtered_chunks.append(conv2d(img, kernel))
    return filtered_chunks


@time_it
def rlee(infolder, outname=None, chi_in=0, psi_in=0, window_size=3, write_flag=True, max_workers=None):
    # File reading and output setup similar to the boxcar filter
    input_filepaths, output_filepaths = get_filter_io_paths(infolder, outname, window_size, filter_type="refined_lee")

    # Process chunks in parallel
    num_outputs = len(output_filepaths)
    # process_chunks_parallel(input_filepaths, output_filepaths, window_size=window_size, write_flag=write_flag,
    #                         processing_func=process_chunk_refined_lee, block_size=(512, 512), max_workers=max_workers,
    #                         num_outputs=num_outputs)
    
    process_chunks_parallel(input_filepaths, output_filepaths, window_size=window_size, write_flag=write_flag,
                        processing_func=process_chunk_rfl, block_size=(512, 512), max_workers=max_workers,
                        num_outputs=num_outputs)

def process_chunk_rfl(chunks, window_size,input_filepaths, *args):

    chunk_arrays = [np.array(ch) for ch in chunks]  
    print("chunk_arrays shape:", np.shape(chunk_arrays))
    vi_c_raw = process_chunk_rfleecpp(chunk_arrays, window_size)
    
    proc_chunks=[]
    for chunk in vi_c_raw:
        proc_chunks.append(np.array(chunk))
    print("vi_c_raw shape:", np.shape(vi_c_raw))
    # print("vi_c_raw len:", type(vi_c_raw[0]))
    # print("vi_c_raw len:", len(vi_c_raw[0]))
    # print("vi_c_raw np.shape(:", np.shape(vi_c_raw[0]))
    print(np.nanmean(proc_chunks[0]), np.nanstd(proc_chunks[0]), np.nanmin(proc_chunks[0]), np.nanmax(proc_chunks[0]))
    return proc_chunks
    