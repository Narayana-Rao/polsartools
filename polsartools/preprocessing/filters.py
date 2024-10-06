import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d
from polsartools.preprocessing.pre_utils import get_filter_io_paths

from polsartools.preprocessing.rflee_filter import process_chunk_refined_lee

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
    process_chunks_parallel(input_filepaths, output_filepaths, window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_refined_lee, block_size=(512, 512), max_workers=max_workers,
                            num_outputs=num_outputs)
