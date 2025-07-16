import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d
from polsartools.preprocess.pre_utils import get_filter_io_paths

@time_it
def gamma_filter(infolder, outname=None, window_size=3, enl=1, write_flag=True, max_workers=None):
    # File reading and output setup similar to the boxcar filter
    input_filepaths, output_filepaths = get_filter_io_paths(infolder, outname, window_size, filter_type="gamma")

    # Process chunks in parallel
    num_outputs = len(output_filepaths)
    process_chunks_parallel(input_filepaths, output_filepaths, window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_gamma, block_size=(512, 512), max_workers=max_workers,
                            num_outputs=num_outputs)

def process_chunk_gamma(chunks, window_size, enl=1, *args):
    """Process Gamma MAP filtering on the chunk."""
    filtered_chunks = []
    for i in range(len(chunks)):
        img = np.array(chunks[i])
        mean_kernel = np.ones((window_size, window_size)) / (window_size * window_size)
        
        local_mean = conv2d(img, mean_kernel)
        local_var = conv2d(img ** 2, mean_kernel) - local_mean ** 2
        gamma_filtered = (local_mean + local_var / (enl + local_var)) * (img - local_mean)
        filtered_chunks.append(gamma_filtered)
    
    return filtered_chunks
