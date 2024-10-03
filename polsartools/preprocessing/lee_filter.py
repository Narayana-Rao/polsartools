import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d

@time_it
def lee_filter(infolder, outname=None, window_size=3, write_flag=True, max_workers=None):
    # File reading and output setup similar to the boxcar filter
    # Identical file structure to boxcar, will use same input/output logic
    input_filepaths, output_filepaths = get_input_output_filepaths(infolder, outname, window_size, filter_type="lee")

    # Process chunks in parallel
    num_outputs = len(output_filepaths)
    process_chunks_parallel(input_filepaths, output_filepaths, window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_lee, block_size=(512, 512), max_workers=max_workers,
                            num_outputs=num_outputs)

def process_chunk_lee(chunks, window_size, *args):
    """Process Lee filtering on the chunk."""
    filtered_chunks = []
    for i in range(len(chunks)):
        img = np.array(chunks[i])
        mean_kernel = np.ones((window_size, window_size)) / (window_size * window_size)
        
        local_mean = conv2d(img, mean_kernel)
        local_var = conv2d(img ** 2, mean_kernel) - local_mean ** 2
        noise_var = np.mean(local_var)  # Assumed noise variance (can be adjusted based on actual noise)
        
        lee_filtered = local_mean + (local_var - noise_var) / (local_var + noise_var) * (img - local_mean)
        filtered_chunks.append(lee_filtered)
    
    return filtered_chunks
