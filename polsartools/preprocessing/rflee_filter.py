import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d

@time_it
def refined_lee_filter(infolder, outname=None, window_size=3, k=1, write_flag=True, max_workers=None):
    # File reading and output setup similar to the boxcar filter
    input_filepaths, output_filepaths = get_input_output_filepaths(infolder, outname, window_size, filter_type="refined_lee")

    # Process chunks in parallel
    num_outputs = len(output_filepaths)
    process_chunks_parallel(input_filepaths, output_filepaths, window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_refined_lee, block_size=(512, 512), max_workers=max_workers,
                            num_outputs=num_outputs)

def process_chunk_refined_lee(chunks, window_size, k=1, *args):
    """Process Refined Lee filtering on the chunk."""
    filtered_chunks = []
    for i in range(len(chunks)):
        img = np.array(chunks[i])
        mean_kernel = np.ones((window_size, window_size)) / (window_size * window_size)
        
        local_mean = conv2d(img, mean_kernel)
        local_var = conv2d(img ** 2, mean_kernel) - local_mean ** 2
        coeff_var = np.sqrt(local_var) / (local_mean + 1e-6)

        # Classify the pixel based on coefficient of variation
        smooth_weight = np.exp(-k * coeff_var)
        refined_lee_filtered = smooth_weight * local_mean + (1 - smooth_weight) * img
        filtered_chunks.append(refined_lee_filtered)
    
    return filtered_chunks
