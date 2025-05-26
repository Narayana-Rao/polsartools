import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import time_it
from polsartools.cprvicpp import process_chunk_cprvicpp
from .cp_infiles import cpc2files
@time_it
def cprvi(infolder, outname=None, chi_in=45, psi_in=0,window_size=1,write_flag=True,max_workers=None):
    input_filepaths = cpc2files(infolder)
    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(infolder, "cprvi.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, 
                write_flag=write_flag,processing_func=process_chunk_cprvi,block_size=(512, 512), 
                max_workers=max_workers,  num_outputs=1, chi_in=chi_in,psi_in=psi_in)

def process_chunk_cprvi(chunks, window_size,chi_in,psi_in):

    chunk_arrays = [np.array(ch) for ch in chunks]  
    # CPP function
    vi_c_raw = process_chunk_cprvicpp( chunk_arrays, window_size, chi_in, psi_in )

    return np.array(vi_c_raw, copy=True)  
    
 