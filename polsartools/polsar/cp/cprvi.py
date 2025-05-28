import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import time_it
from polsartools.cprvicpp import process_chunk_cprvicpp
from .cp_infiles import cpc2files
@time_it
def cprvi(infolder,   chi_in=45, psi_in=0, window_size=1, outType="tif", cog_flag=False, 
          cog_overviews = [2, 4, 8, 16], write_flag=True, 
          max_workers=None,block_size=(512, 512)):
    input_filepaths = cpc2files(infolder)
    output_filepaths = []
    if outType == "bin":
        output_filepaths.append(os.path.join(infolder, "cprvi.bin"))
    else:
        output_filepaths.append(os.path.join(infolder, "cprvi.tif"))

    process_chunks_parallel(input_filepaths, list(output_filepaths), 
                            window_size,
                        write_flag,
                        process_chunk_cprvi,
                        *[chi_in, psi_in],
                        block_size=block_size, 
                        max_workers=max_workers,  
                        num_outputs=len(output_filepaths),
                        cog_flag=cog_flag,
                        cog_overviews=cog_overviews,

                        )
def process_chunk_cprvi(chunks, window_size, *args, **kwargs):
    
    chi_in=args[-2]
    psi_in=args[-1]
    # print(chi_in,psi_in)
    
    chunk_arrays = [np.array(ch) for ch in chunks]  
    # CPP function
    vi_c_raw = process_chunk_cprvicpp( chunk_arrays, window_size, chi_in, psi_in )

    return np.array(vi_c_raw, copy=True)  
    
 