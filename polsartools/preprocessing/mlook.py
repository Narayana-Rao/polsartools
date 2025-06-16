import os
import numpy as np
from polsartools.utils.proc_utils import process_chunks_parallel
from polsartools.utils.utils import time_it, mlook_arr
from polsartools.preprocessing.pre_utils import get_filter_io_paths
from osgeo import gdal
@time_it
def mlook(infolder,  azlks=2, rglks=2, outType="tif", cog_flag=False, cog_overviews = [2, 4, 8, 16], write_flag=True, max_workers=None,block_size=(512, 512)):
    
    window_size = azlks
    input_filepaths, output_filepaths = get_filter_io_paths(infolder, window_size, 
                                                            outType=outType, 
                                                            filter_type="ml")


    window_size=None

    """
    GET MULTI-LOOKED RASTER PROPERTIES
       
    """
    
    dataset = gdal.Open(input_filepaths[0], gdal.GA_ReadOnly)
    if dataset is None:
        raise FileNotFoundError(f"Cannot open {input_filepaths[0]}")

    in_cols = dataset.RasterXSize
    in_rows = dataset.RasterYSize
    in_geotransform = dataset.GetGeoTransform()
    in_projection = dataset.GetProjection()

    # Calculate output size after multilooking
    out_x_size = in_cols // rglks
    out_y_size = in_rows // azlks

    # Calculate new geotransform with updated pixel size
    out_geotransform = list(in_geotransform)
    
    if in_geotransform[0]==0.0 and in_geotransform[1]==1.0 and in_geotransform[2]==-0.0 and in_geotransform[3]==0.0 and in_geotransform[4]==-0.0 and in_geotransform[5]==-1.0:
        out_geotransform[1] *= 1 
        out_geotransform[5] *= 1 
    else:
        out_geotransform[1] *= rglks 
        out_geotransform[5] *= azlks
        
    out_geotransform = tuple(out_geotransform)  # back to immutable
    # print(in_geotransform,out_geotransform)
    dataset = None  # Close GDAL dataset
    


    def closest_multiple(value, base):
        return round(value / base) * base

    # Calculate closest multiples of rglks and azlks to 512
    block_x_size = closest_multiple(block_size[0] , rglks)
    block_y_size = closest_multiple(block_size[1] , azlks)

    # print(block_x_size,block_y_size)
    block_size = (block_x_size, block_y_size)




    # Process chunks in parallel
    num_outputs = len(output_filepaths)

    process_chunks_parallel(input_filepaths, list(output_filepaths), 
                             window_size,
                            write_flag,
                            process_chunk_mlook,
                            *[azlks, rglks],
                            block_size=block_size, 
                            max_workers=max_workers,  
                            num_outputs=len(output_filepaths),
                            cog_flag=cog_flag,
                            cog_overviews=cog_overviews,
                            out_x_size=out_x_size,
                            out_y_size=out_y_size,
                            out_geotransform=out_geotransform,
                            out_projection=in_projection,
                            azlks=azlks,
                            rglks=rglks
                            )
    # process_chunks_parallel(input_filepaths, list(output_filepaths), 
    #                         window_size=window_size, 
    #                         write_flag=write_flag,
    #                         processing_func=process_chunk_mlook,
    #                         block_size=block_size, max_workers=max_workers,  
    #                         num_outputs=num_outputs,
    #                         cog_flag=cog_flag,
    #                         cog_overviews=cog_overviews,
    #                         )

def process_chunk_mlook(chunks, window_size, *args, **kwargs):
    azlks=args[-2]
    rglks=args[-1]
    # print("mlook",azlks,rglks)
    mlook_chunks = []
    for i in range(len(chunks)):
        mlook_chunks.append(mlook_arr(np.array(chunks[i]), azlks,rglks).astype(np.float32))

    return mlook_chunks