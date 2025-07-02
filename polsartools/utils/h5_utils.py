# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 11:44:43 2025

@author: nbhogapurapu
"""

import os
import glob
import numpy as np
import tables
from skimage.util.shape import view_as_blocks
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import sys
from netCDF4 import Dataset
from osgeo import gdal, osr
gdal.UseExceptions()


import warnings
warnings.filterwarnings("ignore", category=tables.exceptions.DataTypeWarning)
warnings.filterwarnings(action='ignore', message='Mean of empty slice')

def get_dtype(numpy_dtype):
    """Map numpy dtype to GDAL data type."""
    dtype = np.dtype(numpy_dtype)  # Ensure it's a dtype object
    return {
        np.uint8: gdal.GDT_Byte,
        np.int16: gdal.GDT_Int16,
        np.uint16: gdal.GDT_UInt16,
        np.int32: gdal.GDT_Int32,
        np.uint32: gdal.GDT_UInt32,
        np.float32: gdal.GDT_Float32,
        np.float64: gdal.GDT_Float64,
        np.complex64: gdal.GDT_CFloat32,
        np.complex128: gdal.GDT_CFloat64
    }.get(dtype.type, gdal.GDT_Unknown)
def mlook_arr(data, az, rg):
    temp = data[0:data.shape[0] - data.shape[0] % az,
                0:data.shape[1] - data.shape[1] % rg]
    blocks = view_as_blocks(temp, block_shape=(az, rg))
    flatten = blocks.reshape(blocks.shape[0], blocks.shape[1], -1)
    return np.nanmean(flatten, axis=2)



def compute_elements(chunks, matrix_type, azlks, rglks, apply_multilook):
    if matrix_type == "C3":
        return compute_c3(chunks, azlks, rglks, apply_multilook)
    elif matrix_type == "S2":
        return compute_s2(chunks)  # You’ll define this next
    elif matrix_type == "C2":
        return compute_c2(chunks, azlks, rglks, apply_multilook)
    elif matrix_type == "T3":
        return compute_t3(chunks, azlks, rglks, apply_multilook)
    elif matrix_type == "T4":
        return compute_t4(chunks, azlks, rglks, apply_multilook)
    elif matrix_type == "C4":
        return compute_c4(chunks, azlks, rglks, apply_multilook)
    else:
        raise ValueError(f"Unsupported matrix type: {matrix_type}")
        
def compute_c3(chunks, azlks, rglks, apply_multilook):
    def opt_mlook(data):
        return mlook_arr(data, azlks, rglks) if apply_multilook else data

    Kl = np.array([chunks["HH"], np.sqrt(2)*0.5*(chunks["HV"]+chunks["VH"]), chunks["VV"]])

    return {
        "C11": opt_mlook(np.real(np.abs(Kl[0])**2)),
        "C12_real": opt_mlook(np.real(Kl[0]*np.conj(Kl[1]))),
        "C12_imag": opt_mlook(np.imag(Kl[0]*np.conj(Kl[1]))),
        "C13_real": opt_mlook(np.real(Kl[0]*np.conj(Kl[2]))),
        "C13_imag": opt_mlook(np.imag(Kl[0]*np.conj(Kl[2]))),
        "C22": opt_mlook(np.real(np.abs(Kl[1])**2)),
        "C23_real": opt_mlook(np.real(Kl[1]*np.conj(Kl[2]))),
        "C23_imag": opt_mlook(np.imag(Kl[1]*np.conj(Kl[2]))),
        "C33": opt_mlook(np.real(np.abs(Kl[2])**2))
    }


def get_chunk_jobs(h5_file, dataset_path, chunk_size_x, chunk_size_y):
    import tables
    jobs = []
    with tables.open_file(h5_file, 'r') as h5:
        arr = h5.get_node(dataset_path)
        height, width = arr.shape
        for y in range(0, height, chunk_size_y):
            for x in range(0, width, chunk_size_x):
                jobs.append({
                    "x_start": x,
                    "x_end": min(x + chunk_size_x, width),
                    "y_start": y,
                    "y_end": min(y + chunk_size_y, height)
                })
    return jobs



def process_and_write_tile(job, h5_file, dataset_paths,  azlks, rglks, matrix_type,apply_multilook, temp_dir,
                           start_x, start_y, xres, yres, epsg):
    os.makedirs(temp_dir, exist_ok=True)
    x0, x1, y0, y1 = job["x_start"], job["x_end"], job["y_start"], job["y_end"]

    chunks = {}
    with tables.open_file(h5_file, 'r') as h5:
        for name, path in dataset_paths.items():
            chunks[name] = h5.get_node(path)[y0:y1, x0:x1]

    # results = compute_c3_elements(chunks, azlks, rglks, apply_multilook)
    results = compute_elements(chunks, matrix_type, azlks, rglks, apply_multilook)

    for name, data in results.items():
        save_tiff(name, data, x0, y0, azlks, rglks, apply_multilook, temp_dir,
                  start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=epsg)

def save_tiff(name, data, x0, y0, azlks, rglks, apply_multilook, temp_dir,
              start_x=0, start_y=0, xres=1.0, yres=1.0, epsg=4326):


    h_out, w_out = data.shape
    out_path = os.path.join(temp_dir, f"{name}_{y0}_{x0}.tif")
    driver = gdal.GetDriverByName('GTiff')
    dst = driver.Create(out_path, w_out, h_out, 1, gdal.GDT_Float32)

    if not apply_multilook:
        azlks = rglks = 1
    res_x = xres * rglks
    res_y = yres * azlks

    gt_x = start_x + (x0 // rglks) * res_x
    gt_y = start_y - (y0 // azlks) * res_y

    dst.SetGeoTransform([gt_x, res_x, 0, gt_y, 0, -res_y])
    srs = osr.SpatialReference(); srs.ImportFromEPSG(epsg)
    dst.SetProjection(srs.ExportToWkt())
    dst.GetRasterBand(1).WriteArray(data)
    dst.GetRasterBand(1).SetNoDataValue(0)
    dst.SetMetadata({
        'AzimuthLooks': str(azlks),
        'RangeLooks': str(rglks),
        'Multilooked': str(apply_multilook)
    })
    dst = None


def mosaic_chunks(element_name, temp_dir, output_dir, chunk_size_x, chunk_size_y, 
                  azlks, rglks, apply_multilook,
                  start_x=0, start_y=0, xres=1.0, yres=1.0, epsg=4326,
                  outType='tif', dtype=np.float32,
                  inshape=None,outshape=None
                  
                  ):

    dtype = get_dtype(dtype)
    os.makedirs(output_dir, exist_ok=True)
    chunk_files = sorted(glob.glob(os.path.join(temp_dir, f"{element_name}_*.tif")))
    if not chunk_files:
        print(f"No chunks found for {element_name}")
        return

    if not apply_multilook:
        azlks = rglks = 1
        
    res_x = xres * rglks
    res_y = yres * azlks
    
    # print(res_x,res_y)
    
    if apply_multilook and (inshape is not None) and (outshape is not None):
        old_width = inshape[0]
        old_height = inshape[1]
        ml_width = outshape[0]
        ml_height = outshape[1]
        res_x = (xres * old_width) / ml_width
        res_y = (yres * old_height) / ml_height
        
        # print(res_x,res_y)
    

    max_x = 0
    max_y = 0
    offsets = []

    for f in chunk_files:
        y0 = int(f.split("_")[-2])
        x0 = int(f.split("_")[-1].split(".")[0])
        x_ml = x0 // rglks
        y_ml = y0 // azlks
        tile = gdal.Open(f)
        h, w = tile.RasterYSize, tile.RasterXSize
        max_x = max(max_x, x_ml + w)
        max_y = max(max_y, y_ml + h)
        offsets.append((f, x_ml, y_ml))
        tile = None

    driver = gdal.GetDriverByName('ENVI') if outType == 'bin' else gdal.GetDriverByName('GTiff')
    options = ['COMPRESS=LZW', 'BIGTIFF=YES', 'TILED=YES'] if outType == 'tif' else []
    out_path = os.path.join(output_dir, f"{element_name}.{outType}")
    dst = driver.Create(out_path, max_x, max_y, 1, dtype, options=options)

    dst.SetGeoTransform([start_x, res_x, 0, start_y, 0, -abs(res_y)])
   
    # print(start_x, start_y, res_x, -abs(res_y))
    srs = osr.SpatialReference(); srs.ImportFromEPSG(epsg)
    dst.SetProjection(srs.ExportToWkt())

    for f, x_ml, y_ml in offsets:
        tile = gdal.Open(f)
        data = tile.GetRasterBand(1).ReadAsArray()
        
        
        dst.GetRasterBand(1).WriteArray(data, x_ml, y_ml)
        tile = None

    dst.GetRasterBand(1).SetNoDataValue(0)
    dst.SetMetadata({
        'AzimuthLooks': str(azlks),
        'RangeLooks': str(rglks),
        'Multilooked': str(apply_multilook)
    })
    dst = None
    print(f"Saved file {out_path}")


def cleanup_temp_files(temp_dir):
    for f in glob.glob(os.path.join(temp_dir, "*.tif")):
        os.remove(f)
    
    os.rmdir(temp_dir)


def h5_polsar(h5_file, dataset_paths, output_dir, temp_dir,
                            azlks, rglks,matrix_type, apply_multilook,
                            chunk_size_x=512, chunk_size_y=512,
                            max_workers=8,
                            start_x=None, start_y=None, xres=1.0, yres=1.0, epsg=4326,
                            outType='tif',dtype=np.float32,
                            inshape=None,outshape=None):

    jobs = get_chunk_jobs(h5_file, dataset_paths["HH"], chunk_size_x, chunk_size_y)

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(process_and_write_tile, job, h5_file, dataset_paths,
                               azlks, rglks, matrix_type, apply_multilook, temp_dir, 
                               start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=epsg) for job in jobs]
        with tqdm(total=len(futures), desc="Processing chunks") as pbar:
            for _ in as_completed(futures):
                pbar.update(1)

    # Discover output band names
    dummy_shape = (azlks * 2, rglks * 2)
    dummy_data = {k: np.zeros(dummy_shape, dtype=np.complex64) for k in dataset_paths.keys()}
    keys = compute_elements(dummy_data, matrix_type, azlks, rglks, apply_multilook).keys()

    for name in keys:
        mosaic_chunks(name, temp_dir, output_dir, chunk_size_x, chunk_size_y,
                      azlks, rglks, apply_multilook, 
                      start_x, start_y, xres, yres, epsg, 
                      outType, dtype,inshape,outshape)

    cleanup_temp_files(temp_dir)



if __name__ == "__main__":
    inFile = r"C:\Users\nbhogapurapu\Desktop\temp\pstdata\NISAR\RSLC_QP.h5"
    freq_band  ='L'
    output_dir = r'C:\Users\nbhogapurapu\Desktop\temp\pstdata\NISAR\RSLC_QP_tab\C3'
    temp_dir = r'C:\Users\nbhogapurapu\Desktop\temp\pstdata\NISAR\RSLC_QP_tab\C3\temp'
    azlks = 22
    rglks = 10
    h5_polsar(
        h5_file=inFile,
        dataset_paths={
            "HH": "/science/LSAR/RSLC/swaths/frequencyA/HH",
            "HV": "/science/LSAR/RSLC/swaths/frequencyA/HV",
            "VH": "/science/LSAR/RSLC/swaths/frequencyA/VH",
            "VV": "/science/LSAR/RSLC/swaths/frequencyA/VV"
        },
        output_dir=output_dir,
        temp_dir=temp_dir,
        azlks=azlks,
        rglks=rglks,
        matrix_type = 'C3',
        apply_multilook=True,
        chunk_size_x=500,
        chunk_size_y=528,
        max_workers=10,
        # start_x=None, start_y=None, xres=0.0002/rglks, yres=0.0002/azlks, epsg=4326,
        start_x=0, start_y=0, xres=1/rglks, yres=1/azlks, epsg=4326,
        outType='tif',
        dtype = np.float32
        
    )
    
    inFile = r"C:\Users\nbhogapurapu\Desktop\temp\pstdata\NISAR\GSLC_QP.h5"
    freq_band  ='L'
    output_dir = r'C:\Users\nbhogapurapu\Desktop\temp\pstdata\NISAR\GSLC_QP_tab\C3'
    temp_dir = r'C:\Users\nbhogapurapu\Desktop\temp\pstdata\NISAR\GSLC_QP_tab\C3\temp'
    azlks = 2
    rglks = 2
    
    
    with tables.open_file(inFile, "r") as h5File:
        base_path = f'/science/{freq_band}SAR/GSLC/grids/frequencyA'
        projection = np.array(h5File.get_node(f'/science/{freq_band}SAR/GSLC/metadata/radarGrid/projection'))
        xCoordinateSpacing = np.array(h5File.get_node(base_path + '/xCoordinateSpacing'))
        yCoordinateSpacing = np.array(h5File.get_node(base_path + '/yCoordinateSpacing'))
    
    
    ds = Dataset(inFile, "r")
    xcoords = ds.groups["science"].groups[f"{freq_band}SAR"].groups["GSLC"] \
                   .groups["grids"].groups["frequencyA"] \
                   .variables["xCoordinates"][:]
    ycoords = ds.groups["science"].groups[f"{freq_band}SAR"].groups["GSLC"] \
                   .groups["grids"].groups["frequencyA"] \
                   .variables["yCoordinates"][:]
    
    inshape = [len(xcoords),len(ycoords)]
    outshape = [len(xcoords)//rglks,len(ycoords)//azlks]
    
    start_x = min(xcoords)
    start_y = max(ycoords)
    
    
    xres = xCoordinateSpacing[0]
    yres = yCoordinateSpacing[0]
    projection =int( projection[0])
    ds=None
    # print(projection,start_x,start_y,xres,yres,inshape,outshape)
    # print(min(ycoords),max(ycoords), min(ycoords)+np.abs(yres)*(inshape[1]-1))
    h5_polsar(
        h5_file=inFile,
        dataset_paths={
            "HH": f"{base_path}/HH",
            "HV": f"{base_path}/HV",
            "VH": f"{base_path}/VH",
            "VV": f"{base_path}/VV"
            
            
        },
        output_dir=output_dir,
        temp_dir=temp_dir,
        azlks=azlks,
        rglks=rglks,
        matrix_type = 'C3',
        apply_multilook=True,
        chunk_size_x=512,
        chunk_size_y=512,
        max_workers=10,
        start_x=start_x, start_y=start_y, xres=xres, yres=yres, epsg=projection,
        outType='tif',
        dtype = np.float32,
        inshape=inshape,
        outshape=outshape
        
    )
    
    
