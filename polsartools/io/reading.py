import os
import glob
import xarray as xr
from osgeo import gdal

def read_geotiffs(folder_path):
    """
    Loads all GeoTIFF files in the specified folder into Xarrays.
    
    Parameters:
        folder_path (str): Path to the folder containing GeoTIFF files.
    
    Returns:
        list of xarray.DataArray: List of DataArrays containing data from each GeoTIFF.
    """
    tif_files = glob.glob(os.path.join(folder_path, '*.tif'))
    
    if not tif_files:
        raise FileNotFoundError("No GeoTIFF files found in the specified folder.")
    
    data_arrays = []
    
    for tif_file in tif_files:
        ds = gdal.Open(tif_file)
        band = ds.GetRasterBand(1)
        data = band.ReadAsArray()
        data_array = xr.DataArray(data, dims=['y', 'x'])
        data_arrays.append(data_array)
    
    return data_arrays, tif_files
