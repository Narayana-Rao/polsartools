import numpy as np
from osgeo import gdal, osr

def write_raster(file_path, data_array):
    """
    Writes an Xarray DataArray to a raster file.
    """
    driver = gdal.GetDriverByName('GTiff')
    rows, cols = data_array.shape
    ds = driver.Create(file_path, cols, rows, 1, gdal.GDT_Float32)
    band = ds.GetRasterBand(1)
    band.WriteArray(data_array.values)
    ds.FlushCache()
