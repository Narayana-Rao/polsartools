import gdal
import xarray as xr

def read_raster(file_path):
    """
    Reads a raster file and returns an Xarray DataArray.
    """
    ds = gdal.Open(file_path)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    return xr.DataArray(data, dims=['y', 'x'])
