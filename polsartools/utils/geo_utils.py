from osgeo import gdal
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import xml.etree.ElementTree as ET
def geocode_grid(input_vrt, output_file, filtype, dx, dy, bbox, dstSRS="EPSG:4326",resampleAlg="near"):
    if filtype=='bin':
        warp_options = gdal.WarpOptions(
            format="ENVI",
            dstSRS=dstSRS,
            geoloc=True,
            xRes=dx,
            yRes=dy,
            resampleAlg=resampleAlg,
            srcNodata=0,
            dstNodata=0,
            outputBounds=(bbox[0],bbox[1],bbox[2],bbox[3]),
        )
    else:
        warp_options = gdal.WarpOptions(
            format="GTiff",
            dstSRS=dstSRS,
            geoloc=True,
            xRes=dx,
            yRes=dy,
            resampleAlg=resampleAlg,
            srcNodata=0,
            dstNodata=0,
            outputBounds=(bbox[0],bbox[1],bbox[2],bbox[3]),
        )
    
    ds = gdal.Warp(
        destNameOrDestDS=output_file,
        srcDSOrSrcDSTab=input_vrt,
        options=warp_options
    )
    ds = None
    

def intp_grid(original_array, target_shape):
    x = np.arange(original_array.shape[0])
    y = np.arange(original_array.shape[1]) 
    interpolator = RegularGridInterpolator((x, y), original_array, method='linear', bounds_error=False, fill_value=None)
    
    x_new = np.linspace(0, original_array.shape[0] - 1, target_shape[0])
    y_new = np.linspace(0, original_array.shape[1] - 1, target_shape[1])

    xx, yy = np.meshgrid(x_new, y_new, indexing='ij')
    target_points = np.stack((xx, yy), axis=-1)

    resized_array = interpolator(target_points)
    return resized_array


def update_vrt(datafile):   
    tree = ET.parse(datafile)
    root = tree.getroot()
    
    meta = ET.SubElement(root, "Metadata", domain="GEOLOCATION")
    for key, value in {
        "X_DATASET": "lon.vrt",
        "Y_DATASET": "lat.vrt",
        "X_BAND": "1",
        "Y_BAND": "1",
        "PIXEL_OFFSET": "0",
        "LINE_OFFSET": "0",
        # "PIXEL_STEP": f"{ml.shape[1]//coordinateX.shape[1]}",
        # "LINE_STEP": f"{ml.shape[0]//coordinateX.shape[0]}"
        
        "PIXEL_STEP": "1",
        "LINE_STEP": "1"
        
        
    }.items():
        ET.SubElement(meta, "MDI", key=key).text = value
    
    tree.write(datafile)

def write_latlon(filename, array):
    rows,cols = array.shape
    driver = gdal.GetDriverByName('ENVI')  # or 'GTiff'
    ds = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
    ds.GetRasterBand(1).WriteArray(array)
    ds = None
