import os
import numpy as np
import tables
from osgeo import gdal, osr

tables.parameters.PYTABLES_SYS_ATTRS = False

def process_nisar_soil_moisture(h5_path, layer_to_extract=None, output_tif_path=None):
    """
    Extracts Soil Moisture (SM) or Soil Moisture Uncertainty (SMU) from a NISAR SME2 HDF5 file 
    using PyTables, and converts it into a georeferenced GeoTIFF using GDAL.
    """
    layer_mapping = {
        "SM": {
            "name": "Soil Moisture",
            "path": "/science/LSAR/SME2/grids/soilMoisture"
        },
        "SMU": {
            "name": "Soil Moisture Uncertainty",
            "path": "/science/LSAR/SME2/grids/soilMoistureUncertainty"
        }
    }
    
    if not layer_to_extract:
        print("==================================================")
        print("     NISAR SME2 Soil Moisture Raster Extractor   ")
        print("==================================================")
        print("No layer was specified. Please choose an option:")
        print("  • Enter layer_to_extract='SM'  -> Extract Soil Moisture")
        print("  • Enter layer_to_extract='SMU' -> Extract Soil Moisture Uncertainty")
        print("==================================================")
        return

    layer_key = layer_to_extract.upper()
    if layer_key not in layer_mapping:
        raise ValueError(f"Invalid layer '{layer_to_extract}'. Please use 'SM' for Soil Moisture or 'SMU' for Uncertainty.")
        
    target_info = layer_mapping[layer_key]
    target_node_path = target_info["path"]
    
    in_dir = os.path.dirname(h5_path) or "."
    base_name = os.path.basename(h5_path).split('.h5')[0]
    out_dir = os.path.join(in_dir, base_name)
    os.makedirs(out_dir, exist_ok=True)
    
    if not output_tif_path:
        output_tif_path = os.path.join(out_dir, f"{layer_key.lower()}_output.tif")
    else:
        if not os.path.isabs(output_tif_path):
            output_tif_path = os.path.join(out_dir, output_tif_path)

    # Open H5 file using PyTables
    with tables.open_file(h5_path, mode='r') as h5file:
        try:
            print(f"\n[INFO] Reading metadata and {target_info['name']} ({layer_key})...")
            
            grid_base = "/science/LSAR/SME2/grids"
            xcoords = h5file.get_node(f"{grid_base}/xCoordinates").read()
            ycoords = h5file.get_node(f"{grid_base}/yCoordinates").read()
            x_spacing = float(h5file.get_node(f"{grid_base}/xCoordinateSpacing").read())
            y_spacing = float(h5file.get_node(f"{grid_base}/yCoordinateSpacing").read())
            projection = int(h5file.get_node(f"{grid_base}/projection").read())
            
            data_node = h5file.get_node(target_node_path)
            data_array = data_node.read()
            
        except tables.NoSuchNodeError as e:
            raise KeyError(f"Required node/dataset missing in the SME2 file: {e}")
        except Exception as e:
            raise RuntimeError(f"Error reading HDF5 dataset or metadata: {e}")
            
    if data_array.ndim > 2:
        data_array = np.squeeze(data_array)
        
    # ==========================================
    # Convert to Georeferenced GeoTIFF via GDAL
    # ==========================================
    height, width = data_array.shape
    
    # Upper-Left Coordinates
    # Half-pixel shift to account for pixel-center coordinates to pixel-corner bounding box
    start_x = np.min(xcoords) - (x_spacing / 2.0)
    start_y = np.max(ycoords) + (abs(y_spacing) / 2.0)
    
    # Affine GeoTransform: (Top Left X, W-E pixel resolution, Rotation 0, Top Left Y, Rotation 0, N-S pixel resolution [Negative])
    geotransform = (start_x, abs(x_spacing), 0.0, start_y, 0.0, -abs(y_spacing))
    
    # Map NumPy data types to GDAL data types
    gdal_datatype_map = {
        np.dtype('float32'): gdal.GDT_Float32,
        np.dtype('float64'): gdal.GDT_Float64,
        np.dtype('int16'): gdal.GDT_Int16,
        np.dtype('uint16'): gdal.GDT_UInt16,
        np.dtype('int32'): gdal.GDT_Int32,
        np.dtype('uint32'): gdal.GDT_UInt32,
        np.dtype('uint8'): gdal.GDT_Byte
    }
    gdal_dtype = gdal_datatype_map[data_array.dtype]

    # Initialize GeoTIFF driver and options
    driver = gdal.GetDriverByName('GTiff')
    options = ['COMPRESS=LZW', 'TILED=YES']
    out_ds = driver.Create(output_tif_path, width, height, 1, gdal_dtype, options)
    
    # Set spatial reference (CRS)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(projection)
    
    # Enable modern Axis Mapping Strategy to ensure correct EPSG coordinates order
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    out_ds.SetProjection(srs.ExportToWkt())
    out_ds.SetGeoTransform(geotransform)
    
    # Write array band and release dataset
    print(f"[INFO] Writing layer to GeoTIFF via GDAL: {output_tif_path}")
    band = out_ds.GetRasterBand(1)
    
    # Handle NaN fill values if present
    if np.issubdtype(data_array.dtype, np.floating):
        band.SetNoDataValue(-9999)
        
    band.WriteArray(data_array)
    band.FlushCache()
    
    # Dereference dataset to force flush write and close file safely
    out_ds = None
    band = None

if __name__ == "__main__":
    h5_file_path = "/media/hm/63740439-1587-4f23-a0bf-254c0b8985c7/NISAR_sample_data/NISAR_L3_PR_SME2_005_156_A_011_4005_DHDH_A_20251121T000949_20251121T001024_X05010_N_F_J_001.h5"
    
    process_nisar_soil_moisture(
        h5_file_path, 
        layer_to_extract="SMU", 
        output_tif_path="soil_moisture_uncertainty.tif"
    )