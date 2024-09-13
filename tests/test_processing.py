import os
import pytest
import numpy as np
import rasterio
from osgeo import gdal
from polsartools.polsar.fullpol.processing import process_raster_in_chunks

@pytest.fixture
def create_dummy_geotiff():
    """
    Fixture to create a dummy single-band GeoTIFF file for testing.
    """
    directory = 'test_data'
    filename = 'test_raster.tif'
    filepath = os.path.join(directory, filename)
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Define data type mapping
    dtype = 'float32'
    dtype_map = {
        'float32': rasterio.float32,
        'int16': rasterio.int16,
        'uint16': rasterio.uint16,
        # Add more mappings if needed
    }
    
    if dtype not in dtype_map:
        raise ValueError(f"Unsupported dtype: {dtype}")

    # Create a dummy GeoTIFF file
    with rasterio.open(filepath, 'w', driver='GTiff', height=512, width=512, count=1, dtype=dtype_map[dtype], crs='EPSG:4326', transform=rasterio.Affine(1, 0, 0, 0, -1, 0)) as dst:
        data = np.random.random((512, 512)).astype(dtype)
        dst.write(data, 1)
    
    yield filepath

    # Clean up
    os.remove(filepath)
    os.rmdir(directory)

def test_process_raster_in_chunks(create_dummy_geotiff):
    """
    Test processing of raster file in chunks.
    """
    input_file = create_dummy_geotiff
    output_file = os.path.join('test_data', 'processed_raster.tif')
    param1 = 2  # Example parameter for processing
    param2 = 5  # Example parameter for processing
    chunk_size = (256, 256)  # Smaller chunks for testing

    process_raster_in_chunks(input_file, output_file, param1, param2, chunk_size)

    # Verify the output file
    with rasterio.open(output_file) as dst:
        assert dst is not None, "Output file could not be opened"
        assert dst.count == 1  # Verify single-band
        assert dst.width == 512
        assert dst.height == 512
        assert dst.dtypes[0] == 'float32'  # Ensure dtype matches

    # Clean up
    os.remove(output_file)
