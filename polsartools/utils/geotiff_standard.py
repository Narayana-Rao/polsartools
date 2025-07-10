from osgeo import gdal
gdal.UseExceptions()
import numpy as np

def save_polarimetric_geotiff(output_path, matrices, metadata, x_size, y_size, projection, geo_transform):
    """
    Saves polarimetric data (S2/C3/C4/T3/T4/T2/C2) into a multi-band GeoTIFF with metadata.
    
    Args:
        output_path (str): Path to save the GeoTIFF file.
        matrices (dict): Dictionary with band names and corresponding numpy arrays.
        metadata (dict): Dictionary containing metadata (sensor, acquisition date, matrix type, etc.).
        x_size, y_size (int): Image dimensions.
        projection (str): Projection string (e.g., EPSG:4326).
        geo_transform (tuple): Georeferencing tuple.

    Saves:
        - Multi-band GeoTIFF with real & imaginary components for complex data
        - Metadata including matrix type flag and band names
    """

    driver = gdal.GetDriverByName("GTiff")
    creation_options = ["TILED=YES", "COPY_SRC_OVERVIEWS=YES", "COMPRESS=LZW"]
    dataset = driver.Create(output_path, x_size, y_size, len(matrices), gdal.GDT_Float32,
                           options=creation_options)

    dataset.SetProjection(projection)
    dataset.SetGeoTransform(geo_transform)

    band_names = []
    for idx, (band_name, matrix) in enumerate(matrices.items(), start=1):
        band = dataset.GetRasterBand(idx)
        band.WriteArray(matrix)
        band.SetDescription(band_name)  # Store band name
        band_names.append(band_name)

    # Adding metadata with placeholders
    dataset.SetMetadata({
        "Sensor": metadata.get("sensor", "Unknown"),
        "Date_Acquired": metadata.get("date_acquired", "N/A"),
        "Matrix_Type": metadata.get("matrix_type", "Unknown"),
        "Polarization_Mode": metadata.get("polarization", "N/A"),
        "Processing_Steps": metadata.get("processing_steps", "N/A"),
        "Azimuth_Looks": str(metadata.get("azlks", "N/A")),
        "Range_Looks": str(metadata.get("rglks", "N/A")),
        "Band_Names": ", ".join(band_names)  # Store band names as metadata
    })

    dataset.FlushCache()
    print(f"Polarimetric GeoTIFF saved as {output_path}")

# Example Call for S2
save_polarimetric_geotiff(
    "s2_polarimetric_data.tif",
    {
        "S11_real": np.random.rand(100, 100),
        "S11_imag": np.random.rand(100, 100),
        "S12_real": np.random.rand(100, 100),
        "S12_imag": np.random.rand(100, 100),
        "S21_real": np.random.rand(100, 100),
        "S21_imag": np.random.rand(100, 100),
        "S22_real": np.random.rand(100, 100),
        "S22_imag": np.random.rand(100, 100),
    },
    {
        "sensor": "NISAR",
        "date_acquired": "2025-06-10",
        "matrix_type": "S2",
        "polarization": "monostatic",
        "processing_steps": "Multi-look applied",
        "azlks": 3,
        "rglks": 4
    },
    100, 100, "EPSG:4326", (0, 10, 0, 0, 0, -10)
)