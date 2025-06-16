import os
import numpy as np
from osgeo import gdal,osr
import matplotlib.pyplot as plt
from polsartools.utils.utils import conv2d,time_it
def read_bin(file):
    ds = gdal.Open(file)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr

def norm_data(data):
    
    data = 10*np.log10(data)
    data[data==-np.inf] = np.nan
    data[data==np.inf] = np.nan
    data = (data - np.nanmin(data)) / (np.nanmax(data) - np.nanmin(data))

    p5 = np.nanpercentile(data, 5)
    p95 = np.nanpercentile(data, 95)
    data = np.clip(data, p5, p95)
    data = (data - p5) / (p95 - p5)

    return data   
def read_and_normalize(file_path):
    """Reads binary data and normalizes it."""
    return norm_data(read_bin(file_path)) if os.path.isfile(file_path) else None

def generate_pauli_rgb(red, green, blue, output_path):
    """Creates an RGBA image and saves it."""
    rgb_uint8 = (np.dstack((red, green, blue)) * 255).astype(np.uint8)
    alpha_channel = np.where(np.all(rgb_uint8 == 0, axis=2), 0, 255).astype(np.uint8)
    rgba_uint8 = np.dstack((rgb_uint8, alpha_channel))
    
    plt.imsave(output_path, rgba_uint8)
    print(f"Pauli RGB image saved as {output_path}")
    
    fig, ax = plt.subplots()
    plt.imshow(rgba_uint8, vmin=0, vmax=255)
    ax.axis('off')
    plt.savefig(output_path.replace(".png", "_thumb.png"), format='png', 
                bbox_inches='tight', pad_inches=0, transparent=True)

def save_rgb_geotiff(red, green, blue, georef_file, output_path):
    driver = gdal.GetDriverByName('GTiff')
    height, width = red.shape
    dataset = driver.Create(output_path, width, height, 4, gdal.GDT_Byte,options=['COMPRESS=DEFLATE','PREDICTOR=2','ZLEVEL=9', 'BIGTIFF=YES','TILED=YES'])

    # Get geotransform and projection from reference file
    ref = gdal.Open(georef_file)
    dataset.SetGeoTransform(ref.GetGeoTransform())
    dataset.SetProjection(ref.GetProjection())
    # Prepare the image bands
    red_uint8 = (red * 255).astype(np.uint8)
    green_uint8 = (green * 255).astype(np.uint8)
    blue_uint8 = (blue * 255).astype(np.uint8)
    # alpha = np.where((red_uint8 == 0) & (green_uint8 == 0) & (blue_uint8 == 0), 0, 255).astype(np.uint8)
    ref_band = ref.GetRasterBand(1).ReadAsArray()
    alpha = np.where(ref_band == 0, 0, 255).astype(np.uint8)


    dataset.GetRasterBand(1).WriteArray(red_uint8)
    dataset.GetRasterBand(2).WriteArray(green_uint8)
    dataset.GetRasterBand(3).WriteArray(blue_uint8)
    dataset.GetRasterBand(4).WriteArray(alpha)

    dataset.FlushCache()
    ref = None
    print(f"Pauli RGB GeoTIFF saved to {output_path}")
    
@time_it
def pauliRGB(infolder,tif_flag=False):
    """
    Generates a Pauli RGB image from full polarimetric datasets (`S2`,  `C4`, `T4`, `C3`, or `T3`).
    
    The function checks for the presence of specific binary data files in the given folder
    and processes the available dataset into an RGB visualization using the Pauli decomposition.
    
    Example Usage:
        ```python
        polsartools.utils.pauliRGB("/path/to/data")
        ```
    Args:
        infolder (str): Path to the folder containing the full polarimetric data files.

    Raises:
        ValueError: If the folder does not contain a valid dataset (`S2`,  `C4`, `T4`, `C3`, or `T3`).

    Output:
        - Saves `PauliRGB.png` and a thumbnail (`PauliRGB_thumb.png`) in `infolder`.
        - The images represent a false-color visualization of polarization characteristics.

    """
    
    # Check for C4 first, ensuring all required files exist
    c4_files = [os.path.join(infolder, f"C{i}.bin") for i in [11, 22, 33, 44]]
    if all(os.path.isfile(file) for file in c4_files):
        C11, C22, C33, C44 = [read_bin(file) for file in c4_files]
        C14_real = read_bin(os.path.join(infolder, "C14_real.bin"))
        red = norm_data(C11 + C44 - 2 * C14_real)
        green = norm_data(C22 + C33)
        blue = norm_data(C11 + C44 + 2 * C14_real)
        georef_file = os.path.join(infolder, "C11.bin")

    # If C4 is missing, check for C3
    elif all(os.path.isfile(os.path.join(infolder, f"C{i}.bin")) for i in [11, 22, 33]):
        C11, C22, C33 = [read_bin(os.path.join(infolder, f"C{i}.bin")) for i in [11, 22, 33]]
        blue = norm_data(0.5 * (C11 + C33 + 2 * read_bin(os.path.join(infolder, "C13_real.bin"))))
        red = norm_data(0.5 * (C11 + C33 - 2 * read_bin(os.path.join(infolder, "C13_real.bin"))))
        green = norm_data(C22)
        georef_file = os.path.join(infolder, "C11.bin")

    # Check for T3 next
    elif all(os.path.isfile(os.path.join(infolder, f"T{i}.bin")) for i in [11, 22, 33]):
        red = read_and_normalize(os.path.join(infolder, "T22.bin"))
        green = read_and_normalize(os.path.join(infolder, "T33.bin"))
        blue = read_and_normalize(os.path.join(infolder, "T11.bin"))
        georef_file = os.path.join(infolder, "T11.bin")

    # Lastly, check for S2
    elif all(os.path.isfile(os.path.join(infolder, f"s{i}.bin")) for i in [11, 12, 22]):
        S11, S12, S22 = [read_bin(os.path.join(infolder, f"s{i}.bin")) for i in [11, 12, 22]]
        blue = norm_data(np.abs((1 / np.sqrt(2)) * (S11 + S22))**2)
        red = norm_data(np.abs((1 / np.sqrt(2)) * (S11 - S22))**2)
        green = norm_data(np.abs((2 / np.sqrt(2)) * S12)**2)
        georef_file = os.path.join(infolder, "s11.bin")

    else:
        raise ValueError("Invalid S2/C3/T3/C4 folder! Provide a valid dataset.")
    
    generate_pauli_rgb(red, green, blue, os.path.join(infolder, "PauliRGB.png"))
    if tif_flag:
        save_rgb_geotiff(red, green, blue, georef_file, os.path.join(infolder, "PauliRGB.tif"))
