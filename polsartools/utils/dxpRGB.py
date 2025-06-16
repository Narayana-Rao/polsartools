import os
import numpy as np
from polsartools.utils.utils import conv2d,time_it
from osgeo import gdal
import matplotlib.pyplot as plt

from .pauliRGB import generate_rgb_png, create_pgw, generate_rgb_tif, norm_data, read_bin #create_prj

@time_it
def dxpRGB(infolder, type = 1 , tif_flag = False):
    """
    Generate false-color RGB visualization from C2 dual-polarimetric SAR data.

    This function creates an RGB image from dual-pol (C2) covariance matrix inputs 
    to enhance polarimetric contrast through various combinations of C11, C22, and 
    C12_real. The output provides a quick-look visualization useful for interpreting 
    scattering mechanisms and surface features.

    Examples
    --------
    >>> # Generate type-1 false-color RGB from C2 folder
    >>> dxpRGB("/path/to/C2_data")

    >>> # Output both PNG and GeoTIFF versions
    >>> dxpRGB("/path/to/C2_data", type=3, tif_flag=True)

    Parameters
    ----------
    infolder : str
        Path to the input folder containing C2 matrix components (C11, C22, C12_real).
    type : int, {1, 2, 3, 4}, default=1
        RGB combination type:
        - 1: Blue=C11, Red=C22, Green=|C11 + C22 - 2·C12_real|
        - 2: Blue=C22, Red=C11, Green=|C11 + C22 - 2·C12_real|
        - 3: Blue=C11, Red=C22, Green=|C11 - C22|
        - 4: Blue=C22, Red=C11, Green=|C22 - C11|
    tif_flag : bool, default=False
        If True, generates a GeoTIFF (.tif) file alongside the PNG image.

    Returns
    -------
    None
        Writes RGB visualization to disk:
        - 'RGB#.png': False-color PNG image (with world file for georeferencing)
        - 'RGB#.tif': Optional GeoTIFF output if `tif_flag=True`

    """
    
    def find_file(prefix):
        for ext in [".bin", ".tif"]:
            path = os.path.join(infolder, f"{prefix}{ext}")
            if os.path.isfile(path):
                return path
        return None

    c11_path = find_file("C11")
    c22_path = find_file("C22")
    c12r_path = find_file("C12_real")

    if c11_path and c22_path:
        C11 = read_bin(c11_path)
        C22 = read_bin(c22_path)
        C12r = read_bin(c12r_path) if c12r_path else None

        if type == 1:
            blue = norm_data(C11)
            red = norm_data(C22)
            green = norm_data(np.abs(C11 + C22 - 2 * C12r))
        elif type == 2:
            blue = norm_data(C22)
            red = norm_data(C11)
            green = norm_data(np.abs(C11 + C22 - 2 * C12r))
        elif type == 3:
            blue = norm_data(C11)
            red = norm_data(C22)
            green = norm_data(np.abs(C11 - C22))
        elif type == 4:
            blue = norm_data(C22)
            red = norm_data(C11)
            green = norm_data(np.abs(C22 - C11))
        else:
            raise ValueError("Invalid type!! Valid types are 1,2,3,4")

        georef_file = c11_path
        output_path = os.path.join(infolder, f'RGB{type}.png')
        generate_rgb_png(red, green, blue,georef_file, output_path)
        create_pgw(georef_file, output_path)
        # create_prj(georef_file, output_path)
        
        print(f"RGB image saved as {output_path}")
        if tif_flag:
            output_path = os.path.join(infolder,  f'RGB{type}.png')
            generate_rgb_tif(red, green, blue, georef_file, output_path)
            print(f"RGB image saved as {output_path}")
    
    
    else:
        raise ValueError("Invalid C2 folder!!")
        
        
    