import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d
from polsartools.utils.convert_matrices import C3_T3_mat
from osgeo import gdal
import matplotlib.pyplot as plt

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

@time_it
def dxpRGB(infolder, type = 1,outname=None, chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    if os.path.isfile(os.path.join(infolder,"C11.bin")) and os.path.isfile(os.path.join(infolder,"C22.bin")):
        
        # blue_ = norm_data(C11_m)
        # red_ = norm_data(C22_m)
        # green_ = norm_data(np.abs(C11_m+C22_m-2*C12_m.real))
        if type == 1:
            blue_ = norm_data(read_bin(os.path.join(infolder,"C11.bin")))
            red_ = norm_data(read_bin(os.path.join(infolder,"C22.bin")))
            g_data = read_bin(os.path.join(infolder,"C11.bin"))+read_bin(os.path.join(infolder,"C22.bin"))-2*read_bin(os.path.join(infolder,"C12_real.bin"))
            green_ = norm_data(np.abs(g_data))
            
        elif type == 2:
            blue_ = norm_data(read_bin(os.path.join(infolder,"C22.bin")))
            red_ = norm_data(read_bin(os.path.join(infolder,"C11.bin")))
            g_data = read_bin(os.path.join(infolder,"C11.bin"))+read_bin(os.path.join(infolder,"C22.bin"))-2*read_bin(os.path.join(infolder,"C12_real.bin"))
            green_ = norm_data(np.abs(g_data))
        elif type == 3:
            blue_ = norm_data(read_bin(os.path.join(infolder,"C11.bin")))
            red_ = norm_data(read_bin(os.path.join(infolder,"C22.bin")))
            g_data = read_bin(os.path.join(infolder,"C11.bin"))-read_bin(os.path.join(infolder,"C22.bin"))
            green_ = norm_data(np.abs(g_data))
        elif type == 4:
            blue_ = norm_data(read_bin(os.path.join(infolder,"C22.bin")))
            red_ = norm_data(read_bin(os.path.join(infolder,"C11.bin")))
            g_data = read_bin(os.path.join(infolder,"C22.bin"))-read_bin(os.path.join(infolder,"C11.bin"))
            green_ = norm_data(np.abs(g_data))
        else:
            raise ValueError("Invalid type!! Valid types are 1,2,3,4")


        rgb_uint8 = (np.dstack((red_,green_,blue_)) * 255) .astype(np.uint8)
        alpha_channel = np.where(np.all(rgb_uint8 == 0, axis=2), 0, 255).astype(np.uint8)
        rgba_uint8 = np.dstack((rgb_uint8, alpha_channel))
        

        plt.imsave(os.path.join(infolder,f"RGB{type}.png"),rgba_uint8)       
        print(f"RGB image saved as {os.path.join(infolder, f'RGB{type}.png')}")
        
        fig,ax = plt.subplots()
        plt.imshow(rgba_uint8,vmin=0,vmax=255)
        ax.axis('off')
        plt.savefig(os.path.join(infolder,f"RGB{type}_thumb.png"), format='png', bbox_inches='tight', pad_inches=0,transparent=True)
    else:
        raise ValueError("Invalid C2 folder!!")
        
        
    