import os
import numpy as np
from osgeo import gdal
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

@time_it
def pauliRGB(infolder, outname=None, chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    if os.path.isfile(os.path.join(infolder,"T11.bin")) and os.path.isfile(os.path.join(infolder,"T22.bin")) and os.path.isfile(os.path.join(infolder,"T33.bin")):
        blue_ = norm_data(read_bin(os.path.join(infolder,"T11.bin")))
        red_ = norm_data(read_bin(os.path.join(infolder,"T22.bin")))
        green_ = norm_data(read_bin(os.path.join(infolder,"T33.bin")))

        rgb_uint8 = (np.dstack((red_,green_,blue_)) * 255) .astype(np.uint8)
        del blue_, red_, green_
        alpha_channel = np.where(np.all(rgb_uint8 == 0, axis=2), 0, 255).astype(np.uint8)
        rgba_uint8 = np.dstack((rgb_uint8, alpha_channel))
        del rgb_uint8, alpha_channel

        plt.imsave(os.path.join(infolder,"PauliRGB.png"),rgba_uint8)       
        print(f"Pauli RGB image saved as {os.path.join(infolder, 'PauliRGB.png')}")
        
        fig,ax = plt.subplots()
        plt.imshow(rgba_uint8,vmin=0,vmax=255)
        ax.axis('off')
        plt.savefig(os.path.join(infolder,"PauliRGB_thumb.png"), format='png', bbox_inches='tight', pad_inches=0,transparent=True)
    elif os.path.isfile(os.path.join(infolder,"C11.bin")) and os.path.isfile(os.path.join(infolder,"C22.bin")) and os.path.isfile(os.path.join(infolder,"C33.bin")):
        blue_ = norm_data(0.5*(read_bin(os.path.join(infolder,'C11.bin'))+read_bin(os.path.join(infolder,'C33.bin'))+2*read_bin(os.path.join(infolder,'C13_real.bin'))))
        red_ = norm_data(0.5*(read_bin(os.path.join(infolder,'C11.bin'))+read_bin(os.path.join(infolder,'C33.bin'))-2*read_bin(os.path.join(infolder,'C13_real.bin'))))
        green_ = norm_data(read_bin(os.path.join(infolder,"C22.bin")))
        
        rgb_uint8 = (np.dstack((red_,green_,blue_)) * 255) .astype(np.uint8)
        del blue_, red_, green_
        alpha_channel = np.where(np.all(rgb_uint8 == 0, axis=2), 0, 255).astype(np.uint8)
        rgba_uint8 = np.dstack((rgb_uint8, alpha_channel))
        del rgb_uint8, alpha_channel
        
        plt.imsave(os.path.join(infolder,"PauliRGB.png"),rgba_uint8)
        print(f"Pauli RGB image saved as {os.path.join(infolder, 'PauliRGB.png')}")
        fig,ax = plt.subplots()
        plt.imshow(rgba_uint8,vmin=0,vmax=255)
        ax.axis('off')
        plt.savefig(os.path.join(infolder,"PauliRGB_thumb.png"), format='png', bbox_inches='tight', pad_inches=0,transparent=True)
    
    elif os.path.isfile(os.path.join(infolder,"s11.bin")) and os.path.isfile(os.path.join(infolder,"s22.bin")) and os.path.isfile(os.path.join(infolder,"s12.bin")):
        b_data =  (1/np.sqrt(2))*(read_bin(os.path.join(infolder,"s11.bin")) + read_bin(os.path.join(infolder,"s22.bin")))
        blue_ = norm_data(np.abs(b_data)**2)
        
        b_data =  (1/np.sqrt(2))*(read_bin(os.path.join(infolder,"s11.bin")) - read_bin(os.path.join(infolder,"s22.bin")))
        red_ = norm_data(np.abs(b_data)**2)
        
        b_data =  (2/np.sqrt(2))*read_bin(os.path.join(infolder,"s12.bin")) 
        green_ = norm_data(np.abs(b_data)**2)
        del b_data

        rgb_uint8 = (np.dstack((red_,green_,blue_)) * 255) .astype(np.uint8)
        del blue_, red_, green_
        alpha_channel = np.where(np.all(rgb_uint8 == 0, axis=2), 0, 255).astype(np.uint8)
        rgba_uint8 = np.dstack((rgb_uint8, alpha_channel))
        del rgb_uint8, alpha_channel

        plt.imsave(os.path.join(infolder,"PauliRGB.png"),rgba_uint8)       
        print(f"Pauli RGB image saved as {os.path.join(infolder, 'PauliRGB.png')}")
        
        fig,ax = plt.subplots()
        plt.imshow(rgba_uint8,vmin=0,vmax=255)
        ax.axis('off')
        plt.savefig(os.path.join(infolder,"PauliRGB_thumb.png"), format='png', bbox_inches='tight', pad_inches=0,transparent=True)

    else:
        raise ValueError("Invalid S2/C3/T3 folder!! \n Please provide a valid S2/C3/T3 folder")
        
        
  