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
def pauliRGB(infolder, outname=None, chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):
    if os.path.isfile(os.path.join(infolder,"T11.bin")) and os.path.isfile(os.path.join(infolder,"T22.bin")) and os.path.isfile(os.path.join(infolder,"T33.bin")):
        blue_ = norm_data(read_bin(os.path.join(infolder,"T11.bin")))
        red_ = norm_data(read_bin(os.path.join(infolder,"T22.bin")))
        green_ = norm_data(read_bin(os.path.join(infolder,"T33.bin")))

        rgb_uint8 = (np.dstack((red_,green_,blue_)) * 255) .astype(np.uint8)

        plt.imsave(os.path.join(infolder,"PauliRGB.png"),rgb_uint8)
        print(f"Pauli RGB image saved as {os.path.join(infolder,"PauliRGB.png")}PauliRGB.png")
        fig,ax = plt.subplots()
        plt.imshow(rgb_uint8,vmin=0,vmax=255)
        ax.axis('off')
        plt.savefig(os.path.join(infolder,"PauliRGB_thumb.png"), format='png', bbox_inches='tight', pad_inches=0)
    elif os.path.isfile(os.path.join(infolder,"C11.bin")) and os.path.isfile(os.path.join(infolder,"C22.bin")) and os.path.isfile(os.path.join(infolder,"C33.bin")):
        blue_ = norm_data(0.5*(read_bin(os.path.join(infolder,'C11.bin'))+read_bin(os.path.join(infolder,'C33.bin'))+2*read_bin(os.path.join(infolder,'C13_real.bin'))))
        red_ = norm_data(0.5*(read_bin(os.path.join(infolder,'C11.bin'))+read_bin(os.path.join(infolder,'C33.bin'))-2*read_bin(os.path.join(infolder,'C13_real.bin'))))
        green_ = norm_data(read_bin(os.path.join(infolder,"C22.bin")))

        rgb_uint8 = (np.dstack((red_,green_,blue_)) * 255) .astype(np.uint8)

        plt.imsave(os.path.join(infolder,"PauliRGB.png"),rgb_uint8)
        print(f"Pauli RGB image saved as {os.path.join(infolder,"PauliRGB.png")}PauliRGB.png")
        fig,ax = plt.subplots()
        plt.imshow(rgb_uint8,vmin=0,vmax=255)
        ax.axis('off')
        plt.savefig(os.path.join(infolder,"PauliRGB_thumb.png"), format='png', bbox_inches='tight', pad_inches=0)
    else:
        raise ValueError("Invalid C3/T3 folder!!")
        
        
        
        
        
        
        
# @time_it
# def pauliRGB(infolder, outname=None, chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):

#     if os.path.isfile(os.path.join(infolder,"T11.bin")):
#         input_filepaths = [
#         os.path.join(infolder,"T11.bin"),
#         os.path.join(infolder,'T12_real.bin'), os.path.join(infolder,'T12_imag.bin'),  
#         os.path.join(infolder,'T13_real.bin'), os.path.join(infolder,'T13_imag.bin'),
#         os.path.join(infolder,"T22.bin"),
#         os.path.join(infolder,'T23_real.bin'), os.path.join(infolder,'T23_imag.bin'),  
#         os.path.join(infolder,"T33.bin"),
#         ]
#     elif os.path.isfile(os.path.join(infolder,"C11.bin")):
#         input_filepaths = [
#         os.path.join(infolder,"C11.bin"),
#         os.path.join(infolder,'C12_real.bin'), os.path.join(infolder,'C12_imag.bin'),  
#         os.path.join(infolder,'C13_real.bin'), os.path.join(infolder,'C13_imag.bin'),
#         os.path.join(infolder,"C22.bin"),
#         os.path.join(infolder,'C23_real.bin'), os.path.join(infolder,'C23_imag.bin'),  
#         os.path.join(infolder,"C33.bin"),
#         ]

#     else:
#         print(f"Invalid C3 or T3 folder!!")

#     output_filepaths = []
#     if outname is None:
#         output_filepaths.append(os.path.join(infolder, "red.bin"))
#         output_filepaths.append(os.path.join(infolder, "green.bin"))
#         output_filepaths.append(os.path.join(infolder, "blue.bin"))
    
#     process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
#             processing_func=process_chunk_pauliRGB,
#             block_size=(512, 512), max_workers=max_workers, 
#             num_outputs=3)

    
#     rgb_uint8 = np.dstack((read_bin(output_filepaths[0]),
#                            read_bin(output_filepaths[1]),
#                            read_bin(output_filepaths[2]))).astype(np.uint8)
    
#     plt.imsave(os.path.join(infolder,'PauliRGB.png'),rgb_uint8) 
#     # # plt.close(0)
#     fig,ax = plt.subplots()
#     plt.imshow(rgb_uint8,vmin=0,vmax=255)
#     ax.axis('off')
#     plt.savefig(os.path.join(infolder,'PauliRGB_thumb.png'), format='png', bbox_inches='tight', pad_inches=0)
#     plt.close()
    
#     os.remove(output_filepaths[0])
#     os.remove(output_filepaths[1])
#     os.remove(output_filepaths[2])
# def process_chunk_pauliRGB(chunks, window_size, input_filepaths, *args):

#     # additional_arg1 = args[0] if len(args) > 0 else None
#     # additional_arg2 = args[1] if len(args) > 1 else None

#     if 'T11' in input_filepaths[0] and 'T22' in input_filepaths[5] and 'T33' in input_filepaths[8]:
#         t11_T1 = np.array(chunks[0])
#         t12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
#         t13_T1 = np.array(chunks[3])+1j*np.array(chunks[4])
#         t21_T1 = np.conj(t12_T1)
#         t22_T1 = np.array(chunks[5])
#         t23_T1 = np.array(chunks[6])+1j*np.array(chunks[7])
#         t31_T1 = np.conj(t13_T1)
#         t32_T1 = np.conj(t23_T1)
#         t33_T1 = np.array(chunks[8])

#         T_T1 = np.array([[t11_T1, t12_T1, t13_T1], 
#                      [t21_T1, t22_T1, t23_T1], 
#                      [t31_T1, t32_T1, t33_T1]])
#         print(T_T1.shape)

#     if 'C11' in input_filepaths[0] and 'C22' in input_filepaths[5] and 'C33' in input_filepaths[8]:
#         C11 = np.array(chunks[0])
#         C12 = np.array(chunks[1])+1j*np.array(chunks[2])
#         C13 = np.array(chunks[3])+1j*np.array(chunks[4])
#         C21 = np.conj(C12)
#         C22 = np.array(chunks[5])
#         C23 = np.array(chunks[6])+1j*np.array(chunks[7])
#         C31 = np.conj(C13)
#         C32 = np.conj(C23)
#         C33 = np.array(chunks[8])
#         C3 = np.array([[C11, C12, C13], 
#                          [C21, C22, C23], 
#                          [C31, C32, C33]])

#         T_T1 = C3_T3_mat(C3)


#     if window_size>1:
#         kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

#         t11f = conv2d(T_T1[0,0,:,:],kernel)
#         t12f = conv2d(np.real(T_T1[0,1,:,:]),kernel)+1j*conv2d(np.imag(T_T1[0,1,:,:]),kernel)
#         t13f = conv2d(np.real(T_T1[0,2,:,:]),kernel)+1j*conv2d(np.imag(T_T1[0,2,:,:]),kernel)
        
#         t21f = np.conj(t12f) 
#         t22f = conv2d(T_T1[1,1,:,:],kernel)
#         t23f = conv2d(np.real(T_T1[1,2,:,:]),kernel)+1j*conv2d(np.imag(T_T1[1,2,:,:]),kernel)

#         t31f = np.conj(t13f) 
#         t32f = np.conj(t23f) 
#         t33f = conv2d(T_T1[2,2,:,:],kernel)

#         T_T1 = np.array([[t11f, t12f, t13f], [t21f, t22f, t23f], [t31f, t32f, t33f]])

#     blue_ = norm_data(np.real(T_T1[0,0,:,:]))
#     red_ = norm_data(np.real(T_T1[1,1,:,:]))
#     green_ = norm_data(np.real(T_T1[2,2,:,:]))

#     rgb_uint8 = (np.dstack((red_,green_,blue_)) * 255) .astype(np.uint8)

#     return rgb_uint8[:,:,0],rgb_uint8[:,:,1],rgb_uint8[:,:,2]