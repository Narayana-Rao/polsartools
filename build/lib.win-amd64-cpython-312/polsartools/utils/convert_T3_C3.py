import os
import numpy as np
from osgeo import gdal
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d
from polsartools.utils.convert_matrices import T3_C3_mat

@time_it
def convert_T3_C3(infolder, outname=None, chi_in=0, psi_in=0, window_size=1,write_flag=True,max_workers=None):

    if os.path.isfile(os.path.join(infolder,"T11.bin")):
        
        ds = gdal.Open(os.path.join(infolder,"T11.bin"))
        rows,cols = ds.RasterYSize, ds.RasterXSize 
        ds = None
        input_filepaths = [
        os.path.join(infolder,"T11.bin"),
        os.path.join(infolder,'T12_real.bin'), os.path.join(infolder,'T12_imag.bin'),  
        os.path.join(infolder,'T13_real.bin'), os.path.join(infolder,'T13_imag.bin'),
        os.path.join(infolder,"T22.bin"),
        os.path.join(infolder,'T23_real.bin'), os.path.join(infolder,'T23_imag.bin'),  
        os.path.join(infolder,"T33.bin"),
        ]
       
    # elif os.path.isfile(os.path.join(infolder,"C11.bin")):
    #     input_filepaths = [
    #     os.path.join(infolder,"C11.bin"),
    #     os.path.join(infolder,'C12_real.bin'), os.path.join(infolder,'C12_imag.bin'),  
    #     os.path.join(infolder,'C13_real.bin'), os.path.join(infolder,'C13_imag.bin'),
    #     os.path.join(infolder,"C22.bin"),
    #     os.path.join(infolder,'C23_real.bin'), os.path.join(infolder,'C23_imag.bin'),  
    #     os.path.join(infolder,"C33.bin"),
    #     ]
        os.makedirs(os.path.join(os.path.dirname(infolder), "C3"), exist_ok=True)
    else:
        raise Exception(f"Invalid T3 folder!!")

    output_filepaths = []
    if outname is None:
        output_filepaths.append(os.path.join(os.path.dirname(infolder), "C3", "C11.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(infolder), "C3", "C12_real.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(infolder), "C3", "C12_imag.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(infolder), "C3", "C13_real.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(infolder), "C3", "C13_imag.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(infolder), "C3", "C22.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(infolder), "C3", "C23_real.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(infolder), "C3", "C23_imag.bin"))
        output_filepaths.append(os.path.join(os.path.dirname(infolder), "C3", "C33.bin"))
    
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
            processing_func=process_chunk_T3C3,
            block_size=(512, 512), max_workers=max_workers, 
            num_outputs=9)
    
    file = open(os.path.join(os.path.dirname(infolder), "C3",'config.txt') ,"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()  

def process_chunk_T3C3(chunks, window_size, input_filepaths, *args):

    # additional_arg1 = args[0] if len(args) > 0 else None
    # additional_arg2 = args[1] if len(args) > 1 else None

    if 'T11' in input_filepaths[0] and 'T22' in input_filepaths[5] and 'T33' in input_filepaths[8]:
        
        t11_T1 = np.array(chunks[0])
        t12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
        t13_T1 = np.array(chunks[3])+1j*np.array(chunks[4])
        t21_T1 = np.conj(t12_T1)
        t22_T1 = np.array(chunks[5])
        t23_T1 = np.array(chunks[6])+1j*np.array(chunks[7])
        t31_T1 = np.conj(t13_T1)
        t32_T1 = np.conj(t23_T1)
        t33_T1 = np.array(chunks[8])

        T3 = np.array([[t11_T1, t12_T1, t13_T1], 
                     [t21_T1, t22_T1, t23_T1], 
                     [t31_T1, t32_T1, t33_T1]])
        T_T1 = T3_C3_mat(T3)
    else:
        raise Exception(f"Invalid T3 folder!!")

    # if 'C11' in input_filepaths[0] and 'C22' in input_filepaths[5] and 'C33' in input_filepaths[8]:
    #     C11 = np.array(chunks[0])
    #     C12 = np.array(chunks[1])+1j*np.array(chunks[2])
    #     C13 = np.array(chunks[3])+1j*np.array(chunks[4])
    #     C21 = np.conj(C12)
    #     C22 = np.array(chunks[5])
    #     C23 = np.array(chunks[6])+1j*np.array(chunks[7])
    #     C31 = np.conj(C13)
    #     C32 = np.conj(C23)
    #     C33 = np.array(chunks[8])
    #     T_T1 = np.array([[C11, C12, C13], 
    #                      [C21, C22, C23], 
    #                      [C31, C32, C33]])


    if window_size>1:
        kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)

        t11f = conv2d(T_T1[0,0,:,:],kernel)
        t12f = conv2d(np.real(T_T1[0,1,:,:]),kernel)+1j*conv2d(np.imag(T_T1[0,1,:,:]),kernel)
        t13f = conv2d(np.real(T_T1[0,2,:,:]),kernel)+1j*conv2d(np.imag(T_T1[0,2,:,:]),kernel)
        
        t21f = np.conj(t12f) 
        t22f = conv2d(T_T1[1,1,:,:],kernel)
        t23f = conv2d(np.real(T_T1[1,2,:,:]),kernel)+1j*conv2d(np.imag(T_T1[1,2,:,:]),kernel)

        t31f = np.conj(t13f) 
        t32f = np.conj(t23f) 
        t33f = conv2d(T_T1[2,2,:,:],kernel)

        T_T1 = np.array([[t11f, t12f, t13f], [t21f, t22f, t23f], [t31f, t32f, t33f]])


    return T_T1[0,0,:,:], T_T1[0,1,:,:].real, T_T1[0,1,:,:].imag, T_T1[0,2,:,:].real, T_T1[0,2,:,:].imag, T_T1[1,1,:,:], T_T1[1,2,:,:].real, T_T1[1,2,:,:].imag, T_T1[2,2,:,:]