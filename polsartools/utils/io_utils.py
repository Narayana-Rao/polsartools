import numpy as np
from osgeo import gdal 
from skimage.util.shape import view_as_blocks
import os
def read_bin(file):
    ds = gdal.Open(file,gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr

def mlook(data,az,rg):
    temp = data[0:data.shape[0]-data.shape[0]%az,0:data.shape[1]-data.shape[1]%rg]
    blocks = view_as_blocks(temp, block_shape=(az, rg))
    flatten = blocks.reshape(blocks.shape[0], blocks.shape[1], -1)
    mean = np.nanmean(flatten, axis=2)
    return mean

def write_s2_bin(file,wdata):
    [cols, rows] = wdata.shape
    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_CFloat32)
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    outdata.FlushCache()
    
def write_T3(T3_list,folder):
    
    out_file = folder +'/T11.bin'
    write_bin(out_file,T3_list[0])
    print("Saved file "+out_file)
    
    out_file = folder +'/T12_real.bin'
    write_bin(out_file,T3_list[1])
    print("Saved file "+out_file)
    out_file = folder +'/T12_imag.bin'
    write_bin(out_file,T3_list[2])
    print("Saved file "+out_file)
    
    out_file = folder +'/T13_real.bin'
    write_bin(out_file,T3_list[3])
    print("Saved file "+out_file)
    out_file = folder +'/T13_imag.bin'
    write_bin(out_file,T3_list[4])
    print("Saved file "+out_file)
    
    out_file = folder +'/T22.bin'
    write_bin(out_file,T3_list[5])
    print("Saved file "+out_file)
    
    out_file = folder +'/T23_real.bin'
    write_bin(out_file,T3_list[6])
    print("Saved file "+out_file)
    out_file = folder +'/T23_imag.bin'
    write_bin(out_file,T3_list[7])
    print("Saved file "+out_file)
    out_file = folder +'/T33.bin'
    write_bin(out_file,T3_list[8])
    print("Saved file "+out_file)
    
    rows, cols = np.shape(T3_list[0])
    file = folder +'/config.txt'
    file = open(file,"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()   

def write_C3(C3_list,folder):
    
    out_file = folder +'/C11.bin'
    write_bin(out_file,C3_list[0])
    print("Saved file "+out_file)
    out_file = folder +'/C12_real.bin'
    write_bin(out_file,C3_list[1])
    print("Saved file "+out_file)
    out_file = folder +'/C12_imag.bin'
    write_bin(out_file,C3_list[2])
    print("Saved file "+out_file)
    
    out_file = folder +'/C13_real.bin'
    write_bin(out_file,C3_list[3])
    print("Saved file "+out_file)
    out_file = folder +'/C13_imag.bin'
    write_bin(out_file,C3_list[4])
    print("Saved file "+out_file)
    
    
    out_file = folder +'/C22.bin'
    write_bin(out_file,C3_list[5])
    print("Saved file "+out_file)
    out_file = folder +'/C23_real.bin'
    write_bin(out_file,C3_list[6])
    print("Saved file "+out_file)
    out_file = folder +'/C23_imag.bin'
    write_bin(out_file,C3_list[7])
    print("Saved file "+out_file)
    out_file = folder +'/C33.bin'
    write_bin(out_file,C3_list[8])
    print("Saved file "+out_file)
    
    rows, cols = np.shape(C3_list[0])
    file = folder +'/config.txt'
    file = open(file,"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()

    
def write_C4(C4_list,folder):
    
    out_file = folder +'/C11.bin'
    write_bin(out_file,C4_list[0])
    print("Saved file "+out_file)
    out_file = folder +'/C12_real.bin'
    write_bin(out_file,C4_list[1])
    print("Saved file "+out_file)
    out_file = folder +'/C12_imag.bin'
    write_bin(out_file,C4_list[2])
    print("Saved file "+out_file)
    

    
    out_file = folder +'/C13_real.bin'
    write_bin(out_file,C4_list[3])
    print("Saved file "+out_file)
    out_file = folder +'/C13_imag.bin'
    write_bin(out_file,C4_list[4])
    print("Saved file "+out_file)
    
    out_file = folder +'/C14_real.bin'
    write_bin(out_file,C4_list[5])
    print("Saved file "+out_file)
    out_file = folder +'/C14_imag.bin'
    write_bin(out_file,C4_list[6])
    print("Saved file "+out_file)

    out_file = folder +'/C22.bin'
    write_bin(out_file,C4_list[7])
    print("Saved file "+out_file)
    out_file = folder +'/C23_real.bin'
    write_bin(out_file,C4_list[8])
    print("Saved file "+out_file)
    out_file = folder +'/C23_imag.bin'
    write_bin(out_file,C4_list[9])
    print("Saved file "+out_file)
    
    out_file = folder +'/C24_real.bin'
    write_bin(out_file,C4_list[10])
    print("Saved file "+out_file)
    out_file = folder +'/C24_imag.bin'
    write_bin(out_file,C4_list[11])
    print("Saved file "+out_file)
    
    
    
    out_file = folder +'/C33.bin'
    write_bin(out_file,C4_list[12])
    print("Saved file "+out_file)
    
    
    out_file = folder +'/C34_real.bin'
    write_bin(out_file,C4_list[13])
    print("Saved file "+out_file)
    out_file = folder +'/C34_imag.bin'
    write_bin(out_file,C4_list[14])
    print("Saved file "+out_file)    
    
    
    out_file = folder +'/C44.bin'
    write_bin(out_file,C4_list[15])
    print("Saved file "+out_file)  
    
    
    
    rows, cols = np.shape(C4_list[0])
    file = folder +'/config.txt'
    file = open(file,"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()

def write_T4(T4_list,folder):
    
    out_file = folder +'/T11.bin'
    write_bin(out_file,T4_list[0])
    print("Saved file "+out_file)
    
    out_file = folder +'/T12_real.bin'
    write_bin(out_file,T4_list[1])
    print("Saved file "+out_file)
    out_file = folder +'/T12_imag.bin'
    write_bin(out_file,T4_list[2])
    print("Saved file "+out_file)
    
    out_file = folder +'/T13_real.bin'
    write_bin(out_file,T4_list[3])
    print("Saved file "+out_file)
    out_file = folder +'/T13_imag.bin'
    write_bin(out_file,T4_list[4])
    print("Saved file "+out_file)
    
    out_file = folder +'/T14_real.bin'
    write_bin(out_file,T4_list[5])
    print("Saved file "+out_file)
    out_file = folder +'/T14_imag.bin'
    write_bin(out_file,T4_list[6])
    print("Saved file "+out_file)

    out_file = folder +'/T22.bin'
    write_bin(out_file,T4_list[7])
    print("Saved file "+out_file)
    out_file = folder +'/T23_real.bin'
    write_bin(out_file,T4_list[8])
    print("Saved file "+out_file)
    out_file = folder +'/T23_imag.bin'
    write_bin(out_file,T4_list[9])
    print("Saved file "+out_file)
    
    out_file = folder +'/T24_real.bin'
    write_bin(out_file,T4_list[10])
    print("Saved file "+out_file)
    out_file = folder +'/T24_imag.bin'
    write_bin(out_file,T4_list[11])
    print("Saved file "+out_file)
    
    
    
    out_file = folder +'/T33.bin'
    write_bin(out_file,T4_list[12])
    print("Saved file "+out_file)    
    
    
    out_file = folder +'/T34_real.bin'
    write_bin(out_file,T4_list[13])
    print("Saved file "+out_file)    
    out_file = folder +'/T34_imag.bin'
    write_bin(out_file,T4_list[14])
    print("Saved file "+out_file)    
    out_file = folder +'/T44.bin'
    write_bin(out_file,T4_list[15])
    print("Saved file "+out_file)   
    
    
    rows, cols = np.shape(T4_list[0])
    file = folder +'/config.txt'
    file = open(file,"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
    file.close()
    
def write_bin_s2(file,wdata,refData):
    
    # ds = gdal.Open(refData)
    [cols, rows] = wdata.shape

    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    # outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    # outdata.SetProjection(ds.GetProjection())##sets same projection as input
    
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    # outdata.GetRasterBand(1).SetNoDataValue(np.NaN)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!!   
    

def write_bin(file,wdata):
    
    # ds = gdal.Open(refData)
    [cols, rows] = wdata.shape

    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    # outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    # outdata.SetProjection(ds.GetProjection())##sets same projection as input
    
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    # outdata.GetRasterBand(1).SetNoDataValue(np.NaN)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!! 



def write_s2_ct(matrix_type, matrixFolder, K, azlks, rglks):
    """
    Processes and saves C3, T3, C4, or T4 matrices.
    
    Args:
        matrix_type (str): One of "C3", "T3", "C4", "T4".
        inFile (str): Input file path.
        K (np.ndarray): Covariance matrix elements.
        azlks (int): Azimuth looks.
        rglks (int): Range looks.

    Saves:
        Binary files for matrix elements in the appropriate directory.
        Configuration file specifying matrix dimensions.
    """
    
    if matrix_type not in ["C3", "T3", "C4", "T4"]:
        raise ValueError("Invalid matrix type. Choose 'C3', 'T3', 'C4', or 'T4'.")

    prefix = matrix_type[0]  # Extract first character ('C' or 'T')
    

    # Compute matrix elements
    matrices = {
        f'{prefix}11': mlook(np.abs(K[0])**2, azlks, rglks).astype(np.float32),
        f'{prefix}22': mlook(np.abs(K[1])**2, azlks, rglks).astype(np.float32),
        f'{prefix}33': mlook(np.abs(K[2])**2, azlks, rglks).astype(np.float32),
        f'{prefix}12_real': mlook(np.real(K[0] * np.conj(K[1])), azlks, rglks).astype(np.float32),
        f'{prefix}12_imag': mlook(np.imag(K[0] * np.conj(K[1])), azlks, rglks).astype(np.float32),
        f'{prefix}13_real': mlook(np.real(K[0] * np.conj(K[2])), azlks, rglks).astype(np.float32),
        f'{prefix}13_imag': mlook(np.imag(K[0] * np.conj(K[2])), azlks, rglks).astype(np.float32),
        f'{prefix}23_real': mlook(np.real(K[1] * np.conj(K[2])), azlks, rglks).astype(np.float32),
        f'{prefix}23_imag': mlook(np.imag(K[1] * np.conj(K[2])), azlks, rglks).astype(np.float32),
    }

    # Extend for C4/T4 processing if K has 4 elements
    if len(K) == 4:
        matrices.update({
            f'{prefix}44': mlook(np.abs(K[3])**2, azlks, rglks).astype(np.float32),
            f'{prefix}14_real': mlook(np.real(K[0] * np.conj(K[3])), azlks, rglks).astype(np.float32),
            f'{prefix}14_imag': mlook(np.imag(K[0] * np.conj(K[3])), azlks, rglks).astype(np.float32),
            f'{prefix}24_real': mlook(np.real(K[1] * np.conj(K[3])), azlks, rglks).astype(np.float32),
            f'{prefix}24_imag': mlook(np.imag(K[1] * np.conj(K[3])), azlks, rglks).astype(np.float32),
            f'{prefix}34_real': mlook(np.real(K[2] * np.conj(K[3])), azlks, rglks).astype(np.float32),
            f'{prefix}34_imag': mlook(np.imag(K[2] * np.conj(K[3])), azlks, rglks).astype(np.float32),
        })

    # Save matrix components
    for key, matrix in matrices.items():
        out_file = os.path.join(matrixFolder, f"{key}.bin")
        write_bin(out_file, matrix)
        print(f"Saved file {out_file}")

    # Write configuration file
    rows, cols = matrices[f'{prefix}11'].shape
    config_file_path = os.path.join(matrixFolder, 'config.txt')
    with open(config_file_path, "w") as file:
        file.write(f'Nrow\n{rows}\n---------\nNcol\n{cols}\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull')
