
import numpy as np
from osgeo import gdal
import os, glob
from polsartools.utils.utils import time_it, mlook_arr
from polsartools.utils.io_utils import write_T3, write_C3
gdal.UseExceptions()
def read_bin(file):
    ds = gdal.Open(file)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr

def read_a2(file):
    
    fp = open(file,mode='rb')
    fp.seek(232)
    ch = int(fp.read(4))
    # print(ch)
    
    fp.seek(236)
    nline = int(fp.read(8))
    # print(nline)
    fp.seek(248)
    npixel = int(fp.read(8))
    # print(npixel)
    
    nrec = 544 + npixel*8
    # print(nrec)
    fp.seek(720)
    data = np.frombuffer(fp.read(int(nrec * nline)), dtype='>f4')
    data = np.array(data).reshape(-1,int(nrec/4)) 
    # print(np.shape(data))
    
    data = data[:,int(544/4):int(nrec/4)] 
    slc = data[:,::2] + 1j*data[:,1::2]
    # print(np.shape(slc))
    del data
    
    return slc

def write_a2_bin(file,wdata):
    
    # ds = gdal.Open(refData)
    [cols, rows] = wdata.shape

    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    # outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    # outdata.SetProjection(ds.GetProjection())##sets same projection as input
    
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    # outdata.GetRasterBand(1).SetNoDataValue(0)##if you want these values transparent
    outdata.FlushCache() ##saves to disk!!

@time_it    
def alos2_fbd_l11(inFolder,azlks=3,rglks=2,calfac_dB=-83):
    """
    Extracts the C2 matrix elements (C11, C22, and C12) from ALOS-2 Fine Beam Dual-Pol (FBD) CEOS data 
    and saves them into respective binary files.

    Example:
    --------
    >>> alos2_fbd_l11("path_to_folder", azlks=5, rglks=3, calfac_dB=-80)
    This will extract the C2 matrix elements from the ALOS-2 Fine Beam Dual-Pol data 
    in the specified folder and save them in the 'C2' directory.
    
    Parameters:
    -----------
    inFolder : str
        The path to the folder containing the ALOS-2 Fine Beam Dual-Pol CEOS data files.

    azlks : int, optional (default=3)
        The number of azimuth looks for multi-looking.

    rglks : int, optional (default=2)
        The number of range looks for multi-looking.
        
    calfac_dB : float, optional (default=-83)
        The calibration factor in dB used to scale the raw radar data. It is applied to 
        the HH and HV polarization data before matrix computation.

    Returns:
    --------
    None
        The function does not return any value. Instead, it creates a folder named `C2` 
        (if not already present) and saves the following binary files:

        - `C11.bin`: Contains the C11 matrix elements.
        - `C22.bin`: Contains the C22 matrix elements.
        - `C12_real.bin`: Contains the real part of the C12 matrix.
        - `C12_imag.bin`: Contains the imaginary part of the C12 matrix.
        - `config.txt`: A text file containing grid dimensions and polarimetric configuration.

    Raises:
    -------
    FileNotFoundError
        If the required ALOS-2 data files (e.g., `IMG-HH` and `IMG-HV`) cannot be found in the specified folder.

    ValueError
        If the calibration factor is invalid or if the files are not in the expected format.


    """
    
    
    hh_file = list(glob.glob(os.path.join(inFolder,'IMG-HH-*-FBDR1.1__A')) + \
        glob.glob(os.path.join(inFolder, 'IMG-HH-*-FBDR1.1__D')))[0]

    hv_file = list(glob.glob(os.path.join(inFolder,'IMG-HV-*-FBDR1.1__A')) + \
        glob.glob(os.path.join(inFolder, 'IMG-HV-*-FBDR1.1__D')))[0]

    calfac_linear = np.sqrt(10 ** ((calfac_dB - 32) / 10))
    
    S11 = read_a2(hh_file).astype(np.complex64)*calfac_linear 
    S12 = read_a2(hv_file).astype(np.complex64)*calfac_linear 
    
    
    C11 = mlook_arr(np.abs(S11)**2,azlks,rglks).astype(np.float32)
    C22 = mlook_arr(np.abs(S12)**2,azlks,rglks).astype(np.float32)
    
    C12 = mlook_arr(S11*np.conjugate(S12),azlks,rglks).astype(np.complex64)
    rows,cols = C11.shape
    
    inFolder = os.path.dirname(hh_file)   
    C2Folder = os.path.join(inFolder,os.path.basename(hh_file).split('-HH-')[1].split('1.1')[0],'C2')

    del S11,S12
    
    os.makedirs(C2Folder,exist_ok=True)
    
    write_a2_bin( os.path.join(C2Folder,'C11.bin'), C11)
    print(f"Saved file {C2Folder}/C11.bin")
    del C11  
    write_a2_bin( os.path.join(C2Folder,'C22.bin'), C22)
    print(f"Saved file {C2Folder}/C22.bin")
    del C22
     
    write_a2_bin( os.path.join(C2Folder,'C12_real.bin'), np.real(C12))
    print(f"Saved file {C2Folder}/C12_real.bin")
    write_a2_bin( os.path.join(C2Folder,'C12_imag.bin'), np.imag(C12))
    print(f"Saved file {C2Folder}/C12_imag.bin")
    del C12
    
    
    file = open(C2Folder +'/config.txt',"w+")
    file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\npp1'%(rows,cols))
    file.close()  
#################################################################################################

def write_s2_bin(file,wdata):
    [cols, rows] = wdata.shape
    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_CFloat32)
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    outdata.FlushCache()
    
@time_it    
def alos2_hbq_l11(inFolder,matrix='T3', azlks=8,rglks=4,sym=True,calfac_dB=-83):

    """
    Extracts single look S2 or Multi-look T3/C3 matrix elements from ALOS-2 Quad-Pol (HBQ) CEOS data 
    and saves them into respective binary files.

    Example:
    --------
    >>> alos2_hbq_l11("path_to_folder", azlks=5, rglks=3, calfac_dB=-80)
    This will extract the T3 matrix elements from the ALOS-2 Fine Beam Dual-Pol data 
    in the specified folder and save them in the 'C2' directory.
    
    Parameters:
    -----------
    inFolder : str
        The path to the folder containing the ALOS-2 Quad-Pol (HBQ) CEOS data folder.
    
    matrix : str, optional (default='T3')
        The matrix type to extract. Can be 'S2','T3' or 'C3'.
        
    azlks : int, optional (default=8)
        The number of azimuth looks for multi-looking.

    rglks : int, optional (default=4)
        The number of range looks for multi-looking.
        
    calfac_dB : float, optional (default=-83)
        The calibration factor in dB used to scale the raw radar data. It is applied to 
        the HH and HV polarization data before matrix computation.

    Returns:
    --------
    None
        The function does not return any value. Instead, it creates a folders named 'S2` or 'C3` or 'T3` 
        (depending on the chosen matrix) and saves the corresponding binary files.

    """
    
    hh_file = list(glob.glob(os.path.join(inFolder,'IMG-HH-*-HBQR1.1__A')) + \
        glob.glob(os.path.join(inFolder, 'IMG-HH-*-HBQR1.1__D')))[0]

    hv_file = list(glob.glob(os.path.join(inFolder,'IMG-HV-*-HBQR1.1__A')) + \
        glob.glob(os.path.join(inFolder, 'IMG-HV-*-HBQR1.1__D')))[0]


    vh_file = list(glob.glob(os.path.join(inFolder,'IMG-VH-*-HBQR1.1__A')) + \
        glob.glob(os.path.join(inFolder, 'IMG-VH-*-HBQR1.1__D')))[0]

    vv_file = list(glob.glob(os.path.join(inFolder,'IMG-VV-*-HBQR1.1__A')) + \
        glob.glob(os.path.join(inFolder, 'IMG-VV-*-HBQR1.1__D')))[0]



    calfac_linear = np.sqrt(10 ** ((calfac_dB - 32) / 10))
    
    S11 = read_a2(hh_file).astype(np.complex64)*calfac_linear 
    S21 = read_a2(hv_file).astype(np.complex64)*calfac_linear 
    S12 = read_a2(vh_file).astype(np.complex64)*calfac_linear 
    S22 = read_a2(vv_file).astype(np.complex64)*calfac_linear 
    
    if matrix=='S2':
        rows,cols = S11.shape
    
        out_dir = os.path.join(inFolder,'S2')
        os.makedirs(out_dir,exist_ok=True)
        
        
        out_file = os.path.join(out_dir,'s11.bin')
        write_s2_bin(out_file,S11)
        print("Saved file "+out_file)
        
        if sym:
            out_file = os.path.join(out_dir,'s12.bin')
            write_s2_bin(out_file,(S12+S21)/2)
            print("Saved file "+out_file)
            
            out_file = os.path.join(out_dir,'s21.bin')
            write_s2_bin(out_file,(S12+S21)/2)
            print("Saved file "+out_file)
        else:
            out_file = os.path.join(out_dir,'s12.bin')
            write_s2_bin(out_file,S12)
            print("Saved file "+out_file)
            
            out_file = os.path.join(out_dir,'s21.bin')
            write_s2_bin(out_file,S21)
            print("Saved file "+out_file)            
        
        out_file = os.path.join(out_dir,'s22.bin')
        write_s2_bin(out_file,S22)
        print("Saved file "+out_file)
        
        file = open(out_dir +'/config.txt',"w+")
        file.write('Nrow\n%d\n---------\nNcol\n%d\n---------\nPolarCase\nmonostatic\n---------\nPolarType\nfull'%(rows,cols))
        file.close() 
    elif matrix=='T3':
        # 3x1 Pauli Coherency Matrix
        Kp = (1/np.sqrt(2))*np.array([S11+S22, S11-S22, S12+S21])

        del S11,S12,S21,S22

        # 3x3 Pauli Coherency Matrix elements
        T11 = mlook_arr(np.abs(Kp[0])**2,azlks,rglks).astype(np.float32)
        T22 = mlook_arr(np.abs(Kp[1])**2,azlks,rglks).astype(np.float32)
        T33 = mlook_arr(np.abs(Kp[2])**2,azlks,rglks).astype(np.float32)

        T12 = mlook_arr(Kp[0]*np.conj(Kp[1]),azlks,rglks).astype(np.complex64)
        T13 = mlook_arr(Kp[0]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)
        T23 = mlook_arr(Kp[1]*np.conj(Kp[2]),azlks,rglks).astype(np.complex64)

        del Kp
        
        
        T3Folder = os.path.join(inFolder,'T3')
        os.makedirs(T3Folder,exist_ok=True)
        
        if not os.path.isdir(T3Folder):
            print("T3 folder does not exist. \nCreating folder {}".format(T3Folder))
            os.mkdir(T3Folder)
            
        # write_T3(np.dstack([T11,T12,T13,np.conjugate(T12),T22,T23,np.conjugate(T13),np.conjugate(T23),T33]),T3Folder)
        write_T3([np.real(T11),np.real(T12),np.imag(T12),np.real(T13),np.imag(T13),
                  np.real(T22),np.real(T23),np.imag(T23),
                  np.real(T33)],T3Folder)
    elif matrix=='C3':
        # Kl- 3-D Lexicographic feature vector
        Kl = np.array([S11, np.sqrt(2)*0.5*(S12+S21), S22])
        del S11, S12, S21, S22

        # 3x3 COVARIANCE Matrix elements

        C11 = mlook_arr(np.abs(Kl[0])**2,azlks,rglks).astype(np.float32)
        C22 = mlook_arr(np.abs(Kl[1])**2,azlks,rglks).astype(np.float32)
        C33 = mlook_arr(np.abs(Kl[2])**2,azlks,rglks).astype(np.float32)

        C12 = mlook_arr(Kl[0]*np.conj(Kl[1]),azlks,rglks).astype(np.complex64)
        C13 = mlook_arr(Kl[0]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)
        C23 = mlook_arr(Kl[1]*np.conj(Kl[2]),azlks,rglks).astype(np.complex64)


        C3Folder = os.path.join(inFolder,'C3')
        os.makedirs(C3Folder,exist_ok=True)
        
        
        if not os.path.isdir(C3Folder):
            print("C3 folder does not exist. \nCreating folder {}".format(C3Folder))
            os.mkdir(C3Folder)
        
        # write_C3(np.dstack([C11,C12,C13,np.conjugate(C12),C22,C23,np.conjugate(C13),np.conjugate(C23),C33]),C3Folder)
        write_C3([np.real(C11),np.real(C12),np.imag(C12),np.real(C13),np.imag(C13),
                  np.real(C22),np.real(C23),np.imag(C23),
                  np.real(C33)],C3Folder)
        
        
    else:
        raise ValueError('Invalid matrix type. Valid types are "S2", "T3" and "C3"')
        