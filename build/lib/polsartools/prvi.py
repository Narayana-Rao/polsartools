from osgeo import gdal
import numpy as np
import os 
import warnings
warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d


def grvi(T3_folder,window_size=1,write_flag=None):

    T11 = read_bin(T3_folder+"/T11.bin")
    T22 = read_bin(T3_folder+"/T22.bin")
    T33 = read_bin(T3_folder+"/T33.bin")

    T12_i = read_bin(T3_folder+'/T12_imag.bin')
    T12_r = read_bin(T3_folder+'/T12_real.bin')
    T13_i = read_bin(T3_folder+'/T13_imag.bin')
    T13_r = read_bin(T3_folder+'/T13_real.bin')
    T23_i = read_bin(T3_folder+'/T23_imag.bin')
    T23_r = read_bin(T3_folder+'/T23_real.bin')
    T12 = T12_r + 1j*T12_i
    T13 = T13_r + 1j*T13_i
    T23 = T23_r + 1j*T23_i

    T3 = np.dstack((T11,T12,T13,np.conj(T12),T22,T23,np.conj(T13),np.conj(T23),T33))
   
    DOP = np.zeros([np.size(T3,0),np.size(T3,1)])
    prvi = np.zeros([np.size(T3,0),np.size(T3,1)])
    "Special Unitary Matrix"
    D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
    
    for i in range(np.size(T3,0)):

        self.pBar.emit((i/np.size(T3,0))*100)
        for j in range(np.size(T3,1)):
            det = np.abs(np.linalg.det(T3[i,j,:].reshape((3, 3))))
            trace = np.abs(np.trace(T3[i,j,:].reshape((3, 3))))
            if trace==0:
                DOP[i][j]=0
            # elif((27*det/trace**3)>1.0):
            #     DOP[i][j]=0
            else:
                DOP[i][j] = np.sqrt(1-((27*det)/trace**3))
            
            tempT3 =  np.reshape(T3[i,j,:],(3,3))
            C3 = np.matmul(np.matmul((D.T),tempT3),D).flatten()
            
            prvi[i][j] =(1- DOP[i][j])*C3[4]*0.5 # (1-dop)*vh


    if write_flag:
        infile = T3_folder+'/T11.bin'
        """Write files to disk"""
        if os.path.exists(T3_folder+'/T11.bin'):
            infile = iFolder+'/T11.bin'
        elif os.path.exists(T3_folder+'/C11.bin'):
            infile = T3_folder+'/C11.bin'

        ofilegrvi = T3_folder+'/PRVI.bin'
        write_bin(ofilegrvi,prvi,infile)
                    

    return prvi
