from osgeo import gdal
import numpy as np
import os 
# import warnings
from shutil import copyfile
# warnings.filterwarnings('ignore')

from .basic_func import read_bin, write_bin, conv2d

def T3_C3(T3_stack):
    nrows = np.size(T3_stack,0)
    ncols = np.size(T3_stack,1)
    
    D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
    TT1 = np.dstack((T3_stack[:,:,0].flatten(), T3_stack[:,:,1].flatten(), T3_stack[:,:,2].flatten(),
                     T3_stack[:,:,3].flatten(), T3_stack[:,:,4].flatten(), T3_stack[:,:,5].flatten(),
                     T3_stack[:,:,6].flatten(), T3_stack[:,:,7].flatten(), T3_stack[:,:,8].flatten()))
    TT1 = TT1[0,:,:]
    TT1 = np.reshape(TT1,(TT1.shape[0],3,3))
    TT1 = np.matmul(np.matmul((D.T),TT1),D)
    TT1 = np.reshape(TT1,(TT1.shape[0],9))

    return np.reshape(TT1,(nrows,ncols,9)) 

def C3_T3(C3_stack):
    nrows = np.size(C3_stack,0)
    ncols = np.size(C3_stack,1)
    
    D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
    CC1 = np.dstack((C3_stack[:,:,0].flatten(), C3_stack[:,:,1].flatten(), C3_stack[:,:,2].flatten(),
                     C3_stack[:,:,3].flatten(), C3_stack[:,:,4].flatten(), C3_stack[:,:,5].flatten(),
                     C3_stack[:,:,6].flatten(), C3_stack[:,:,7].flatten(), C3_stack[:,:,8].flatten()))
    CC1 = CC1[0,:,:]
    CC1 = np.reshape(CC1,(CC1.shape[0],3,3))
    CC1 = np.matmul(np.matmul((D),CC1),D.T)
    CC1 = np.reshape(CC1,(CC1.shape[0],9))

    return np.reshape(CC1,(nrows,ncols,9)) 

def C3_compact_C2(C3,psi,chi):


	C11 = C3[:,:,0]
	C12 = C3[:,:,1]
	C13 = C3[:,:,2]

	C21 = C3[:,:,3]
	C22 = C3[:,:,4]
	C23 = C3[:,:,5]

	C31 = C3[:,:,6]
	C32 = C3[:,:,7]
	C33 = C3[:,:,8]


	#%%
	psi = psi*np.pi/180
	chi = chi*np.pi/180
	 
	CP11 = 0.5*((1+np.cos(2*psi)*np.cos(2*chi))*C11+
	             0.5*(1-np.cos(2*psi)*np.cos(2*chi))*C22+
	             (1/np.sqrt(2))*(np.sin(2*psi)*np.cos(2*chi))*(C12+np.conj(C12))+
	                (1j/np.sqrt(2))*np.sin(2*chi)*(C12-np.conj(C12)) 
	           )

	CP12 = 0.5*((1/np.sqrt(2))*(1+np.cos(2*psi)*np.cos(2*chi))*C3[:,:,1]+
	            (1/np.sqrt(2))*(1-np.cos(2*psi)*np.cos(2*chi))*C3[:,:,5]+
	            (np.sin(2*psi)*np.cos(2*chi))*(C3[:,:,2]+0.5*C3[:,:,4])+
	            1j*np.sin(2*chi)*(C3[:,:,2]-0.5*C3[:,:,4])
	            )

	CP22 = 0.5*(0.5*(1+np.cos(2*psi)*np.cos(2*chi))*C3[:,:,4]+
	                (1-np.cos(2*psi)*np.cos(2*chi))*C3[:,:,8]+
	              (1/np.sqrt(2))*(np.sin(2*psi)*np.cos(2*chi))*(C3[:,:,5]+np.conj(C3[:,:,5]))+
	                (1j/np.sqrt(2))*np.sin(2*chi)*(C3[:,:,5]-np.conj(C3[:,:,5])) 
	            )

	return np.dstack((CP11, CP12, np.conjugate(CP12),CP22))


def C3_c21(C3_stack):
	return np.dstack((np.real(C3_stack[:,:,0]),np.real((C3_stack[:,:,1])/np.sqrt(2))+1j*np.imag((C3_stack[:,:,1])/np.sqrt(2)), np.real((C3_stack[:,:,1])/np.sqrt(2))-1j*np.imag((C3_stack[:,:,1])/np.sqrt(2)), np.real((C3_stack[:,:,4])/(2))))

def C3_c22(C3_stack):
	return np.dstack((np.real(C3_stack[:,:,8]),np.real((C3_stack[:,:,5])/np.sqrt(2))+1j*np.imag((C3_stack[:,:,5])/np.sqrt(2)), np.real((C3_stack[:,:,5])/np.sqrt(2))-1j*np.imag((C3_stack[:,:,5])/np.sqrt(2)), np.real((C3_stack[:,:,4])/(2))))
	
def C3_c23(C3_stack):
	return np.dstack((np.real(C3_stack[:,:,0]),np.real((C3_stack[:,:,2])/np.sqrt(2))+1j*np.imag((C3_stack[:,:,2])/np.sqrt(2)), np.real((C3_stack[:,:,2])/np.sqrt(2))-1j*np.imag((C3_stack[:,:,2])/np.sqrt(2)), np.real((C3_stack[:,:,8])/(2))))
	

def C3_c2(inFolder):

	if os.path.isfile(inFolder+'/C11.bin'):
	    
	    C3_stack = load_C3(inFolder)
	    # T3_stack = C3_T3(C3_stack)
	    
	    if os.path.isdir(inFolder+'\\C2_pp1'):
	        pass
	    else:
	        os.mkdir(inFolder+'\\C2_pp1')
	        
	    infile = inFolder+'/C11.bin'
	    out_file =  inFolder+'\\C2_pp1'+'/C11.bin'
	    write_bin(out_file,np.real(C3_stack[:,:,0]),infile)
	    # infile = inFolder+'/C11.bin'
	    out_file =  inFolder+'\\C2_pp1'+'/C12_real.bin'
	    write_bin(out_file,np.real((C3_stack[:,:,1])/np.sqrt(2)),infile)
	    out_file =  inFolder+'\\C2_pp1'+'/C12_imag.bin'
	    write_bin(out_file,np.imag((C3_stack[:,:,1])/np.sqrt(2)),infile)
	    out_file =  inFolder+'\\C2_pp1'+'/C22.bin'
	    write_bin(out_file,np.real((C3_stack[:,:,4])/(2)),infile)    
	    
	    # copy config file
	    copyfile(inFolder+'/config.txt', inFolder+'\\C2_pp1\\config.txt')
	    
	    print('completed'+inFolder+'\\C2_pp1')
	    # ################## PP2 ###############
	    if os.path.isdir(inFolder+'\\C2_pp2'):
	        pass
	    else:
	        os.mkdir(inFolder+'\\C2_pp2')
	        
	    # infile = inFolder+'/C11.bin'
	    out_file =  inFolder+'\\C2_pp2'+'/C11.bin'
	    write_bin(out_file,np.real(C3_stack[:,:,8]),infile)
	    # infile = inFolder+'/C11.bin'
	    out_file =  inFolder+'\\C2_pp2'+'/C12_real.bin'
	    write_bin(out_file,np.real((C3_stack[:,:,5])/np.sqrt(2)),infile)
	    out_file =  inFolder+'\\C2_pp2'+'/C12_imag.bin'
	    write_bin(out_file,np.imag((C3_stack[:,:,5])/np.sqrt(2)),infile)
	    out_file =  inFolder+'\\C2_pp2'+'/C22.bin'
	    write_bin(out_file,np.real((C3_stack[:,:,4])/(2)),infile)  
	    
	    # copy config file
	    copyfile(inFolder+'/config.txt', inFolder+'\\C2_pp2\\config.txt')
	    
	    print('completed'+inFolder+'\\C2_pp2')
	    ################## PP3 ###############
	    if os.path.isdir(inFolder+'\\C2_pp3'):
	        pass
	    else:
	        os.mkdir(inFolder+'\\C2_pp3')
	        
	    # infile = inFolder+'/C11.bin'
	    out_file =  inFolder+'\\C2_pp3'+'/C11.bin'
	    write_bin(out_file,np.real(C3_stack[:,:,0]),infile)
	    # infile = inFolder+'/C11.bin'
	    out_file =  inFolder+'\\C2_pp3'+'/C12_real.bin'
	    write_bin(out_file,np.real(C3_stack[:,:,2]),infile)
	    out_file =  inFolder+'\\C2_pp3'+'/C12_imag.bin'
	    write_bin(out_file,np.imag(C3_stack[:,:,2]),infile)
	    out_file =  inFolder+'\\C2_pp3'+'/C22.bin'
	    write_bin(out_file,np.real((C3_stack[:,:,8])),infile)  
	    
	    # copy config file
	    copyfile(inFolder+'/config.txt', inFolder+'\\C2_pp3\\config.txt')
	else:
		pass

def S2_T3(S11,S12,S22):

	# Kp- 3-D Pauli feature vector
	Kp = np.array([S11+S22, S11-S22, S12])

	# 3x3 Pauli Coherency Matrix elements

	T11 = np.abs(Kp[0])**2
	T22 = np.abs(Kp[1])**2
	T33 = np.abs(Kp[2])**2

	T12 = Kp[0]*np.conj(Kp[1])
	T13 = Kp[0]*np.conj(Kp[2])
	T23 = Kp[1]*np.conj(Kp[2])
	return np.dstack((T11,T12,T13,np.conjugate(T12),T22,T23,np.conjugate(T13),T23,T33))











