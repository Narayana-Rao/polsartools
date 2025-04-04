import numpy as np
# import process_arrays,refined_lee  # The compiled module
import polsartools
from functools import wraps
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import time

def time_it(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        try:
            result = func(*args, **kwargs)  # Call the original function
            end_time = time.time()
            print(f"Execution time for {func.__name__}: {end_time - start_time:.2f} seconds")
            return result
        except Exception as e:
            # If an exception occurs, re-raise it but don't print the execution time
            raise e
    return wrapper
def process_chunk_cprvi(chunks, window_size,input_filepaths,chi_in,psi_in):

    kernel = np.ones((window_size,window_size),np.float32)/(window_size*window_size)
    c11_T1 = np.array(chunks[0])
    c12_T1 = np.array(chunks[1])+1j*np.array(chunks[2])
    c21_T1 = np.conj(c12_T1)
    c22_T1 = np.array(chunks[3])

    ncols,nrows = np.shape(c11_T1)
    
    fp22 = np.zeros((ncols,nrows))
    l_lambda = np.zeros((ncols,nrows))

    wsi=wsj=window_size

    inci=int(np.fix(wsi/2)) # Up & down movement margin from the central row
    incj=int(np.fix(wsj/2)) # Left & right movement from the central column
    # % Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

    starti=int(np.fix(wsi/2)) # Starting row for window processing
    startj=int(np.fix(wsj/2)) # Starting column for window processing

    stopi= int(nrows-inci)-1 # Stop row for window processing
    stopj= int(ncols-incj)-1 # Stop column for window processing

    for ii in np.arange(startj,stopj+1):

        for jj in np.arange(starti,stopi+1):
            
            C11c = np.nanmean(c11_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            C12c = np.nanmean(c12_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            C21c = np.nanmean(c21_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample
            C22c = np.nanmean(c22_T1[ii-inci:ii+inci+1,jj-incj:jj+incj+1])#i sample

            C0 = np.array([[C11c,C12c], [C21c, C22c]]);
            
            # %% GD_VI -- VV-VH/VV-HH
            if np.isnan(np.real(C0)).any() or np.isinf(np.real(C0)).any() or np.isneginf(np.real(C0)).any():
                C0 = np.array([[0,0],[0,0]])

            # Stokes Parameter
            s0 = C11c + C22c;
            s1 = C11c - C22c;
            s2 = (C12c + C21c);

            if (chi_in >= 0):
                s3 = (1j*(C12c - C21c)); # The sign is according to RC or LC sign !!
            if (chi_in < 0):
                s3 = -(1j*(C12c - C21c)); # The sign is according to RC or LC sign !!
            
            k11 = s0
            k12 = 0
            k13 = s2
            k14 = 0
            k21 = k12 
            k22 = 0
            k23 = 0
            k24 = s1
            k31 = k13
            k32 = k23
            k33 = 0 
            k34 = 0;
            k41 = k14; 
            k42 = k24; 
            k43 = k34; 
            k44 = s3;

            K_T = 0.5*np.array([[k11,k12,k13,k14], [k21,k22,k23,k24], 
                [k31, k32, k33, k34], [k41,k42,k43,k44]])       

            # Stokes vector child products
            SC = ((s0)-(s3))/2;
            OC = ((s0)+(s3))/2;
            
            min_sc_oc = min(SC,OC);
            max_sc_oc = max(SC,OC);                
            
            K_depol = np.array([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]);
            
            # GD_DEPOL
            
            num1 = np.matmul(K_T.T,K_depol);
            num = np.trace(num1);
            den1 = np.sqrt(abs(np.trace(np.matmul(K_T.T,K_T))));
            den2 = np.sqrt(abs(np.trace(np.matmul(K_depol.T,K_depol))));
            den = den1*den2;
            
            temp_aa = np.real(2*np.arccos(num/den)*180/np.pi);
            GD_t1_depol = np.real(temp_aa/180);
            
                                                           
            l_lambda[ii,jj] = (3/2)*GD_t1_depol;
            
            #GD_VI -- RH-RV/LH-LV
            
            fp22[ii,jj] = (min_sc_oc/max_sc_oc);
            
            # print(ii,jj)
            # print(num,den1,den2,den)
            # print(C11c,C12c,'SC',SC,'OC',OC,'GD',GD_t1_depol, fp22[ii,jj], 'll',l_lambda[ii][jj])

    vi_c = np.real((1 - l_lambda)*np.power(fp22, 2*l_lambda))

    return vi_c.astype(np.float32)




# Example input: List of 2D NumPy arrays
arr1 = np.random.rand(100, 100).astype(np.float32)
arr2 = np.random.rand(100, 100).astype(np.float32)
arr3 = np.random.rand(100, 100).astype(np.float32)
arr4 = np.random.rand(100, 100).astype(np.float32)
# arr2 = np.random.rand(200, 200).astype(np.float32)




num_nans = 500  # Number of NaNs to insert
nan_indices = np.random.choice(arr1.size, num_nans, replace=False) 

nan_indices_2d = np.unravel_index(nan_indices, arr1.shape)

arr1[nan_indices_2d] = np.nan
arr2[nan_indices_2d] = np.nan
arr3[nan_indices_2d] = np.nan
arr4[nan_indices_2d] = np.nan




arr_list = [arr1, arr2, arr3, arr4]

# Process the arrays
# processed_arr_list = process_arrays.process_arrays(arr_list)

# # Print the results
# for i, processed_arr in enumerate(processed_arr_list):
#     print(f"Processed array {i}:\n{processed_arr}")

print("Input array statistics:")
for i, arr in enumerate(arr_list):
    print(f"Array {i}: mean={np.nanmean(arr)}, min={np.nanmin(arr)}, max={np.nanmax(arr)}")

# processed_arr_list = testpip.refined_lee.refined_lee(arr_list,7)

import time
t0 = time.time()

cprvi_cpp_out = polsartools.cprvicpp.process_chunk_cprvicpp(arr_list,5,['sdfgfd0','sdfgfd1','sdfgfd2','sdfgfd3'],-45,0)
print(f"Output array statistics: mean={np.nanmean(cprvi_cpp_out)}, min={np.nanmin(cprvi_cpp_out)}, max={np.nanmax(cprvi_cpp_out)}")
print('Time for cprvi_cpp:', time.time() - t0)
t1 = time.time()
cprvi_py_out = process_chunk_cprvi(arr_list,5,['sdfgfd0','sdfgfd1','sdfgfd2','sdfgfd3'],-45,0)
print(f"Output array statistics: mean={np.nanmean(cprvi_py_out)}, min={np.nanmin(cprvi_py_out)}, max={np.nanmax(cprvi_py_out)}")
print('Time for cprvi_py:', time.time() - t1)

print('Done!!')

are_equal = np.allclose(cprvi_cpp_out, cprvi_py_out,rtol=1e-5, atol=1e-8)
print(are_equal) 
# print(cprvi_cpp_out)
# print(cprvi_py_out)

# print(cprvi_cpp_out-cprvi_py_out)
# Print the results
# for i, processed_arr in enumerate(processed_arr_list):
#     print(f"Processed array {i}:\n{processed_arr}")
# print(processed_arr_list.shape)
# for i, processed_arr in enumerate(processed_arr_list):
    # print(f"Processed array {i}: mean={np.nanmean(processed_arr)}, min={np.nanmin(processed_arr)}, max={np.nanmax(processed_arr)}")