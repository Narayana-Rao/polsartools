import os
import numpy as np
from polsartools.utils.utils import process_chunks_parallel, time_it, conv2d
from polsartools.utils.convert_matrices import C3_T3_mat

def get_filter_io_paths(infolder, outname, window_size, filter_type=None):
    # Determine the input filepaths based on available matrix types
    input_filepaths = []
    matrix_type = None  # To identify if it's C2, C3, T2, or T3

    if os.path.isfile(os.path.join(infolder, "T11.bin")) and os.path.isfile(os.path.join(infolder, "T33.bin")):
        matrix_type = "T3"
        input_filepaths = [
            os.path.join(infolder, "T11.bin"),
            os.path.join(infolder, 'T12_real.bin'), os.path.join(infolder, 'T12_imag.bin'),
            os.path.join(infolder, 'T13_real.bin'), os.path.join(infolder, 'T13_imag.bin'),
            os.path.join(infolder, "T22.bin"),
            os.path.join(infolder, 'T23_real.bin'), os.path.join(infolder, 'T23_imag.bin'),
            os.path.join(infolder, "T33.bin"),
        ]
    elif os.path.isfile(os.path.join(infolder, "C11.bin")) and os.path.isfile(os.path.join(infolder, "C33.bin")):
        matrix_type = "C3"
        input_filepaths = [
            os.path.join(infolder, "C11.bin"),
            os.path.join(infolder, 'C12_real.bin'), os.path.join(infolder, 'C12_imag.bin'),
            os.path.join(infolder, 'C13_real.bin'), os.path.join(infolder, 'C13_imag.bin'),
            os.path.join(infolder, "C22.bin"),
            os.path.join(infolder, 'C23_real.bin'), os.path.join(infolder, 'C23_imag.bin'),
            os.path.join(infolder, "C33.bin"),
        ]
    elif os.path.isfile(os.path.join(infolder, "C11.bin")) and os.path.isfile(os.path.join(infolder, "C22.bin")) \
            and not os.path.isfile(os.path.join(infolder, "C33.bin")):
        matrix_type = "C2"
        input_filepaths = [
            os.path.join(infolder, "C11.bin"),
            os.path.join(infolder, 'C12_real.bin'), os.path.join(infolder, 'C12_imag.bin'),
            os.path.join(infolder, "C22.bin"),
        ]
    elif os.path.isfile(os.path.join(infolder, "T11.bin")) and os.path.isfile(os.path.join(infolder, "T22.bin")) \
            and not os.path.isfile(os.path.join(infolder, "T33.bin")):
        matrix_type = "T2"
        input_filepaths = [
            os.path.join(infolder, "T11.bin"),
            os.path.join(infolder, 'T12_real.bin'), os.path.join(infolder, 'T12_imag.bin'),
            os.path.join(infolder, "T22.bin"),
        ]
    else:
        raise ValueError("Invalid C2/C3 or T2/T3 folder!")

    # Determine the output folder and filepaths
    output_filepaths = []
    if outname is None:
        outFolder = os.path.join(os.path.dirname(infolder), os.path.basename(infolder) + f"_{window_size}x{window_size}")
        os.makedirs(outFolder, exist_ok=True)
        
        # Only use the first letter of the matrix type (C for C3, T for T3)
        matrix_prefix = matrix_type[0]

        if matrix_type == "C3" or matrix_type == "T3":
            output_filepaths = [
                os.path.join(outFolder, f"{matrix_prefix}11.bin"),
                os.path.join(outFolder, f"{matrix_prefix}12_real.bin"), os.path.join(outFolder, f"{matrix_prefix}12_imag.bin"),
                os.path.join(outFolder, f"{matrix_prefix}13_real.bin"), os.path.join(outFolder, f"{matrix_prefix}13_imag.bin"),
                os.path.join(outFolder, f"{matrix_prefix}22.bin"),
                os.path.join(outFolder, f"{matrix_prefix}23_real.bin"), os.path.join(outFolder, f"{matrix_prefix}23_imag.bin"),
                os.path.join(outFolder, f"{matrix_prefix}33.bin"),
            ]
        elif matrix_type == "C2" or matrix_type == "T2":
            output_filepaths = [
                os.path.join(outFolder, f"{matrix_prefix}11.bin"),
                os.path.join(outFolder, f"{matrix_prefix}12_real.bin"), os.path.join(outFolder, f"{matrix_prefix}12_imag.bin"),
                os.path.join(outFolder, f"{matrix_prefix}22.bin"),
            ]
    
    return input_filepaths, output_filepaths

@time_it
def boxcar(infolder, outname=None, chi_in=0, psi_in=0, window_size=3, write_flag=True, max_workers=None):
    
    input_filepaths, output_filepaths = get_filter_io_paths(infolder, outname, window_size, filter_type="boxcar")

    # Process chunks in parallel
    num_outputs = len(output_filepaths)
    process_chunks_parallel(input_filepaths, list(output_filepaths), window_size=window_size, write_flag=write_flag,
                            processing_func=process_chunk_boxcar, block_size=(512, 512), max_workers=max_workers, 
                            num_outputs=num_outputs)

def process_chunk_boxcar(chunks, window_size, input_filepaths, *args):
    filtered_chunks = []
    for i in range(len(chunks)):
        img = np.array(chunks[i])
        kernel = np.ones((window_size, window_size), np.float32) / (window_size * window_size)
        filtered_chunks.append(conv2d(img, kernel))
    return filtered_chunks

    # # Handle T3 case (9 elements)
    # if 'T11' in input_filepaths[0] and 'T22' in input_filepaths[5] and 'T33' in input_filepaths[8]:
    #     t11_T1 = np.array(chunks[0])
    #     t12_T1 = np.array(chunks[1]) + 1j * np.array(chunks[2])
    #     t13_T1 = np.array(chunks[3]) + 1j * np.array(chunks[4])
    #     t21_T1 = np.conj(t12_T1)
    #     t22_T1 = np.array(chunks[5])
    #     t23_T1 = np.array(chunks[6]) + 1j * np.array(chunks[7])
    #     t31_T1 = np.conj(t13_T1)
    #     t32_T1 = np.conj(t23_T1)
    #     t33_T1 = np.array(chunks[8])

    #     T_T1 = np.array([[t11_T1, t12_T1, t13_T1],
    #                      [t21_T1, t22_T1, t23_T1],
    #                      [t31_T1, t32_T1, t33_T1]])

    # # Handle C3 case (9 elements)
    # elif 'C11' in input_filepaths[0] and 'C22' in input_filepaths[5] and 'C33' in input_filepaths[8]:
    #     C11 = np.array(chunks[0])
    #     C12 = np.array(chunks[1]) + 1j * np.array(chunks[2])
    #     C13 = np.array(chunks[3]) + 1j * np.array(chunks[4])
    #     C21 = np.conj(C12)
    #     C22 = np.array(chunks[5])
    #     C23 = np.array(chunks[6]) + 1j * np.array(chunks[7])
    #     C31 = np.conj(C13)
    #     C32 = np.conj(C23)
    #     C33 = np.array(chunks[8])

    #     T_T1 = np.array([[C11, C12, C13],
    #                      [C21, C22, C23],
    #                      [C31, C32, C33]])

    # # Handle T2 case (6 elements)
    # elif 'T11' in input_filepaths[0] and 'T22' in input_filepaths[3]:
    #     t11_T1 = np.array(chunks[0])
    #     t12_T1 = np.array(chunks[1]) + 1j * np.array(chunks[2])
    #     t21_T1 = np.conj(t12_T1)
    #     t22_T1 = np.array(chunks[3])

    #     T_T1 = np.array([[t11_T1, t12_T1],
    #                      [t21_T1, t22_T1]])

    # # Handle C2 case (6 elements)
    # elif 'C11' in input_filepaths[0] and 'C22' in input_filepaths[3]:
    #     C11 = np.array(chunks[0])
    #     C12 = np.array(chunks[1]) + 1j * np.array(chunks[2])
    #     C21 = np.conj(C12)
    #     C22 = np.array(chunks[3])

    #     T_T1 = np.array([[C11, C12],
    #                      [C21, C22]])

    # # Perform boxcar filtering for window size > 1
    # if window_size > 1:
    #     kernel = np.ones((window_size, window_size), np.float32) / (window_size * window_size)

    #     filtered_chunks = []
    #     for i in range(T_T1.shape[0]):
    #         for j in range(T_T1.shape[1]):
    #             filtered_chunks.append(conv2d(T_T1[i, j, :, :], kernel))

    #     # Return filtered outputs for further writing
    #     return filtered_chunks
    # else:
    #     return T_T1.flatten()
