import os

def get_filter_io_paths(infolder,  window_size, filter_type=None):
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
    # if outname is None:
    # outFolder = os.path.join(os.path.dirname(infolder), os.path.basename(infolder) + f"_{window_size}x{window_size}")
    if filter_type =='rlee':
        outFolder = os.path.join(os.path.dirname(infolder)+ f"_{filter_type}_{window_size}x{window_size}", os.path.basename(infolder) )
    elif filter_type == 'boxcar':
        outFolder = os.path.join(os.path.dirname(infolder)+ f"_boxcar_{window_size}x{window_size}", os.path.basename(infolder) )
    else:
        outFolder = os.path.join(os.path.dirname(infolder)+ f"_{window_size}x{window_size}", os.path.basename(infolder) )
    
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
