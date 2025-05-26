import os

def dcpt2files(infolder):
    if os.path.isfile(os.path.join(infolder,"T11.bin")):
        input_filepaths = [
            os.path.join(infolder, "T11.bin"), 
            os.path.join(infolder, "T12_real.bin"),
            os.path.join(infolder, "T12_imag.bin"),
            os.path.join(infolder, "T22.bin")
        ]
    else:
        raise(f"Invalid T2 folder!!")
    return input_filepaths