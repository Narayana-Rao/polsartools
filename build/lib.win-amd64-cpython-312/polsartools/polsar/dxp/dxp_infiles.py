import os

def dxpc2files(infolder):
    if os.path.isfile(os.path.join(infolder,"C11.bin")):
        input_filepaths = [
        os.path.join(infolder, "C11.bin"), 
        os.path.join(infolder, "C12_real.bin"),
        os.path.join(infolder, "C12_imag.bin"),
        os.path.join(infolder, "C22.bin")
    ]
    else:
        print(f"Invalid C2 folder!!")
        
    return input_filepaths