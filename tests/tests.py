import polsartools as pst
from pathlib import Path
import os,shutil


"""
pytest -v tests/tests.py

"""

def utils_processing(T3_folder, window_size=5):
    pst.utils.pauliRGB(T3_folder)
    
def filters_processing(T3_folder, window_size=5):
    pst.rlee(T3_folder, window_size=window_size)

def cp_processing(compact_c2, chi_in=45, window_size=3):
    """ Decompositions """
    pst.mf3cc(compact_c2, chi_in=chi_in, window_size=window_size)
    pst.misomega(compact_c2, chi_in=chi_in, psi_in=0, window_size=window_size)
    """ Descriptors """
    pst.cprvi(compact_c2, chi_in=chi_in, window_size=window_size)
    pst.dopcp(compact_c2, chi_in=chi_in, window_size=window_size)


def full_pol_processing(full_T3, window_size=3):
    
    """ Decompositions"""
    pst.halphafp(full_T3, window_size=window_size)
    pst.neufp(full_T3, window_size=window_size)
    pst.nnedfp(full_T3, window_size=window_size)    
    pst.mf3cf(full_T3, window_size=window_size)
    pst.mf4cf(full_T3, window_size=window_size)
    
    """ Descriptors """   
    pst.dopfp(full_T3, window_size=window_size)
    pst.grvi(full_T3, window_size=window_size)
    pst.rvifp(full_T3, window_size=window_size)
    pst.prvifp(full_T3, window_size=window_size)
    pst.shannon_h_fp(full_T3, window_size=window_size)
    

def dual_cross_pol_processing(dxp_C2, window_size=3):
    """ Descriptors """
    pst.dprvi(dxp_C2, window_size=window_size)
    pst.rvidp(dxp_C2, window_size=window_size)
    pst.prvidp(dxp_C2, window_size=window_size)
    pst.dopdp(dxp_C2, window_size=window_size)



# current_dir = os.getcwd()
# print(current_dir)

# T3_folder = os.path.join(current_dir,'sample_data/full_pol/T3')
# compact_c2 = os.path.join(current_dir,'sample_data/compact_pol/C2_RHV')
# full_T3 = os.path.join(current_dir,'sample_data/full_pol/T3')
# dxp_C2 = os.path.join(current_dir,'sample_data/dual_pol/C2_VVVH')


T3_folder = './tests/sample_data/full_pol/T3'
compact_c2 = './tests/sample_data/compact_pol/C2_RHV'
full_T3 = './tests/sample_data/full_pol/T3'
dxp_C2 = './tests/sample_data/dual_pol/C2_VVVH'

window_size  = 5



# Tests for refined_lee_filter function
def test_filters_processing():
    # We expect this to run without any exceptions or errors
    filters_processing(T3_folder,window_size)
    
    outFolder = os.path.join(os.path.dirname(T3_folder)+ f"_rlee_{window_size}x{window_size}", os.path.basename(T3_folder) )
    
    output_files = [
        os.path.join(outFolder, 'T11.tif'),
        os.path.join(outFolder, 'T22.tif'),
        os.path.join(outFolder, 'T33.tif'),

    ]

    for file_path in output_files:
        assert os.path.exists(file_path), f"{file_path} was not created"
        assert os.path.getsize(file_path) > 0, f"{file_path} is empty"
    
    for file_path in output_files:
        os.remove(file_path)
        print(f"Deleted {file_path}")
    
    shutil.rmtree(os.path.dirname(outFolder)) 
    


def test_cp_processing():
    # Run the function to process
    cp_processing(compact_c2)

    # Check for the existence and size of the output files
    output_files = [
        os.path.join(compact_c2, 'cprvi.tif'),
        os.path.join(compact_c2, 'dopcp.tif'),
        os.path.join(compact_c2, 'Ps_mf3cc.tif'),
        os.path.join(compact_c2, 'Pd_mf3cc.tif'),
        os.path.join(compact_c2, 'Pv_mf3cc.tif'),
        os.path.join(compact_c2, 'Theta_CP_mf3cc.tif'),
        os.path.join(compact_c2, 'Ps_miSOmega.tif'),
        os.path.join(compact_c2, 'Pd_miSOmega.tif'),
        os.path.join(compact_c2, 'Pv_miSOmega.tif')
    ]

    for file_path in output_files:
        assert os.path.exists(file_path), f"{file_path} was not created"
        assert os.path.getsize(file_path) > 0, f"{file_path} is empty"
    
    for file_path in output_files:
        os.remove(file_path)
        print(f"Deleted {file_path}")


def test_utils_processing():
    utils_processing(full_T3)

    # Check for the existence and size of the output files
    output_files = [
            os.path.join(full_T3,'PauliRGB.png'),
            os.path.join(full_T3,'PauliRGB_thumb.png'),
    ]

    for file_path in output_files:
        assert os.path.exists(file_path), f"{file_path} was not created"
        assert os.path.getsize(file_path) > 0, f"{file_path} is empty"
    
    for file_path in output_files:
        os.remove(file_path)
        print(f"Deleted {file_path}")
    
def test_dual_cross_pol_processing():
    dual_cross_pol_processing(dxp_C2)
    
    # Check for the existence and size of the output files
    output_files = [
        os.path.join(dxp_C2,'dprvi.tif'),
        os.path.join(dxp_C2,'rvidp.tif'),
        os.path.join(dxp_C2,'prvidp.tif'),
        os.path.join(dxp_C2,'dopdp.tif'),

    ]

    for file_path in output_files:
        assert os.path.exists(file_path), f"{file_path} was not created"
        assert os.path.getsize(file_path) > 0, f"{file_path} is empty"
    
    for file_path in output_files:
        os.remove(file_path)
        print(f"Deleted {file_path}")

def test_full_pol_processing():
    full_pol_processing(full_T3)

    # Check for the existence and size of the output files
    output_files = [
        
            os.path.join(full_T3,'H_fp.tif'),os.path.join(full_T3,'alpha_fp.tif'),os.path.join(full_T3,'anisotropy_fp.tif'),os.path.join(full_T3,'e1_norm.tif'),os.path.join(full_T3,'e2_norm.tif'),os.path.join(full_T3,'e3_norm.tif'), 
            
            os.path.join(full_T3,'Neu_delta_mod.tif'), os.path.join(full_T3,'Neu_delta_pha.tif'), os.path.join(full_T3,'Neu_psi.tif'), os.path.join(full_T3,'Neu_tau.tif'),            

            os.path.join(full_T3,'NNED_odd.tif'), os.path.join(full_T3,'NNED_dbl.tif'), os.path.join(full_T3,'NNED_vol.tif'),
            
            os.path.join(full_T3,'Ps_mf3cf.tif'),os.path.join(full_T3,'Pd_mf3cf.tif'),os.path.join(full_T3,'Pv_mf3cf.tif'),os.path.join(full_T3,'Theta_FP_mf3cf.tif'),
            
            os.path.join(full_T3,'Ps_mf4cf.tif'),os.path.join(full_T3,'Pd_mf4cf.tif'),os.path.join(full_T3,'Pv_mf4cf.tif'),os.path.join(full_T3,'Pc_mf4cf.tif'),os.path.join(full_T3,'Theta_FP_mf4cf.tif'),os.path.join(full_T3,'Tau_FP_mf4cf.tif'),
            
            os.path.join(full_T3,'H_Shannon.tif'), os.path.join(full_T3,'HI_Shannon.tif'), os.path.join(full_T3,'HP_Shannon.tif'),
            
            os.path.join(full_T3,'dop_fp.tif'),
            os.path.join(full_T3,'grvi.tif'),
            os.path.join(full_T3,'rvifp.tif'),
            os.path.join(full_T3,'prvi_fp.tif'),
    ]

    for file_path in output_files:
        assert os.path.exists(file_path), f"{file_path} was not created"
        assert os.path.getsize(file_path) > 0, f"{file_path} is empty"
    
    for file_path in output_files:
        os.remove(file_path)
        print(f"Deleted {file_path}")