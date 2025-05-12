import polsartools as pst

def main():

    # Define the file paths for compact-pol covariance matrix 
    compact_c2 = r'./sample_data/compact_pol/C2_RHV'
    chi_in = 45  
    window_size = 3  # Window size for processing

    """ Decompositions for Compact Polarimetric Data """
    # Perform Model free 3-component decomposition
    pst.mf3cc(compact_c2, chi_in=chi_in, window_size=window_size)
    # Perform Modified i-SOmega decomposition
    pst.misomega(compact_c2, chi_in=chi_in, psi_in=0, window_size=window_size)

    """ Descriptors for Compact Polarimetric Data """
    # Compute Compact Polarietric Radara Vegetation Index (CpRVI)
    pst.cprvi(compact_c2, chi_in=chi_in, window_size=window_size)
    # Compute Degree of Polarization 
    pst.dopcp(compact_c2, chi_in=chi_in, window_size=window_size)

    # Define file path for full-pol covarance (C3) or Coherency matrix (T3)
    full_T3 = r'./sample_data/full_pol/T3'
    window_size = 3  # Define window size

    """ Decompositions for Full Polarimetric Data """
    # Perform H-alpha FP decomposition
    pst.halphafp(full_T3, window_size=window_size)
    # Perform Neu FP decomposition
    pst.neufp(full_T3, window_size=window_size)
    # Perform Non-negative Eigenvalue decomposition
    pst.nnedfp(full_T3, window_size=window_size)
    # Perform Model free 3-component decomposition
    pst.mf3cf(full_T3, window_size=window_size)
    # Perform Model free 4-component decomposition
    pst.mf4cf(full_T3, window_size=window_size)

    """ Descriptors for Full Polarimetric Data """
    # Compute Barakat degree of polarization  (DoP)
    pst.dopfp(full_T3, window_size=window_size)
    # Compute Generalizerd volume based Radar Vegetation Index
    pst.grvi(full_T3, window_size=window_size)
    # Compute Radar Vegetation Index (RVI)
    pst.rvifp(full_T3, window_size=window_size)
    # Compute Polarimetric Radar Vegetation Index (PRVI)
    pst.prvifp(full_T3, window_size=window_size)

    # Define file path for dual-pol covariance matrix (C2) for corss-pol configuration (VV+VH) or (HH+HV)
    dxp_C2 = r'./sample_data/dual_pol/C2_VVVH'
    window_size = 3  # Define window size

    """ Descriptors for Dual-Polarimetric Data """
    # Compute Dual Polarimetric Radar Vegetation Index (DpRVI)
    pst.dprvi(dxp_C2, window_size=window_size)
    # Compute Radar Vegetation Index for dual-pol data (RVIdp)
    pst.rvidp(dxp_C2, window_size=window_size)
    # Compute Polarimetric Radar Vegetation Index for dual-pol data (PRVIdp)
    pst.prvidp(dxp_C2, window_size=window_size)
    # Compute Barakat degree of polarization for dual-pol data (DoPdp)
    pst.dopdp(dxp_C2, window_size=window_size)

    #%% Refined Lee Polarimetric Speckle Filter
    T3_folder = r'./sample_data/full_pol/T3'
    # Apply the polarimetric refined Lee speckle filter (C3/T3/C2)
    pst.rlee(T3_folder, window_size=5)
    # Generate a Pauli RGB image from the full polarimetric data (C3/T3)
    pst.utils.pauliRGB(T3_folder)

# Run the main function if this script is executed directly
if __name__ == "__main__":
    main()
