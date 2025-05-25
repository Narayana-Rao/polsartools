<p align="center">
  <img src="logo.png" alt=""/>
</p>

# PolSARtools PyPI package
[![Build](https://github.com/Narayana-Rao/polsartools/actions/workflows/ci.yml/badge.svg)](https://github.com/Narayana-Rao/polsartools/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/polsartools/badge/?version=latest)](https://polsartools.readthedocs.io/en/latest/?badge=latest)
[![image](https://img.shields.io/pypi/v/polsartools.svg)](https://pypi.python.org/pypi/polsartools)
[![GitHub commits](https://img.shields.io/github/commits-since/Narayana-Rao/polsartools/0.6.1.svg)](https://GitHub.com/Narayana-Rao/polsartools/commit/)
[![Downloads](https://static.pepy.tech/badge/polsartools)](https://pepy.tech/project/polsartools)
[![License: GPL 3.0](https://img.shields.io/badge/License-GPL_3.0-green.svg)](https://opensource.org/licenses/gpl-license)


## General Information

This package generates derived SAR parameters (viz. polarimetric descriptors, vegetation indices, polarimetric decomposition parameters) from different SAR sensors or input polarimetric matrix (C3, T3,C2, T2). 

<!-- The input data needs to be in [PolSARpro](https://earth.esa.int/web/polsarpro/home)/[ENVI](https://www.l3harrisgeospatial.com/Software-Technology/ENVI) format (\*.bin and \*.hdr).  -->

## Installation
### 1. Install **`gdal`** Package

```bash
conda install gdal -c conda-forge
```

### 2. Install **`polsartools`** Package

You can install it using `pip`:

```bash
pip install polsartools
```

If you encounter an error like "function not found" or need the latest version, you can install the updated version directly from GitHub:

```bash
pip install git+https://github.com/Narayana-Rao/polsartools.git#egg=polsartools
```

## Available functionalities:
* Supported sensors
  * UAVSAR (GRD, MLC)
  * NISAR (RSLC, GSLC)
  * RADARSAT-2 (Full-pol)
  * ALOS-2 (Fine Beam Dual-pol (FBD) Level 1.1 CEOS)
  * Chandrayaan-II DFSAR (Full-pol)

 * Full-pol :
	* H-Alpha decomposition
    * Shannon Entropy parameters
	* Non-negative Eigen value decomposition
	* Neumann Decomposition 
   * Model free 4-Component decomposition for full-pol data (MF4CF)[[11]](references.md#11)
   * Model free 3-Component decomposition for full-pol data (MF3CF)[[4]]references.md(#4)
	* Radar Vegetation Index (RVI) [[8]](references.md#8) 
   * Generalized volume Radar Vegetation Index (GRVI) [[2]](references.md#2)
   * Polarimetric Radar Vegetation Index (PRVI) [[1]](references.md#1)
   * Degree of Polarization (DOP) [[10]](references.md#10) 

* Compact-pol : 
  * Model free 3-Component decomposition for compact-pol data (MF3CC) [[4]](references.md#4)
  * Improved S-Omega decomposition for compact-pol data (iS-Omega) [[7]](references.md#7)
  * Compact-pol Radar Vegetation Index (CpRVI)  [[6]](references.md#6)
  * Degree of Polarization (DOP)  [[10]](references.md#10) 

* Dual-pol:
  * H-Alpha parameters
  * Shannon Entropy parameters
  * Dual-pol Radar Vegetation Index (DpRVI) [[5]](references.md#5)
  * Dual-pol Radar Vegetation Index for GRD data (DpRVIc) [[12]](references.md#12)
  * Radar Vegetation Index (RVI) [[9]](references.md#9)
  * Degree of Polarization (DOP) [[10]](references.md#10) 
  * Polarimetric Radar Vegetation Index (PRVI) [[1]](references.md#1)
  * Dual-pol descriptors [[13]](references.md#13)
  * Model free 3-Component decomposition for dual-copol data (MF3CD)

* Polarimetric speckle filters:
  * boxcar
  * refine lee

## Example Usage

#### Sample use cases and notebooks are provided at [polsartools-notebooks](https://github.com/Narayana-Rao/polsartools-notebooks) repo. Detailed documentation is available at [polsartools.readthedocs.io](https://polsartools.readthedocs.io/en/latest/) 
<!-- 
```python
import polsartools as pst

def main():

    #Extract data 
    # Extract C2 data from NISAR GSLC data
    pst.nisar_gslc(r"NISAR_L2_PR_GSLC_001_030_A_019_2800_SHNA_A_20081012T060911_20081012T060925_D00404_N_F_J_001.h5")
    pst.nisar_rslc(r"NISAR_L1_PR_RSLC_001_030_A_019_2000_SHNA_A_20081012T060910_20081012T060926_D00402_N_F_J_001.h5")
    
    annFile = r"./winnip_31606_12049_004_120627_L090_CX_03_grd/winnip_31606_12049_004_120627_L090_CX_03.ann"
    #The following function extracts C3 matrix from UAVSAR GRD data
    pst.uavsar_grd(annFile)
    
    #The following function extracts C3 matrix from UAVSAR MLC data
    pst.uavsar_mlc(annFile)
    
    Extracted_fbd_path = r""
    #The following function extracts C2 matrix from ALOS-2 FBD level 1.1 CEOS data
    pst.alos2_fbd_l11(Extracted_fbd_path)

    Extracted_cy2_path = r""
    #The following function extracts T3/C3 matrix from Chandrayaan-2 DFSAR data
    pst.chyaan2_fp(Extracted_cy2_path,matrix='T3')


    
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
    
    dcp_T2 = r'./sample_data/dual_copol/T2'
    # Derive the dual-copol model-free 3-component decomposition powers from T2 matrix
    pst.mf3cdp(dcp_T2, window_size=window_size)


    #%% Refined Lee Polarimetric Speckle Filter
    T3_folder = r'./sample_data/full_pol/T3'
    # Apply the polarimetric refined Lee speckle filter (C3/T3/C2)
    pst.rlee(T3_folder, window_size=5)
    # Generate a Pauli RGB image from the full polarimetric data (C3/T3)
    pst.utils.pauliRGB(T3_folder)

# Run the main function if this script is executed directly
if __name__ == "__main__":
    main()

``` -->



## Contributing

We welcome contributions! Whether it's fixing bugs, adding new features, or improving documentation, your help is greatly appreciated.

### How to Contribute
1. **Fork the repository** - Fork this repository to your GitHub account.

2. **Clone your fork** - Clone the repository to your local machine:

    ```bash
    git clone https://github.com/Narayana-Rao/polsartools.git
    ```

3. Create a branch - Create a new branch for your changes:

    ```bash 
    git checkout -b feature-branch
    ```
4. **Make changes** - Implement your changes or additions.

5. **Test your changes** - Run the tests to ensure that your changes don’t break anything.

6. **Commit and push** - Commit your changes and push them to your fork:
    ```bash
        git commit -am "Description of changes"
        git push origin feature-branch
    ```
7. **Create a Pull Request** - Open a pull request to the main repository with a clear description of the changes.

<!-- For more detailed guidelines on contributing, see the CONTRIBUTING.md (if available). -->


## Bug Reporting

If you encounter a bug or issue, please follow these steps to report it:

1. Check the existing issues: Before submitting a new bug report, check if the issue has already been reported in the [Issues section](https://github.com/Narayana-Rao/polsartools/issues).

2. Submit a bug report: If the issue hasn’t been reported, please open a new issue and include the following information:
    * A clear description of the problem.
    * Steps to reproduce the issue.
    * Expected vs actual behavior.
    * Any error messages or stack traces.
    * Relevant code snippets or files if possible.
    * Version of `polsartools` and Python you're using.

[Click here to report a bug](https://github.com/Narayana-Rao/polsartools/issues/new?template=bug_report.md)


## Feature Requests

We’re always open to suggestions for new features or improvements!

1. **Check existing feature requests:** Please make sure the feature request hasn't already been made in the [Issues section](https://github.com/Narayana-Rao/polsartools/issues).

2. **Submit a feature request:** If it hasn’t been requested already, please open a new issue with the following information:
      * A clear description of the feature.
      * Why you think this feature would be beneficial.
      * Any specific use cases or examples.

[Click here to request a feature](https://github.com/Narayana-Rao/polsartools/issues/new?template=feature_request.md)


## Cite

If you use **PolsarTools** in your research or projects, please cite it as follows:


> Bhogapurapu, N., Dey, S., Mandal, D., Bhattacharya, A. and Rao, Y.S., 2021. PolSAR tools: A QGIS plugin for generating SAR descriptors. Journal of Open Source Software, 6(60), p.2970. doi:  [10.21105/joss.02970](https://doi.org/10.21105/joss.02970)  


```bibtex
@article{bhogapurapu2021polsar,
  title={PolSAR tools: A QGIS plugin for generating SAR descriptors},
  author={Bhogapurapu, Narayanarao and Dey, Subhadip and Mandal, Dipankar and Bhattacharya, Avik and Rao, YS},
  journal={Journal of Open Source Software},
  volume={6},
  number={60},
  pages={2970},
  year={2021},
  doi= {10.21105/joss.02970}
}

```