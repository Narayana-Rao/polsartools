<p align="center">
  <img src="logo.png" alt=""/>
</p>

# PolSARtools PyPI package
[![image](https://img.shields.io/pypi/v/polsartools.svg)](https://pypi.python.org/pypi/polsartools)
[![Downloads](https://static.pepy.tech/badge/polsartools)](https://pepy.tech/project/polsartools)
[![Documentation Status](https://readthedocs.org/projects/polsartools/badge/?version=latest)](https://polsartools.readthedocs.io/en/latest/?badge=latest)
[![License: GPL 3.0](https://img.shields.io/badge/License-GPL_3.0-green.svg)](https://opensource.org/licenses/gpl-license)


<!-- [![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FNarayana-Rao%2Fpolsartools&count_bg=%2379C83D&title_bg=%23555555&icon=cliqz.svg&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com) -->

> **Cite:** Bhogapurapu, N., Dey, S., Mandal, D., Bhattacharya, A. and Rao, Y.S., 2021. PolSAR tools: A QGIS plugin for generating SAR descriptors. Journal of Open Source Software, 6(60), p.2970. doi:  [10.21105/joss.02970](https://doi.org/10.21105/joss.02970)  
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

## General Information
-------------------
This package generates derived SAR parameters (viz. vegetation indices, polarimetric decomposition parameters) from input polarimetric matrix (C3, T3, C2, T2). The input data needs to be in [PolSARpro](https://earth.esa.int/web/polsarpro/home)/[ENVI](https://www.l3harrisgeospatial.com/Software-Technology/ENVI) format (\*.bin and \*.hdr). 

## Installation
```
pip install polsartools
```
 
#### Prerequesites  ```gdal, Numpy```
 

 
#### gdal installation error fix ```conda install gdal```
 
## Example Usage

#### More sample use cases are provided in [examples](https://github.com/Narayana-Rao/polsartools/tree/main/examples) folder

```python
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

```


## Available functionalities:
-----------------------------
 * Full-pol :
	* H-Alpha decomposition
	* Non-negative Eigen value decomposition
	* Neumann Decomposition 
    * Model free 4-Component decomposition for full-pol data (MF4CF)[[11]](#11)
    * Model free 3-Component decomposition for full-pol data (MF3CF)[[4]](#4)
	* Radar Vegetation Index (RVI) [[8]](#8) 
    * Generalized volume Radar Vegetation Index (GRVI) [[2]](#2)
    * Polarimetric Radar Vegetation Index (PRVI) [[1]](#1)
    * Degree of Polarization (DOP) [[10]](#10) 

  * Compact-pol : 
    * Model free 3-Component decomposition for compact-pol data (MF3CC) [[4]](#4)
    * Improved S-Omega decomposition for compact-pol data (iS-Omega) [[7]](#7)
    * Compact-pol Radar Vegetation Index (CpRVI)  [[6]](#6)
    * Degree of Polarization (DOP)  [[10]](#10) 

  * Dual-pol:
	* Dual-pol Radar Vegetation Index (DpRVI) [[5]](#5)
    * Dual-pol Radar Vegetation Index for GRD data (DpRVIc) [[12]](#12)
	* Radar Vegetation Index (RVI) [[9]](#9)
    * Degree of Polarization (DOP) [[10]](#10) 
    * Polarimetric Radar Vegetation Index (PRVI) [[1]](#1)
    * Dual-pol descriptors [[13]](#13)


## Contributions
1) Contribute to the software

    [Contribution guidelines for this project](help/CONTRIBUTING.md)


2) Report issues or problems with the package
	
	Please raise your issues here : <https://github.com/Narayana-Rao/polsartools/issues>

<!-- 3) Seek support

	Please write to us: <bnarayanarao@iitb.ac.in>  -->

## References
-------------
<a id="1">[1]</a> 
Chang, J.G., Shoshany, M. and Oh, Y., 2018. Polarimetric Radar Vegetation Index for Biomass Estimation in Desert Fringe Ecosystems. IEEE Transactions on Geoscience and Remote Sensing, 56(12), pp.7102-7108.

<a id="2">[2]</a> 
Ratha, D., Mandal, D., Kumar, V., McNairn, H., Bhattacharya, A. and Frery, A.C., 2019. A generalized volume scattering model-based vegetation index from polarimetric SAR data. IEEE Geoscience and Remote Sensing Letters, 16(11), pp.1791-1795.

<a id="3">[3]</a> 
Mandal, D., Kumar, V., Ratha, D., J. M. Lopez-Sanchez, A. Bhattacharya, H. McNairn, Y. S. Rao, and K. V. Ramana, 2020. Assessment of rice growth conditions in a semi-arid region of India using the Generalized Radar Vegetation Index derived from RADARSAT-2 polarimetric SAR data, Remote Sensing of Environment, 237: 111561.

<a id="4">[4]</a> 
Dey, S., Bhattacharya, A., Ratha, D., Mandal, D. and Frery, A.C., 2020. Target Characterization and Scattering Power Decomposition for Full and Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing.

<a id="5">[5]</a> 
Mandal, D., Kumar, V., Ratha, D., Dey, S., Bhattacharya, A., Lopez-Sanchez, J.M., McNairn, H. and Rao, Y.S., 2020. Dual polarimetric radar vegetation index for crop growth monitoring using sentinel-1 SAR data. Remote Sensing of Environment, 247, p.111954.

<a id="6">[6]</a> 
Mandal, D., Ratha, D., Bhattacharya, A., Kumar, V., McNairn, H., Rao, Y.S. and Frery, A.C., 2020. A Radar Vegetation Index for Crop Monitoring Using Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing, 58 (9), pp. 6321-6335.

<a id="7">[7]</a> 
V. Kumar, D. Mandal, A. Bhattacharya, and Y. S. Rao, 2020. Crop Characterization Using an Improved Scattering Power Decomposition Technique for Compact Polarimetric SAR Data. International Journal of Applied Earth Observations and Geoinformation, 88: 102052.

<a id="8">[8]</a> 
Kim, Y. and van Zyl, J.J., 2009. A time-series approach to estimate soil moisture using polarimetric radar data. IEEE Transactions on Geoscience and Remote Sensing, 47(8), pp.2519-2527.

<a id="9">[9]</a> 
Trudel, M., Charbonneau, F. and Leconte, R., 2012. Using RADARSAT-2 polarimetric and ENVISAT-ASAR dual-polarization data for estimating soil moisture over agricultural fields. Canadian Journal of Remote Sensing, 38(4), pp.514-527.

<a id="10">[10]</a> 
Barakat, R., 1977. Degree of polarization and the principal idempotents of the coherency matrix. Optics Communications, 23(2), pp.147-150.

<a id="11">[11]</a> S. Dey, A. Bhattacharya, A. C. Frery, C. Lopez-Martinez and Y. S. Rao, "A Model-free Four Component Scattering Power Decomposition for Polarimetric SAR Data," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 2021. doi: [10.1109/JSTARS.2021.3069299](https://doi.org/10.1109/JSTARS.2021.3069299). 

<a id="12">[12]</a> Bhogapurapu, N., Dey, S., Mandal, D., Bhattacharya, A., Karthikeyan, L., McNairn, H. and Rao, Y.S., 2022. Soil moisture retrieval over croplands using dual-pol L-band GRD SAR data. Remote Sensing of Environment, 271, p.112900. 

<a id="13">[13]</a>Bhogapurapu, N., Dey, S., Bhattacharya, A., Mandal, D., Lopez-Sanchez, J.M., McNairn, H., López-Martínez, C. and Rao, Y.S., 2021. Dual-polarimetric descriptors from Sentinel-1 GRD SAR data for crop growth assessment. ISPRS Journal of Photogrammetry and Remote Sensing, 178, pp.20-35.


