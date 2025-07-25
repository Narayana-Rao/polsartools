<p align="center">
  <img src="logo.png" alt=""/>
</p>

## A python package for processing Polarimetric Synthetic Aperture Radar (PolSAR) data.

[![Build](https://github.com/Narayana-Rao/polsartools/actions/workflows/ci.yml/badge.svg)](https://github.com/Narayana-Rao/polsartools/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/polsartools/badge/?version=latest)](https://polsartools.readthedocs.io/en/latest/?badge=latest)

[![image](https://img.shields.io/pypi/v/polsartools.svg)](https://pypi.python.org/pypi/polsartools)
[![](https://static.pepy.tech/badge/polsartools)](https://pepy.tech/project/polsartools)

[![image](https://anaconda.org/bnarayanarao/polsartools/badges/version.svg)](https://anaconda.org/bnarayanarao/polsartools)
[![image](https://anaconda.org/bnarayanarao/polsartools/badges/downloads.svg)](https://anaconda.org/bnarayanarao/polsartools/files)

[![GitHub commits](https://img.shields.io/github/commits-since/Narayana-Rao/polsartools/0.8.svg)](https://GitHub.com/Narayana-Rao/polsartools/commit/)
[![License: GPL 3.0](https://img.shields.io/badge/License-GPL_3.0-green.svg)](https://opensource.org/licenses/gpl-license)



<!-- ## Jointly Developed By

| [![MRSLab](docs/contributions/mrslab_logo.png)](http://mrslab.in) | [![MiRSL](docs/contributions/mirsl_logo.png)](https://www.umass.edu/microwave-remote-sensing) |
|:--:|:--:|
| **Microwave Remote Sensing Lab (MRSLab),** <br> Indian Institute of Technology Bombay, India | **Microwave Remote Sensing Laboratory (MiRSL),** <br> University of Massachusetts Amherst, USA | 
 -->

<h2>Jointly Developed By</h2>

<table align="center">
  <tr>
    <td align="center" style="padding: 20px;">
      <a href="http://mrslab.in">
        <img src="docs/contributions/mrslab_iitb.png" alt="MRSLab" width="200"><br>
        <strong>Microwave Remote Sensing Lab (MRSLab)<br>Indian Institute of Technology Bombay, India</strong>
      </a>
    </td>
    <td align="center" style="padding: 20px;">
      <a href="https://www.umass.edu/microwave-remote-sensing">
        <img src="docs/contributions/mirsl_umass.png" alt="MIRSL" width="250"><br>
        <strong>Microwave Remote Sensing Laboratory (MiRSL)<br>University of Massachusetts Amherst, USA</strong>
      </a>
    </td>
  </tr>
</table>

## 💠General Information
This package generates derived SAR parameters (viz. polarimetric descriptors, vegetation indices, polarimetric decomposition parameters) from various SAR sensors or input polarimetric matrix (S2, C3, T3, Sxy, C2, T2). 

## 💠Installation
1. **Install `gdal` Package**
	
	```bash
	conda install gdal -c conda-forge
	```

2. **Install `polsartools` Package**

	You may choose any of the following options 
	  
	  - a. `pip` (stable release):
	
		```bash
		pip install polsartools
		```
	
	 - b. `conda` (stable release)
	
		```bash
		conda install polsartools -c bnarayanarao
		```
	
	 - c. GitHub (Weekly Build)
	
		```bash
		pip install git+https://github.com/Narayana-Rao/polsartools.git#egg=polsartools
		```
		Use this if you encounter errors like: `AttributeError: module 'polsartools' has no attribute 'xyzabc'`
		
		>**Note for Windows users:** 
		> If installing via GitHub (option c), make sure to install Microsoft C++ build tools first. 
		Download here: https://visualstudio.microsoft.com/visual-cpp-build-tools


## 💠Example Usage

Sample use cases and notebooks are provided at [polsartools-notebooks](https://github.com/Narayana-Rao/polsartools-notebooks) repo. Detailed documentation is available at [polsartools.readthedocs.io](https://polsartools.readthedocs.io/en/latest/) 

## 💠Available functionalities:
* Supported sensors
  * Airborne
    * UAVSAR (GRD, MLC)
    * ISRO ASAR (H5)
    * E-SAR (GTC)
  * Spaceborne
    * NISAR (RSLC, GSLC)
    * RISAT-1 (Compact-pol)
    * RADARSAT-2 (Full-pol)
    * ALOS-2 (Fine Beam Dual-pol (FBD), Quad-pol (HBQ) Level 1.1 CEOS)
  * Planetary 
    * Chandrayaan-II DFSAR (Full-pol)


 * Full-pol :
	* H-A-Alpha decomposition
    * Yamaguchi 4-Component decomposition
    * Shannon Entropy parameters
	* Non-negative Eigen value decomposition
	* Neumann Decomposition 
   * Model free 4-Component decomposition for full-pol data (MF4CF)[[11]](references.md#11)
   * Model free 3-Component decomposition for full-pol data (MF3CF)[[4]](references.md#4)
	* Radar Vegetation Index (RVI) [[8]](references.md#8) 
   * Generalized volume Radar Vegetation Index (GRVI) [[2]](references.md#2)
   * Polarimetric Radar Vegetation Index (PRVI) [[1]](references.md#1)
   * Degree of Polarization (DOP) [[10]](references.md#10) 

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

* Compact-pol : 
  * Model free 3-Component decomposition for compact-pol data (MF3CC) [[4]](references.md#4)
  * Improved S-Omega decomposition for compact-pol data (iS-Omega) [[7]](references.md#7)
  * Compact-pol Radar Vegetation Index (CpRVI)  [[6]](references.md#6)
  * Degree of Polarization (DOP)  [[10]](references.md#10) 

* Polarimetric speckle filters:
  * boxcar
  * refine lee

* Other polarimetric functions
  * Stokes parameters
  * multi-look PolSAR matrix
  * Convert matrices
  * Pauli RGB from full-pol data
  * Psuedo RGB from Dual/compact pol data
  * Polarimetric signatures
  * Full-pol and Dual-pol H-Alpha plot
  * Full-pol H-A-Alpha plot
  * H-Alpha clusters (Full-pol)


## 💠Contributing

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


## 💠Bug Reporting

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


## 💠Feature Requests

We’re always open to suggestions for new features or improvements!

1. **Check existing feature requests:** Please make sure the feature request hasn't already been made in the [Issues section](https://github.com/Narayana-Rao/polsartools/issues).

2. **Submit a feature request:** If it hasn’t been requested already, please open a new issue with the following information:
      * A clear description of the feature.
      * Why you think this feature would be beneficial.
      * Any specific use cases or examples.

[Click here to request a feature](https://github.com/Narayana-Rao/polsartools/issues/new?template=feature_request.md)


## 💠Cite

If you use **`PolSARtools`** in your research or projects, please cite it as follows:


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
