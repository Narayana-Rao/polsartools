---
title: 'PolSAR tools: A python package to process and analyse Polarimetric Synthetic Aperture Radar (PolSAR) data'
tags:
  - SAR
  - PolSAR
  - Speckle filtering
  - polarimetric signature
  - Vegetation indices
  - Polarimetric decompositions

authors:
  - name: Narayanarao Bhogapurapu
    orcid: 0000-0002-6496-7283
    affiliation: "1, 2" 

  - name: Paul Siqueira
    orcid: 0000-0001-5781-8282
    affiliation: "1, 2"  
  
  - name: Avik Bhattacharya
    orcid: 0000-0001-6720-6108
    affiliation: "3" 




affiliations:
 - name: Microwave Remote Sensing Laboratory, University of Massachusetts Amherst, USA
   index: 1
 - name: Division of Geological and Planetary Sciences, California Institute of Technology, USA
   index: 2
 - name: Microwave Remote Sensing Lab, Indian Institute of Technology Bombay, India
   index: 3

date: 4 July 2025
bibliography: refs.bib

---

# Intoduction
Synthetic Aperture Radar has been proven to be a highly reliable earth observation (EO) remote sensing technique due to its all-weather and day-night capabilities. Since its initial development and use on an airborne platform seven decades ago, it has significantly evolved in terms of spatial resolution, observing frequencies and applications. It took nearly four decades since then to have a full/quad polarimetric SAR (PolSAR) instrument in space for earth observation [@ulaby1981microwave;@jordan2002sir]. Since then, PolSAR has become a critical data source in a wide variety of earth observation applications [@lee2017polarimetric;@cloude2010polarisation;@van2011synthetic]. The wide reach of SAR data applications has become possible by several airborne and spaceborne SAR instruments and thier open data. As the world moving towards open data and open science, there are serval operational SAR satellite missions (e.g., Sentinel-1, ALOS-2/4, EOS-04, BIOMASS, etc) that have been collecting and distributing data openly. Additionally, the poposed future missions such as NISAR[@rosen2017global], ROSE-L[@davidson2021rose] and Sentinel-1 constellation that are planned to collect Petabytes of EO data. For an efficient use of these PolSAR data, it is essential to process it and derive interpretable polarimetric parameters to support further downstream applications. 


# Statement of need

With the increasing number of PolSAR missions and openly available datasets due to open data policies, there is an ever increasing demand for processing tools. In addition, the recent migration of PolSAR datasets into cloud-based platforms (E.g. NASA-ESA multi-mission algorithm and analysis platform (MAAP) [@albinet2019joint]) requires processing tools that can be directly used on cloud-based platforms such as [JupyterHub on Kubernetes] (https://z2jh.jupyter.org/en/stable/) and Google colab (https://colab.research.google.com/). Currently, there are limited open-source tools available for PolSAR data processing ([PolSARpro](https://earth.esa.int/web/polsarpro/home), [SNAP](https://step.esa.int/main/toolboxes/snap/), and PolSAR-tools QGIS plug-in [@bhogapurapu2021polsar]). Although researchers widely use these tools, their direct use cases on a cloud native automatic processing environment are limited due to their primary architecture of a GUI-based approach. There are workarounds to implement these tools using Python wrappers and custom batch scripts. However, setting up automated workflows require significant programming skills, which limits the users. Therefore, the current work developed a Python package, `polsartools`, for PolSAR data processing, available through Python Package Index ([PyPI](https://pypi.org/)) and [Anaconda Package repository](https://anaconda.org/anaconda/repo). Further, to build upon the vision of [PolSARpro](https://earth.esa.int/web/polsarpro/home), the current package also implements several analysis functions that can be used for a teaching/demonstrating tutorials of polarimetric SAR applications through [Jupyter notebooks](https://github.com/Narayana-Rao/polsartools-tutorials). 

A typical polarimetric SAR data processing workflow contains the following steps: 
- Load data and extract polarimetric matrix (Scattering matrix, [$S$], Covariance matrix, [$\langle C \rangle$] or Coherence matrix [$\langle T \rangle$]).
- A second order [$\langle C \rangle$] or [$\langle T \rangle$] matrix is in general obtained through multi-looking (spatial averaging) single look complex (SLC) [$S$] elements. 
- Optional speckle filtering.
- Computing derived polarimetric  parameters (decomposition parameters or other descriptors).
- Analysis of the derived parameters and thier applications.

Based on the above processing steps, the functionalities of `polsartools` package can be broadly categorized into Processing & Analysis. Processing functions generate several SAR polarimetric parameters in a raster format, while analysis functions generate plots and quicklooks from the PolSAR data. \autoref{fig:flowchart} presents the core processing architecture of the `polsartools` package. PolSAR tools currently supports data from several SAR sensors, including spaceborne and airborne sensors. The package designed to support various forms of PolSAR datasets viz. full/quad-, dual-, compact- and hybrid polarimtry. 

![Schematic of core processing flow of polsartools package \label{fig:flowchart}](figures/flowchart.png){width=70%}

# Usage

## Installation
PolSAR tools package can be installed using prebuilt wheel (`.whl`) or Anaconda (`.conda`) packages. Alternatively it can be build and installed from source code in the GitHub repository. 

### Pre-built packages
The prebuilt packages are provided for `linux`, `windows`, and `OSX-mac silicon` operating systems for various Python versions (from 3.6 to 3.13). These pre-build packages can be installed through either of the following ways.

- Wheel 

      pip install polsartools

- Conda

      conda install polsartools -c bnarayanarao


### Build and install from source 

To build and install from source, one can use the [`pip`](https://pypi.org/project/pip/) package manager. To build the  `polsartools` from source windows users need to install as [Microsoft C++ build tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/). It can installed from source using the following command.

    pip install git+https://github.com/Narayana-Rao/polsartools.git#egg=polsartools

Any specific stable version (for e.g. v0.8) can also be built and installed using the following command.

    pip install git+https://github.com/Narayana-Rao/polsartools.git@0.8#egg=polsartools



<!-- ### Testing the installation -->

<!-- For the developers who want to contribute this opensource effort through their own code and test the integrated functions Extra packages are required. For documentation update sphinx, pydata, pytests -->

## I/O

The PolSAR data are in general distributed in either raw binary with metadata, HDF5 or tiff with metadata formats. The current implementaion of `polsartools` supports all these three input data formats which are sensor specific and extracts any of the polarimetric matrix elements $S2,~Sxy,~C4,~T4,~C3,~T3,~C2,~T2$. The output data  format is by default set to widely used GeoTiff format with optional arguments for Cloud Optimized Geotiff (COG). Additionally, it also supports raw binary format as output for better compatability with other softwares like PolSARpro.

For the analysis functions, the default output is the plot display and also returns axis handle for further manipulation. Additionally, aplot plot name argument `pname` can be provided to export the default output to a `.png` or any other matplotlib compatable output extension with `dpi=300`. 


## Tutorials and docs

PolSAR tools package structured in modular way for efficiency. Following is the list of main and sub modules in the package.

- `sensors`: handles importing the PolSAR data from different sensors and extract PolSAR matrix
- `preprocessing`: proprocess and prepare the data for Polarimetric functions
- `polsar`: consists core polarimetric functions with submodules for different PolSAR forms as listed below
  - `fp` full/quad polarimetry functions module
  - `cp` compact or hybrid polarimetry functions module
  - `dxp` dual cross polarimetry functions module
  - `dcp` dual co polarimetry functions module
  - `others` other miscillenous polarimetry functions module
- `analysis`: consists analysis functions
- `utils`: extra utilities  and helper functions 

Although the functions accessesd through each module (E.g.`polsartools.sensors.sensor_name.data_format(*args)`), for an easy usage, all function are given unique identfiers such that they can directly accessed from the parent package. For example, `polsartools.example_function(*args)`.

A comprehensive collection of Jupyter notebooks for each sensor specific data processing and other polarimetric functionalities are provided in the polsartools-tutorials git repository: [https://github.com/Narayana-Rao/polsartools-tutorials](https://github.com/Narayana-Rao/polsartools-tutorials).

A detailed documentaion of all the available functions is provided at: [polsartools.readthedocs.io](https://polsartools.readthedocs.io/en/latest/).



<!-- 
# Acknowledgements
The author would like to  -->

# References

