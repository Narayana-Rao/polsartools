---
title: 'PolSAR tools: A python package to process and analyse Polarimetric Synthetic Aperture Radar (PolSAR) data'
tags:
  - SAR
  - PolSAR
  - Speckle filtering
  - polarimetric signature
  - Polarimetric decompositions
  - Vegetation indices


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

# Introduction
Synthetic Aperture Radar has been proven to be a highly reliable Earth observation (EO) remote sensing technique due to its all-weather and day-night capabilities. It has evolved in spatial resolution, observing frequencies, and applications since its initial development and use on an airborne platform seven decades ago. It took nearly four decades since then to have a full/quad polarimetric SAR (PolSAR) instrument in space for earth observation [@ulaby1981microwave;@jordan2002sir]. Since then, PolSAR has become a critical data source in various earth observation applications [@lee2017polarimetric;@cloude2010polarisation;@van2011synthetic]. The widespread use of SAR data applications has become possible through several airborne and spaceborne SAR instruments and their open data. As the world moves towards open data and science, several operational SAR satellite missions (e.g., Sentinel-1, ALOS-2/4, EOS-04, BIOMASS, etc) have been collecting and distributing data openly. Further, the proposed future missions, such as NISAR [@rosen2017global], ROSE-L [@davidson2021rose], and the Sentinel-1 Next Generation (S1 NG) [@geudtner2021copernicus], are planned to collect Petabytes of EO data. For an efficient use of these PolSAR data, it is essential to process them and derive interpretable polarimetric parameters to support further downstream applications. 


# Statement of need

With the increasing number of PolSAR missions and openly available datasets due to open data policies, there is an ever-increasing demand for processing tools. In addition, the recent migration of PolSAR datasets into cloud-based platforms (E.g., NASA-ESA multi-mission algorithm and analysis platform (MAAP) [@albinet2019joint]) requires processing tools that can be directly used on cloud-based platforms such as [JupyterHub on Kubernetes](https://z2jh.jupyter.org/en/stable/), [Google Colab](https://colab.research.google.com/), [AWS SageMaker](https://aws.amazon.com/sagemaker/), and [Azure ML](https://azure.microsoft.com/en-us/products/machine-learning/). Currently, there are limited open-source tools available for PolSAR data processing ([PolSARpro](https://earth.esa.int/web/polsarpro/home), [SNAP](https://step.esa.int/main/toolboxes/snap/), and PolSAR-tools QGIS plug-in [@bhogapurapu2021polsar]). Although researchers widely use these tools, their direct use cases on a cloud native automatic processing environment are limited due to their primary architecture of a GUI-based approach. There are workarounds to implement these tools using Python wrappers and custom batch scripts. However, setting up automated workflows requires significant programming skills, limiting users. Therefore, the current work developed a Python package, `polsartools`, for PolSAR data processing, available through Python Package Index ([PyPI](https://pypi.org/)) and [Anaconda Package repository](https://anaconda.org/anaconda/repo). Further, to build upon the vision of [PolSARpro](https://earth.esa.int/web/polsarpro/home), the current package also implements several analysis functions that can be used for teaching/demonstrating tutorials of polarimetric SAR applications through [Jupyter notebooks](https://github.com/Narayana-Rao/polsartools-tutorials). 

A typical polarimetric SAR data processing workflow contains the following steps: 
- Load data and extract polarimetric matrix (Scattering matrix, [$S$], Covariance matrix, [$\langle C \rangle$] or Coherence matrix [$\langle T \rangle$]).
- A second order [$\langle C \rangle$] or [$\langle T \rangle$] matrix is in general obtained through multi-looking (spatial averaging) single look complex (SLC) [$S$] elements. 
- Optional speckle filtering.
- Computing derived polarimetric parameters (decomposition parameters or other descriptors).
- Analysis of the derived parameters and their applications.

Based on the above processing steps, the functionalities of the `polsartools` package can be broadly categorized into Processing & Analysis. Processing functions generate the SAR polarimetric parameters in a raster format, while analysis functions generate plots and quicklooks from the PolSAR data. \autoref{fig:flowchart} presents the core processing architecture of the `polsartools` package. The `polsartools` package currently supports data from several SAR sensors, including spaceborne and airborne sensors. The package is also designed to support PolSAR datasets from different polarimetric modes, viz., full/quad-, dual-, compact-, and hybrid polarimetry. 

![Schematic of core processing flow of polsartools package \label{fig:flowchart}](figures/flowchart.pdf){width=70%}

# Usage

## Installation
PolSAR tools package can be installed using prebuilt [Wheel](https://packaging.python.org/en/latest/glossary/#term-Wheel) (`.whl`) or [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/concepts/packages.html#conda-file-format) (`.conda`) packages. Alternatively, it can be built and installed from source code in the GitHub repository ([polsartools](https://github.com/Narayana-Rao/polsartools)). 

### Prebuilt packages
The prebuilt packages are provided for `linux`, `windows`, and `OSX-mac silicon` operating systems for various Python versions (from 3.6 to 3.13). Users can install these prebuilt packages in either of the following ways.

- Wheel 

       pip install polsartools

- Conda

      conda install polsartools -c bnarayanarao


### Build and install from source 

One can use the [`pip`](https://pypi.org/project/pip/) package manager to build and install from source. To build the  `polsartools` from source, Windows users need to install [Microsoft C++ build tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/). It can be installed from source using the following command.

    pip install git+https://github.com/Narayana-Rao/polsartools.git#egg=polsartools

Any specific stable version (e.g., v0.8) can also be built and installed using the following command.

    pip install git+https://github.com/Narayana-Rao/polsartools.git@0.8#egg=polsartools



<!-- ### Testing the installation -->

<!-- For the developers who want to contribute to this open-source effort through their own code and test the integrated functions, Extra packages are required. For documentation update sphinx, pydata, pytests -->

## I/O

The PolSAR data are generally distributed in raw binary with metadata, HDF5, or TIFF with metadata formats. Current implementation of `polsartools` supports all three input data formats, which are sensor-specific, and extracts any of the polarimetric matrix elements $S2,~Sxy,~C4,~T4,~C3,~T3,~C2,~T2$. By default, the output data format is set to the widely used GeoTiff format with optional arguments for Cloud Optimized Geotiff (COG). It also supports raw binary format as output for better compatibility with other software like PolSARpro.

For the analysis functions, the default output is to display a plot, and returns the axis handle for further manipulation. Additionally, a plot name argument `pname` can be provided to export the default output to a `.png` or any other [matplotlib](https://matplotlib.org/)-compatible output extension. 


## Tutorials and docs

The `polsartools` package is structured in a modular way for efficiency. The main and submodules of the package are listed below.

- `sensors`: handles importing the PolSAR data from different sensors and extracting the PolSAR matrices
- `preprocess`: includes preprocessing and preparing the data for Polarimetric functions
- `polsar`: consists of core polarimetric functions with submodules for different PolSAR forms as listed below
  - `fp` full/quad polarimetry functions module
  - `cp` compact or hybrid polarimetry functions module
  - `dxp` dual cross-polarimetry functions module
  - `dcp` dual co-polarimetry functions module
  - `others` other miscellaneous polarimetry functions module
- `analysis`: consists of analysis functions
- `utils`: extra utilities and helper functions 

Although the functions can be accessed through each module (E.g., `polsartools.sensors.sensor_name.data_format(*args)`), for easy usage, all functions are given unique identifiers such that these functions can be directly accessed from the parent package. For example, `polsartools.example_function(*args)`.

A detailed documentation of all the functions is available at: [polsartools.readthedocs.io](https://polsartools.readthedocs.io/en/latest/). Additionally, a collection of Jupyter notebooks for each sensor-specific data processing and other polarimetric functionalities is provided in the polsartools-tutorials git repository: [https://github.com/Narayana-Rao/polsartools-tutorials](https://github.com/Narayana-Rao/polsartools-tutorials).

<!-- 
# Acknowledgements
The author would like to  -->

# References

