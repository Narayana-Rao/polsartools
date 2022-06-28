General Information
===================

This plugin generates derived SAR parameters (viz. vegetation indices, polarimetric decomposition parameters) from input polarimetric matrix (C3, T3, C2, T2). The input data needs to be in [PolSARpro](https://earth.esa.int/web/polsarpro/home)/[ENVI](https://www.l3harrisgeospatial.com/Software-Technology/ENVI) format (\*.bin and \*.hdr). It requires [numpy](https://numpy.org/), [matplotlib](https://matplotlib.org/) python libraries pre-installed.

Installation
-------------------
> **__Note:__** SAR tools requires QGIS version >=3.0.

* The easiest way (requires internet connection) : 
	- Open QGIS -> Plugins -> Manage and Install Plugins... -> select ```All``` tab -> search for ```SAR tools``` --> select and install plugin
* Alternative way (offline installation) : 
	- Go to [releases](https://github.com/Narayana-Rao/SAR-tools/releases) of this repository -> select desired version -> download the ```.zip``` file.
	- Open QGIS -> Plugins -> Manage and Install Plugins... -> ```install from ZIP``` tab --> select the downloaded zip --> install plugin (ignore warnings, if any).
 

Available functionalities:
-------------------
  * Indices :
  	* Radar Vegetation Index (RVI) (Full-pol and dual-pol)
  	* Generalized volume Radar Vegetation Index (GRVI)
  	* Polarimetric Radar Vegetation Index (PRVI) (Full-pol and dual-pol) 
  	* Dual-pol Radar Vegetation Index (DpRVI)
  	* Degree of Polarization (DOP) (Full-pol, dual-pol, and compact-pol)
  	* Compact-pol Radar Vegetation Index (CpRVI)
  
  * Polarimetric Decompositions : 
  	* Model free 3-Component decomposition for full-pol data (MF3CF).
  	* Model free 3-Component decomposition for compact-pol data (MF3CC) 
  	* Improved S-Omega decomposition for compact-pol data (iS-Omega)