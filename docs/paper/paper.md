---
title: 'PolSAR tools: A python package to process polarimetric SAR (PolSAR) data'
tags:
  - SAR
  - PolSAR
  - pip
  - conda
  - Speckle filtering
  - polarimetric signature
  - Vegetation indices
  - Polarimetric decompositions

authors:
  - name: Narayanarao Bhogapurapu
    orcid: 0000-0002-6496-7283
    affiliation: 1 
    affiliation: 2 
affiliations:
 - name: Microwave Remote Sensing Laboratory, University of Massachusetts Amherst, USA
   index: 1
 - name: Division of Geological and Planetary Sciences, California Institute of Technology, USA
   index: 2

date: 4 July 2025
bibliography: refs.bib

---

# Statement of need

The demand for processing tools increases with the increasing number of ***Synthetic Aperture Radar (SAR)*** satellite missions and datasets. However, to process SAR data, a minimal number of free tools are available ([PolSARpro](https://earth.esa.int/web/polsarpro/home), [SNAP](https://step.esa.int/main/toolboxes/snap/)) that consolidate all necessary pre-processing steps. Bearing this in mind, there is a need to develop specific tools for the remote sensing user community to derive polarimetric descriptors like vegetation indices and decomposition parameters. With current 


 input polsar data -> functions -> (two type of functions output a raster/s (of requsted parameter) or a plot/image of png (for quicklooks, polarimetric plots for analysis, like h-alpha plot, polarimetric signatures etc)); now in core function given SAR dataset will be read through a chunk by chunk manner in parllel wiht cpu-1 workers all chunks processed parllel and temporaly all output chunks are wrtten to geotiffs and then mosaiced together to get required output parameters. 

<!-- <p align="center">
  <img src="figures/flowchart.pdf" alt=""/>
</p> -->

![core processing flow](figures/flowchart.png)

# References