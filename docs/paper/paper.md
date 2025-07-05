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
    affiliation: "1, 2" 

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


The functionalities can be boradly catogerized into Processing & Analysis. Processing functions genrate several SAR polarimetric parameters in a raster format while analysis functions genrates plots and quick looks from the PolSAR data. A schematic processing flow is presented in Figure \autoref{fig:flowchart}. 

![Schematic of core processing flow of polsartools package \label{fig:flowchart}](figures/flowchart.png){width=60%}



# References