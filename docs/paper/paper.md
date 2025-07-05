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


# Summary
Synthetic Aperture Radar has been proven a more reliable as an earth obsvarion remotesening technique due to it allweather and day-night capabilties. Since its inital development and use on an airbone platforms seven decades ago it has significantly evolved interms of spatil resolution and observing frequencies. It took nearly four decades since then to have a full/quad polarimetric SAR (PolSAR) instrument in space for earth observation[`@ulaby1981microwave;@jordan2002sir`]. Since then PolSAR has become a crtical data source in a wide variety of earth obsertvation applications [`@lee2017polarimetric;@cloude2010polarisation;@van2011synthetic`]. On the other hand there are several operational (eg. Sentinel-1 ALOS-2 PALSAR-2, EOS-04, SAOCOM, BIOMASS, etc) and proposed future (NISAR, ROSE-L etc) SAR satellite missions have been collecting/planned to collect Petabytes of data and distributing openely. For an efficeint use of PolSAR data processing it to derive interprtable polarimetric paramters is essential for any further downstream applications. 


# Statement of need

<!-- The demand for processing tools increases with the increasing number of ***Synthetic Aperture Radar (SAR)*** satellite missions and datasets. However, to process SAR data, a minimal number of free tools are available ([PolSARpro](https://earth.esa.int/web/polsarpro/home), [SNAP](https://step.esa.int/main/toolboxes/snap/)) that consolidate all necessary pre-processing steps. Bearing this in mind, there is a need to develop specific tools for the remote sensing user community to derive polarimetric descriptors like vegetation indices and decomposition parameters. With current  -->
PolSAR tools QGIS plugin `[@bhogapurapu2021polsar]`

The functionalities can be boradly catogerized into Processing & Analysis. Processing functions genrate several SAR polarimetric parameters in a raster format while analysis functions genrates plots and quick looks from the PolSAR data. A schematic processing flow is presented in Figure \autoref{fig:flowchart}. 

![Schematic of core processing flow of polsartools package \label{fig:flowchart}](figures/flowchart.png){width=60%}

both pip conda packages

<!-- 
# Acknowledgements
The author would like to  -->

# References

