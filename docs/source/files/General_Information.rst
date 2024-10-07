General Information
===================

This package generates derived SAR parameters (viz. vegetation indices, polarimetric decomposition parameters) from input polarimetric matrix (C3, T3, C2, T2). The input data needs to be in `PolSARpro`_/`ENVI`_ format (\*.bin and \*.hdr). It requires `gdal`_, `numpy`_, `scipy`_, `matplotlib`_ python libraries pre-installed.

Installation
-------------------

.. note::

    polsartools tools requires python >=3.6.

pip install polsartools
 

Available functionalities
--------------------------
  1. **Full-pol** 

    * Model free 3-Component decomposition for full-pol data (MF3CF) 
    * Radar Vegetation Index (RVI) 
    * Generalized volume Radar Vegetation Index (GRVI) 
    * Polarimetric Radar Vegetation Index (PRVI) 
    * Degree of Polarization (DOP) 

  2. **Compact-pol**

    * Model free 3-Component decomposition for compact-pol data (MF3CC) 
    * Improved S-Omega decomposition for compact-pol data (iS-Omega) 
    * Compact-pol Radar Vegetation Index (CpRVI) 
    * Degree of Polarization (DOP)  

  3. **Dual-pol**

    * Dual-pol Radar Vegetation Index (`DpRVI <functions/dual_pol/DpRVI.html>`_) 
    * Radar Vegetation Index (`RVI <functions/dual_pol/RVI_dp.html>`_) 
    * Degree of Polarization (`DOP <functions/dual_pol/DOP_dp.html>`_)
    * Polarimetric Radar Vegetation Index (`PRVI <functions/dual_pol/PRVI_dp.html>`_) 


Functions description
---------------------

Description and the details of all the core functions of this plugin are available here: (`Functions description <functions_description.html>`_)

Contributions
-------------

1) Contribute to the software

    `Contribution guidelines for this project  <https://github.com/Narayana-Rao/polsartools/tree/main/docs/CONTRIBUTING.md>`_


2) Report issues or problems with the software
  
  Please raise your issues here : https://github.com/Narayana-Rao/polsartools/issues

3) Seek support

  Please write to us: bnarayanarao@iitb.ac.in


.. _PolSARpro: https://earth.esa.int/web/polsarpro/home
.. _ENVI: https://www.l3harrisgeospatial.com/Software-Technology/ENVI
.. _gdal: https://gdal.org/en/latest/
.. _scipy: https://scipy.org/
.. _numpy: https://numpy.org/
.. _matplotlib: https://matplotlib.org/
.. _releases: https://github.com/Narayana-Rao/SAR-tools/releases