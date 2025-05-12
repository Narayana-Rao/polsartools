General Information
===================

This package generates derived SAR parameters  (polarimetric decomposition parameters and descriptors)  from input polarimetric matrix (C3, T3, C2, T2). The input data needs to be in `PolSARpro`_/`ENVI`_ format (\*.bin and \*.hdr). It requires `gdal`_, `numpy`_, `scipy`_, `matplotlib`_ python libraries pre-installed.

Installation
-------------------

To install **polsartools**, use the following command:

.. code-block:: bash

    pip install polsartools

.. note::

    polsartools tools requires python >=3.6.


Prerequisites
---------------
Ensure that the required dependencies are installed:

.. code-block:: bash

    conda install gdal simplekml numpy scipy 

Example Usage
--------------

More sample use cases are provided in the `examples`_ folder.

.. _examples: https://github.com/Narayana-Rao/polsartools/tree/main/examples

.. code-block:: python

    import polsartools as pst

    def main():
        # Define the file paths for compact-pol covariance matrix 
        compact_c2 = r'./sample_data/compact_pol/C2_RHV'
        chi_in = 45  
        window_size = 3  # Window size for processing

        """ Decompositions for Compact Polarimetric Data """
        pst.mf3cc(compact_c2, chi_in=chi_in, window_size=window_size)
        pst.misomega(compact_c2, chi_in=chi_in, psi_in=0, window_size=window_size)

        """ Descriptors for Compact Polarimetric Data """
        pst.cprvi(compact_c2, chi_in=chi_in, window_size=window_size)
        pst.dopcp(compact_c2, chi_in=chi_in, window_size=window_size)

        full_T3 = r'./sample_data/full_pol/T3'
        window_size = 3  

        """ Decompositions for Full Polarimetric Data """
        pst.halphafp(full_T3, window_size=window_size)
        pst.neufp(full_T3, window_size=window_size)
        pst.nnedfp(full_T3, window_size=window_size)
        pst.mf3cf(full_T3, window_size=window_size)
        pst.mf4cf(full_T3, window_size=window_size)

        """ Descriptors for Full Polarimetric Data """
        pst.dopfp(full_T3, window_size=window_size)
        pst.grvi(full_T3, window_size=window_size)
        pst.rvifp(full_T3, window_size=window_size)
        pst.prvifp(full_T3, window_size=window_size)

        dxp_C2 = r'./sample_data/dual_pol/C2_VVVH'
        window_size = 3  

        """ Descriptors for Dual-Polarimetric Data """
        pst.dprvi(dxp_C2, window_size=window_size)
        pst.rvidp(dxp_C2, window_size=window_size)
        pst.prvidp(dxp_C2, window_size=window_size)
        pst.dopdp(dxp_C2, window_size=window_size)

        """ Refined Lee Polarimetric Speckle Filter """
        T3_folder = r'./sample_data/full_pol/T3'
        pst.rlee(T3_folder, window_size=5)
        pst.utils.pauliRGB(T3_folder)

    if __name__ == "__main__":
        main()


Available functionalities
--------------------------
  1. **Full-pol** 

    * H-Alpha decomposition 
    * Non-negative Eigen value decomposition
    * Neumann Decomposition
    * Model free 3-Component decomposition for full-pol data (`MF3CF <functions/full_pol/MF3CF.html>`_) 
    * Radar Vegetation Index (`RVIfp <functions/full_pol/RVI_fp.html>`_)
    * Generalized volume Radar Vegetation Index (`GRVI <functions/full_pol/GRVI.html>`_)
    * Polarimetric Radar Vegetation Index (`PRVIfp <functions/full_pol/PRVI_fp.html>`_)
    * Degree of Polarization (`DOPfp <functions/full_pol/DOP_fp.html>`_)

  2. **Compact-pol**

    * Model free 3-Component decomposition for compact-pol data (MF3CC) 
    * Improved S-Omega decomposition for compact-pol data (iS-Omega) 
    * Compact-pol Radar Vegetation Index (CpRVI) 
    * Degree of Polarization (DOP)  

  3. **Dual-pol**

    * Dual-pol Radar Vegetation Index (`DpRVI <functions/dual_pol/DpRVI.html>`_) 
    * Radar Vegetation Index (`RVIdp <functions/dual_pol/RVI_dp.html>`_) 
    * Degree of Polarization (`DOPdp <functions/dual_pol/DOP_dp.html>`_)
    * Polarimetric Radar Vegetation Index (`PRVIdp <functions/dual_pol/PRVI_dp.html>`_) 


Functions description
---------------------

Description and the details of all the core functions of this package are available here: (`Functions description <functions_description.html>`_)


Contributions
-------------

1) **Contribute to the software**:
   See the `Contribution guidelines`_ for this project.

.. _Contribution guidelines: help/CONTRIBUTING.md

2) **Report issues**:
   Please raise your issues here: `<https://github.com/Narayana-Rao/polsartools/issues>`_



.. _PolSARpro: https://earth.esa.int/web/polsarpro/home
.. _ENVI: https://www.l3harrisgeospatial.com/Software-Technology/ENVI
.. _gdal: https://gdal.org/en/latest/
.. _scipy: https://scipy.org/
.. _numpy: https://numpy.org/
.. _matplotlib: https://matplotlib.org/
.. _releases: https://github.com/Narayana-Rao/SAR-tools/releases