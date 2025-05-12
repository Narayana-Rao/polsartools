
Example Usage
==============

More sample use cases are provided in the `examples`_ folder. Detailed documentation is available at `polsartools.readthedocs.io <https://polsartools.readthedocs.io/en/latest/>`_ 

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

