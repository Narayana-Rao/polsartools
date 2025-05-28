Compact-pol Radar Vegetation Index (``cprvi``)
==============================================
.. This functionality computes the compact-pol radar vegetation index for compact polarimetric SAR data. The required input and the computed output are as follows:

.. .. code-block:: python

..         input : input_C2_folder, window_size
..         output: CpRVI.bin


.. autofunction:: polsartools.cprvi
   :noindex:



The formulation of the CpRVI is as follows:

.. math::
    
    \text{CpRVI}=\left(1-\dfrac{3}{2}\text{GD}_{\text{ID}}\right)\Big(\frac{p}{q}\Big)^{2(\frac{3}{2}\text{GD}_{\text{ID}})}\\
    p=\text{min\{SC,OC\}},q=\text{max\{SC,OC\}}\\
    \text{SC}=\frac{S_0-S_3}{2};\qquad{}\text{OC}=\frac{S_0+S_3}{2};\\
    S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\
    S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)}\\


where, :math:`\text{GD}_\text{ID}` is the geodesic distance between Kennaugh matrices (:math:`\mathbf{K}`) of the observed and the ideal depolarizer, :math:`p, q` are minimum and maximum values of :math:`\text{SC}` and :math:`\text{OC}` which are functions of stocks parameters (:math:`S_0`, :math:`S_1`, :math:`S_2`, and :math:`S_3`). A detailed explanation of CpRVI is available in [[6]](#6).