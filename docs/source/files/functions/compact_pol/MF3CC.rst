Model Free 3-Component decomposition (``mf3cc``)
=================================================
..  This functionality computes the model free 3 component scattering power decomposition for compact polarimetric SAR data. The required input and the computed output are as follows:

.. .. code-block:: python

..     input : input_C2_folder, window_size, tau
..     output: Ps_CP.bin, Pd_CP.bin, Pv_CP.bin, Theta_CP.bin


.. autofunction:: polsartools.mf3cc
   :noindex:



The formulation of the scattering powers (:math:`P_s` : Surface, :math:`P_d`: Double bounce, :math:`P_v`: volume) is as follows:

.. math::

    P_{d}^{\text{CP}}=\frac{m_{\text{FP}}{S_0}}{2}{\left(1-\sin2\theta_{\text{CP}}\right)};\\P_{v}^{\text{CP}}={S_0}\left(1-m_{\text{CP}}\right);\\P_{s}^{\text{CP}}=\frac{m_{\text{CP}}{S_0}}{2}\left(1+\sin2\theta_{\text{CP}}\right)

where :math:`m_\text{CP}` is degree of polarization; :math:`\theta_\text{CP}` : scattering type parameter; :math:`S_0, S_3`, are Stokes parameters. The derivation of these parameters in-terms of covariance matrix (C2) elements is as shown below. Further details can be obtained from [[4]](#4)

.. math::

    m_{\text{CP}}=\sqrt{1-\frac{4|\mathbf{C2}|}{\big(\mathrm{Trace}(\mathbf{C2})\big)^2}};\qquad{}\tan\theta_{\text{CP}}=\frac{m_{\text{CP}}{S_0}\left(\text{OC}-\text{SC}\right)}{\text{OC}\times\text{SC}+m_{\text{CP}}^{2}{S_0}^{2}}\\
    S_0=\text{C11+C22};\qquad{}\qquad{}\qquad{}\qquad{}\qquad{}S_1=\text{C11-C22};\\
    S_2=\text{C12+C21};\qquad{}\qquad{}\qquad{}\qquad{}S_3=\pm\text{j(C12-C21)};\\
    \text{SC}=\frac{S_0-S_3}{2};\qquad{}\qquad{}\qquad{}\qquad{}\qquad{}\text{OC}=\frac{S_0+S_3}{2};