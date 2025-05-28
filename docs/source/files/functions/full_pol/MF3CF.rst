Model Free 3-Component decomposition (``mf3cf``)
================================================
.. This functionality computes the model free 3 component scattering power decomposition for full polarimetric SAR data. The required input and the computed output are as follows:

.. .. code-block:: python
        
..         input : input_T3/C3_folder, window_size
..         output: Ps_FP.bin, Pd_FP.bin, Pv_FP.bin, Theta_FP.bin


.. autofunction:: polsartools.mf3cf
   :noindex:


The formulation of the scattering powers (:math:`P_s` : Surface, :math:`P_d` : Double bounce, :math:`P_v` : volume) is as follows:

.. math::
    
    P_{d}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}{\left(1-\sin2\theta_{\text{FP}}\right)}\\P_{v}^{\text{FP}}={\text{Span}}\left(1-m_{\text{FP}}\right)\\P_{s}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}\left(1+\sin2\theta_{\text{FP}}\right)

where :math:`m_\text{FP}` is degree of polarization, :math:`\theta_\text{FP}` scattering type parameter, Span is the sum of the diagonal elements os coherence matrix (T3).  The derivation of these parameters in-terms of coherancey matrix (T3) elements is as shown below. Further details can be obtained from [[4]](#4)

.. math::

    m_{\text{FP}}=\sqrt{1-\frac{27|\mathbf{T3}|}{\big(\mathrm{Trace}(\mathbf{T3})\big)^3}};\qquad{}\tan\theta_{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}\left(T_{11}-T_{22}-T_{33}\right)}{T_{11}\left(T_{22}+T_{33}\right)+m_{\text{FP}}^{2}{\text{Span}}^{2}}
    
    \text{Span}=T_{11}+T_{22}+T_{33}
    
