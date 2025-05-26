Dual-pol Radar Vegetation Index (``dprvi``)
============================================
.. This functionality computes the dual polarimetric radar vegetation index for dual polarimetric (HH | HV), (VV | VH) SAR data. The required input and the computed output are as follows:

.. .. code-block:: python

..         input : input_C2_folder, window_size
..         output: DpRVI.bin

.. autofunction:: polsartools.dprvi
   :noindex:

The formulation of DpRVI is as follows:

.. math::

    \text{DpRVI} = 1- \text{DOP}_{dp}\Big(\frac{\lambda_1}{\lambda_1+\lambda_2}\Big)

where,

.. math::
	\text{DOP}_{dp} = \sqrt{1-\frac{4\times \text{det ([C2])}}{\text{(Trace [C2])}^2}}

:math:`\text{[C2]}` is co-variance matrix,  and :math:`\lambda_1, \lambda_2` are the eigen values of :math:`\langle\mathbf{[C2]}\rangle` matrix in descending order. Further details on DpRVI can be obtained from [[5]](#5)


