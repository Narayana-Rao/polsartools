Generalized volume based Radar Vegetation Index (``grvi``)
===========================================================
.. This functionality computes the generalized volume based radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:

.. .. code-block:: python

..     input : input_T3/C3_folder, window_size
..     output : GRVI.bin

.. autofunction:: polsartools.grvi
   :noindex:


The formulation of GRVI is as follows:

.. math::

    \text{GRVI} = \left(1 - \text{GD}_{\text{GV}}\right)\Big(\frac{p}{q}\Big)^{2\,\text{GD}_{\text{GV}}}, \quad 0\le \text{GRVI} \le 1

where, :math:`\text{GD}_{\text{GV}}` is the geodesic distance between Kennaugh :math:`(\mathbf{K})` matrices of the observed and the generalized volume scattering model, :math:`p,q` are minimum and maximum value of distances between :math:`\mathbf{K}` matrices of the observed and elementary targets respectively. A detailed explanation of GRVI is available in.




