Compact-pol data format
=======================

Compact-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). 
A typical file structure of C2 matrix is as follows:

+-----------------------------+
|       C2 matrix files       |
+==============+==============+
| C11.bin      | C11.hdr      |
+--------------+--------------+
| C12_real.bin | C12_real.hdr |
+--------------+--------------+
| C12_imag.bin | C12_imag.hdr |
+--------------+--------------+
| C22.bin      | C22.hdr      |
+--------------+--------------+
