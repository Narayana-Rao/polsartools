Dual-pol data format
=====================

Dual-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). 
A typical file structures of C2 matrix is as follows:

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


Dual-co-pol functionalities require the SAR data in the form of 2x2 coherency matrix (T2). 
A typical file structures of T2 matrix is as follows:

+-----------------------------+
|       T2 matrix files       |
+==============+==============+
| T11.bin      | T11.hdr      |
+--------------+--------------+
| T12_real.bin | T12_real.hdr |
+--------------+--------------+
| T12_imag.bin | T12_imag.hdr |
+--------------+--------------+
| T22.bin      | T22.hdr      |
+--------------+--------------+