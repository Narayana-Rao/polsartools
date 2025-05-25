Full-pol data format
====================

Full-pol functionalities require the SAR data in the form of scattering (S2), covariance (C3) or coherency matrix (T3). 
A typical file structures of T3 and C3 matrices are as follows: 

+-----------------------------+-----------------------------+
|       C3 matrix files       |       T3 matrix files       |
+==============+==============+==============+==============+
| C11.bin      | C11.hdr      | T11.bin      | T11.hdr      |
+--------------+--------------+--------------+--------------+
| C12_real.bin | C12_real.hdr | T12_real.bin | T12_real.hdr |
+--------------+--------------+--------------+--------------+
| C12_imag.bin | C12_imag.hdr | T12_imag.bin | T12_imag.hdr |
+--------------+--------------+--------------+--------------+
| C13_real.bin | C13_real.hdr | T13_real.bin | T13_real.hdr |
+--------------+--------------+--------------+--------------+
| C13_imag.bin | C13_imag.hdr | T13_imag.bin | T13_imag.hdr |
+--------------+--------------+--------------+--------------+
| C22.bin      | C22.hdr      | T22.bin      | T22.hdr      |
+--------------+--------------+--------------+--------------+
| C23_real.bin | C23_real.hdr | T23_real.bin | T23_real.hdr |
+--------------+--------------+--------------+--------------+
| C23_imag.bin | C23_imag.hdr | T23_imag.bin | T23_imag.hdr |
+--------------+--------------+--------------+--------------+
| C33.bin      | C33.hdr      | T33.bin      | T33.hdr      |
+--------------+--------------+--------------+--------------+