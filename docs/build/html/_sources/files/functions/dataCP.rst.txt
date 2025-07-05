Compact-pol data format
=======================

Compact-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). 
A typical file structure of C2 matrix is as follows:

+----------------+--------------------------+------------------+
| Element        | .bin/.hdr Files          | .tif Files       |
+================+==========================+==================+
| C11            | C11.bin / C11.hdr        | C11.tif          |
+----------------+--------------------------+------------------+
| C12 (Real)     | C12_real.bin / .hdr      | C12_real.tif     |
+----------------+--------------------------+------------------+
| C12 (Imag)     | C12_imag.bin / .hdr      | C12_imag.tif     |
+----------------+--------------------------+------------------+
| C22            | C22.bin / C22.hdr        | C22.tif          |
+----------------+--------------------------+------------------+