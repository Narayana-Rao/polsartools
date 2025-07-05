Dual-pol data format
=====================

C2 – Dual-pol Covariance Matrix
---------------------------------

Dual-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). 
A typical file structures of C2 matrix is as follows:

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

T2 – Dual-co-pol Coherency Matrix
---------------------------------

Dual-co-pol functionalities require the SAR data in the form of 2x2 coherency matrix (T2). 
A typical file structures of T2 matrix is as follows:

+----------------+--------------------------+------------------+
| Element        | .bin/.hdr Files          | .tif Files       |
+================+==========================+==================+
| T11            | T11.bin / T11.hdr        | T11.tif          |
+----------------+--------------------------+------------------+
| T12 (Real)     | T12_real.bin / .hdr      | T12_real.tif     |
+----------------+--------------------------+------------------+
| T12 (Imag)     | T12_imag.bin / .hdr      | T12_imag.tif     |
+----------------+--------------------------+------------------+
| T22            | T22.bin / T22.hdr        | T22.tif          |
+----------------+--------------------------+------------------+