# PolSARtools PyPI package (beta version!)
 
 [![Documentation Status](https://readthedocs.org/projects/polsartools/badge/?version=latest)](https://polsartools.readthedocs.io/en/latest/?badge=latest)

#### Installation
```
pip install polsartools
```
 
#### Prerequesites 
 
 ```gdal, Numpy```
 
#### gdal installation error fix
 
 ```conda install -c conda-forge gdal```

#### Example
```python
import polsartools as pst

T3_folder = r'../T3'
windows_size=3

ps,pd,pv,pc,tfp,taufp = pst.mf4cf(T3_folder,window_size=window_size,write_flag=1)

```

#### sample use case is provided in [tests](https://github.com/Narayana-Rao/polsartools/tests)
