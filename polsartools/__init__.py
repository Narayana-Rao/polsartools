"""
A package for processing Polarimetric Synthetic Aperture Radar (PolSAR) data.
"""
# polsartools/__init__.py

import warnings
warnings.filterwarnings("ignore")

__version__ = '0.5'  

# Import submodules
from . import polsar
from . import preprocessing
from . import utils 
from . import sensors
# Import functions from the submodules for direct access
from .polsar.fp import grvi, rvifp, mf3cf, mf4cf, dopfp, prvifp
from .polsar.cp import cprvi, dopcp, misomega, mf3cc
from .polsar.dxp import dprvi, dopdp, prvidp, rvidp, halphadp
from .sensors.uavsar import uavsar_grd,uavsar_mlc
__all__ = [
    'uavsar_grd', 'uavsar_mlc', #import from different sensors
    'grvi', 'rvifp', 'mf3cf', 'mf4cf', 'dopfp', 'prvifp', 'neufp', 'nnedfp', # Full-pol
    'cprvi', 'dopcp', 'misomega', 'mf3cc',                 # Compact-pol
    'dprvi', 'dopdp', 'prvidp', 'rvidp', 'halphadp'        # Dual-cross-pol
]