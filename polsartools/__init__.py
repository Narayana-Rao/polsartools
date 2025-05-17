# polsartools/__init__.py

import warnings
warnings.filterwarnings("ignore")

__version__ = "0.6.1"  

# Import submodules
from . import polsar
from . import preprocessing
from . import utils 
from . import sensors
# Import functions from the submodules for direct access
from .polsar.fp import grvi, rvifp, mf3cf, mf4cf, dopfp, prvifp,nnedfp, neufp,halphafp,yam4cfp
from .polsar.cp import cprvi, dopcp, misomega, mf3cc
from .polsar.dxp import dprvi, dopdp, prvidp, rvidp, halphadp
from .polsar.dcp import mf3cd
from .sensors.uavsar import uavsar_grd,uavsar_mlc
from .sensors.nisar import nisar_gslc
from .utils import convert_T3_C3,convert_C3_T3
from .preprocessing.filters import boxcar,rlee

## CPP functions
# import polsartools.refined_lee
# import polsartools.moving_average
# import polsartools.sum_arrays
# import polsartools.cprvicpp #process_chunk_cprvicpp
__all__ = [
    'uavsar_grd', 'uavsar_mlc', 'nisar_gslc', #import data from sensors
    # 'rl', 'boxcar', #import filters
    'convert_T3_C3', 'convert_C3_T3', 'pauliRGB',#import utils
    'grvi', 'rvifp', 'mf3cf', 'mf4cf', 'dopfp', 'prvifp', 'neufp', 'nnedfp', 'halphafp', 'yam4cfp', # Full-pol
    'cprvi', 'dopcp', 'misomega', 'mf3cc',                 # Compact-pol
    'dprvi', 'dopdp', 'prvidp', 'rvidp', 'halphadp'        # Dual-cross-pol
    'mf3cd'                                         # Dual-co-pol
    
]