# polsartools/__init__.py

import warnings
warnings.filterwarnings("ignore")

__version__ = "0.7"  

# Import submodules
# from . import polsar
# from . import preprocessing
# from . import utils 
# from . import sensors
# Import functions from the submodules for direct access
from .polsar.fp import grvi, rvifp, mf3cf, mf4cf, dopfp, prvifp,nnedfp, neufp,halphafp,yam4cfp,shannon_h_fp,halpha_plot_fp
from .polsar.cp import cprvi, dopcp, misomega, mf3cc
from .polsar.dxp import dprvi, dopdp, prvidp, rvidp, halphadp, shannon_h_dp,halpha_plot_dp
from .polsar.dcp import mf3cd
from .sensors.uavsar import uavsar_grd,uavsar_mlc
from .sensors.nisar import nisar_gslc,nisar_rslc
from .sensors.alos2 import alos2_fbd_l11
from .sensors.chyaan2 import chyaan2_fp
from .sensors.rs2_fp import rs2_fp
from .utils import convert_T3_C3,convert_C3_T3
from .preprocessing.filters import boxcar,rlee
# from .utils.pauliRGB import read_bin
## CPP functions
# import polsartools.refined_lee
# import polsartools.moving_average
# import polsartools.sum_arrays
# import polsartools.cprvicpp #process_chunk_cprvicpp
__all__ = [
    'uavsar_grd', 'uavsar_mlc', 'nisar_gslc', 'nisar_rslc', 'alos2_fbd_l11', 'chyaan2_fp','rs2_fp', #import data from sensors
    'rlee', 'boxcar', #import filters
    'convert_T3_C3', 'convert_C3_T3', 'pauliRGB',#import utils
    'grvi', 'rvifp', 'mf3cf', 'mf4cf', 'dopfp', 'prvifp', 'neufp', 'nnedfp', 'halphafp', 'shannon_h_fp','yam4cfp', 'halpha_plot_fp', # Full-pol
    'cprvi', 'dopcp', 'misomega', 'mf3cc',                 # Compact-pol
    'dprvi', 'dopdp', 'prvidp', 'rvidp', 'halphadp', 'shannon_h_dp', 'halpha_plot_dp',      # Dual-cross-pol
    'mf3cd'                                         # Dual-co-pol
    
]