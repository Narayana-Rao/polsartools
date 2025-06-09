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
from .polsar.dxp import dprvi, dopdp, prvidp, rvidp, halphadp, shannon_h_dp,halpha_plot_dp,dprvic, dp_desc
from .polsar.dcp import mf3cd
from .sensors.uavsar import uavsar_grd,uavsar_mlc
from .sensors.nisar import nisar_gslc,nisar_rslc
from .sensors.alos2 import alos2_fbd_l11,alos2_hbq_l11
from .sensors.chyaan2 import chyaan2_fp
from .sensors.rs2_fp import rs2_fp
from .sensors.isro_asar import isro_asar
from .utils import convert_T3_C3,convert_C3_T3, convert_S2_CT
from .utils import fp_sign
from .preprocessing.filters import boxcar,rlee
# from .utils.pauliRGB import read_bin
## CPP functions
# import polsartools.refined_lee
# import polsartools.moving_average
# import polsartools.sum_arrays
# import polsartools.cprvicpp #process_chunk_cprvicpp
__all__ = [
    # SENSORS
    'uavsar_grd', 'uavsar_mlc','isro_asar',  
    'nisar_gslc', 'nisar_rslc',
    'alos2_fbd_l11','alos2_hbq_l11', 'chyaan2_fp','rs2_fp',  
    #
    'fp_sign'
    # SPECKEL FILTERS
    'rlee', 'boxcar',
    # UTILS
    'convert_T3_C3', 'convert_C3_T3', 'pauliRGB', 'convert_S2_CT', 
    # FULL-POL
    'grvi', 'rvifp', 'mf3cf', 'mf4cf', 'dopfp', 'prvifp', 'neufp', 'nnedfp', 
    'halphafp', 'shannon_h_fp','yam4cfp', 'halpha_plot_fp', 
    # COMPACT-POL
    'cprvi', 'dopcp', 'misomega', 'mf3cc',                 
    # DUAL-CROSS-POL
    'dprvi', 'dopdp', 'prvidp', 'rvidp', 'halphadp', 
    'shannon_h_dp', 'halpha_plot_dp',    
    'dprvic','dp_desc',
    # DUAL-CO-POL
    'mf3cd'      
    
]