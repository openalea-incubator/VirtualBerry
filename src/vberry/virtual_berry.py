import os
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from vberry.parameters import params, controls

dir = os.path.dirname(__file__)
r_source = ro.r['source']
r_source(dir + '/virtual_berry.R')

def virtual_berry(data, parameters=params, controls=controls):
    """ call R virtual berry with panda daframe and parameter dicts"""
    with (ro.default_converter + pandas2ri.converter).context():
        r_df = ro.conversion.get_conversion().py2rpy(data)
    res = ro.r['modgrowth'](DATA=r_df, **parameters, **controls)
    
    with (ro.default_converter + pandas2ri.converter).context():
        pd_df = ro.conversion.get_conversion().rpy2py(res)
        
    return(pd_df)
