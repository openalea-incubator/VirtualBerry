import os
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from virtual_berry import virtual_berry

import pandas


data = pandas.read_csv('DATA_input.txt', sep='\t')
virtual_berry(data)
