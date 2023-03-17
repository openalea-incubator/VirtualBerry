from vberry.virtual_berry import virtual_berry
from vberry.parameters import params, controls
from vberry.data_samples import input_data
import pandas



# [option] edit params control as needed
data_path = input_data()
data = pandas.read_csv(data_path, sep='\t')
#
res = virtual_berry(data, params, controls)