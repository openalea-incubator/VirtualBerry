import os


datadir = os.path.join(os.path.dirname(__file__), 'data')


def input_data():
    return os.path.join(datadir, 'DATA_input.txt')