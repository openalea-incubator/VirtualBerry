# Virtual Berry

Simulation of grapevine, virtual berry growth and development, model described in [this paper](https://github.com/openalea-incubator/VirtualBerry/blob/main/doc/2008_Dai_ActaHort.pdf)


## Installation 

At first you make sure that you have installed the package manager conda in your system: follow the instructions at https://docs.conda.io/en/latest/miniconda.html

Then we create a conda environment:

    conda create -n vberry -c conda-forge python=3.8 pytest ipython
	conda activate vberry
 
If you have an active R then you can skip the R installation (check that package desolve is installed). Otherwise

	conda install -c conda-forge R r-desolve
	
Note that if you have an existing R installation, R version should be <4.2 to be compatible with rpy2.

Find the path to R by typing in a console :

	R.home()
	
Define an environment variable R_HOME pointing to this path.

Install rpy2:

	conda install -c conda-forge pandas rpy2
	

Clone project and run setup:

    git clone https://github.com/openalea-incubator/VirtualBerry.git
    cd VirtualBerry
	python setup.py develop

## Usage

A basic example of use is in example/tutorial.py

