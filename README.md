# Virtual Berry

Simulation of grapevine vberry growth and development


If you have an existing R installation, R version should be <4.2 to be compatible with rpy2

## Installation from source

    conda create -n vberry -c conda-forge python=3.8 pytest ipython
	conda activate vberry
	# skip R installation is you have an active R
	conda install -c conda-forge R r-desolve pandas rpy2
	
find the path to R by typing in a console :

	R RHOME
	
define an environment variable R_HOME pointing to this path

clone project and run setup

    git clone https://github.com/openalea-incubator/VirtualBerry.git
    cd VirtualBerry
	python setup.py develop

## Usage

A basic example of use is in example/tutorial.py

