# thermofun-jupyter
ThermoFun examples using Jupyter notebooks

## Try ThermoFun in your browser

Using Jupyter Lab Notebooks and Binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/thermohub/thermofun-jupyter/master)

Wait until the Jupyter Lab Notebook server starts (~1 min) then double click on any `how-to-...` tutorial notebook.

More information on Jupyter Notebooks: [Jupyter Documentation](https://jupyter.readthedocs.io/en/latest/index.html)

## Run ThermoFun locally

Run Jupyter Notebooks examples/tutorials on jupyter notebook/lab running locally on your computer. For this you first need to have [Conda](https://conda.io/docs/) or [Miniconda](https://conda.io/miniconda.html) installed (follow the instructions on their websites).

In the Linux/MacOS terminal / windows conda command prompt or terminal change to the directory thermo-fun jupyter (that you cloned from this repository) and execute the following commands: 

```
conda env create -f environment.yml
conda activate thermofun-jupyter
conda install conda-forge::jupyterlab
``` 

First command will create the `thermofun-jupyter` conda environment congaing all the necessary libraries to run thermofun and associated codes, thermohubclient, chemicalfun, xgems, etc. The second command will activate the environment and your terminal line should show `(thermofun-jupyter)`. The third command installs jupyter lab in this environment.

To start jupyter lab finally execute:

```
jupyter lab
``` 

## ThermoFun examples

This repository contains examples of uses cases and functionality of ThermoFun. Additional information about the methods and data used in ThermoFun you can find [here](https://thermohub.org/thermofun/thermofun) or on the github [repository](https://github.com/thermohub/thermofun).

Start with how-to-use-thermofun-examples.ipynb or look in the folder how-to-use-thermofun-examples. 
