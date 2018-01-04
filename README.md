# ISW_21cm
This program is a part of Master research advised by Prof. Marilena Loverde at Stony Brook University. The goal of this program is constraining cosmological parameters by Fisher analysis using the cross-correlation functions of Integrated Sachs-Wolfe effect and 21cm brightness temperature fluctuations at dark ages. A paper containing results of this program is in preparation.

## Getting started
* This program is written with python 2.7.
* 'default' in /source/path.py should be changed to be a path of 'ISW_21cm' folder in the machine.
* Should do 'make' in both 'class_syn' and 'class_new' folders.

## List of folders

* `class_syn`
  - CLASS in syncronouse gauge
* `class_new`
  - CLASS in newtonian gauge
* `data`
  - This is a folder to save data from CLASS such as transfer functions, CMB power spectrums, and free electron fraction.   
* `params_Yp`
  - Sets of cosmological parameters for Fisher analysis. Yp is also one of parameters.
* `params_Yp_BBN`
  - Sets of cosmological parameters for Fisher analysis. Yp is fixed by BBN constraint. See Eq. (70) in Planck 2015 results. XIII.
* `result`
  - In this folder, the fisher matrix and the constraints from Fisher analysis are saved.
* `result_Yp`
  - The auto- and cross-correlation of 21cm fluctuations and ISW effect are saved using parameters in params_Yp.
* `result_Yp_BBN`
  - The auto- and cross-correlation of 21cm fluctuations and ISW effect are saved using parameters in params_Yp_BBN.
* `source`
  - Contains main source files.

## Source files
* `source/`
  - `path.py`: this defines paths of every folder.

  - `selection_function21cm.py`: this defines selection functions for the obversed frequency of 21 cm fluctuations.
  
  - `run.py`: this runs CLASS code to get transfer functions, CMB power spectrums, and free electron fraction. It has two functions each of which is for Fisher analysis with CMB or 21cm.
  
  - `cl_21.py`: this calculates auto- and cross-correlation functions of 21cm using data from CLASS.
  
  - `fisher.py`: this does Fisher analysis for CMB and 21cm.
  
  - `fisher_which.py`: this determines which Fisher analysis will be done (CMB or 21cm) and which conditions for Yp (BBN constraint or not) will be considered.
  
  - `tools.py`: this contains functions to interpolate recombination coefficient from table given in HyRec code. Its original file is hydrogen.c of HyRec code.
  
  - `updateparams.py`: this updates CLASS ini file (params_prac_.ini) using parameters in `params_Yp` and `params_Yp_BBN`.

## How to run

* `run_fisher.py`: choose between CMB and 21cm and also condition for Yp. Result (fisher matrix and constraints of cosmological parameters) will be saved in `result` folder.

