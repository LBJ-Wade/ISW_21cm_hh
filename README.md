# ISW_21cm
This program is a part of Master research advised by Prof. Marilena Loverde at Stony Brook University. The goal of this program is constraining cosmological parameters by Fisher analysis using the cross-correlation functions of Integrated Sachs-Wolfe effect and 21cm brightness temperature fluctuations at dark ages.

## Getting started
* This program is written with python 2.7.
* 'default' in /source/path.py should be changed to be a path of 'ISW_21cm' folder in the machine.
* Should do 'make' in both 'class_syn' and 'class_new' folder.

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

## Structure and source files.
* `source/`
  - `path.py`:

## How to run

* `run_fisher.py`

