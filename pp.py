from updateparams import *
import numpy as np
import os
from subprocess import call
from path import *

infile = 'err_bar_plot/Neff2.dat'
params_list = np.loadtxt(infile)[0:,]

oldfile = "params_prac_.ini"
newfile = "params_prac2_.ini"

os.chdir ('class_syn')
setparamsfile_CLASS (params_list, oldfile, newfile, Yp_BBN = True)
call ('./class {0} '.format (newfile, "cl_ref.pre"), shell = True)
