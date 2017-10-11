import numpy as np
import os
from subprocess import call
from path import *
from updateparams import *

def run ():
	# Load cosmological parameters
	params_list = np.loadtxt (params_input)[0:,]
		
	# Run CLASS (syncronous gauge)
	os.chdir(path_CLASS_syn)
	oldfile = "params_prac_.ini"
	newfile = "params_prac2_.ini"
	setparamsfile_CLASS (params_list, oldfile, newfile)
	call ('./class  {0}' .format(newfile), shell = True)	
	os.chdir ("output")
	outfile_syn = path_data + "/delta_syn.dat"
	os.system ("cp delta.dat {0}".format (outfile_syn))
	os.system ("rm *")
	os.chdir ("../")
	os.chdir ("../")
	
	# Run CLASS (newtonian gauge)
	os.chdir(path_CLASS_new)
	oldfile = "params_prac_.ini"
	newfile = "params_prac2_.ini"
	setparamsfile_CLASS (params_list, oldfile, newfile)
	call('./class  {0}' .format(newfile), shell = True)	
	os.chdir ("output")
	outfile_new = path_data + "/delta_new.dat"
	os.system ("cp delta.dat {0}".format (outfile_new))
	#os.system("cp delta.dat {0}/delta_new.dat".format (path_data))
	os.system ("rm *")
	os.chdir ("../")
	os.chdir ("../")
	
	# Run HYREC	
	os.chdir(path_HYREC)
	oldfile = "input_prac.dat"
	newfile = "input_prac2.dat"
	outfile_HYREC = "output_prac.dat"
	setparamsfile_HYREC (params_list, oldfile, newfile)
	os.system ("gcc -lm -O3 hyrectools.c helium.c hydrogen.c history.c hyrec.c -o hyrec")
	call("./hyrec < {0} > {1}" .format (newfile, outfile_HYREC), shell = True)	
	
	# Calculate 21cm fluctuation coefficients
	outfile_21 = path_data + "/transfer21.txt"
	os.system ("python c_z.py {0} {1} {2}".format (outfile_syn, outfile_HYREC, "transfer21.txt"))
	os.system ("cp transfer21.txt {0}".format (outfile_21))

	# Calculate C_l^{21,ISW}
	outfile_cl21 = path_result + "/cl21T.txt"
	os.system ("python cl_21.py {0} {1} {2} {3}".format (outfile_syn, outfile_new, outfile_21, outfile_cl21)
	)
run ()

