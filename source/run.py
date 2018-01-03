import numpy as np
import os
from subprocess import call
from path import *
from updateparams import *

def run_cmb (params_list, tag, Yp_BBN = True):

	# Run CLASS (syncronous gauge)
	os.chdir (path_CLASS_syn)
	oldfile = "params_prac_.ini"
	newfile = "params_prac2_.ini"
	clfile = "params_prac_cl.dat"
	clout = path_data + "/cl_" + tag + ".dat"
	setparamsfile_CLASS (params_list, oldfile, newfile, Yp_BBN)
	call ('./class {0} {1}'.format (newfile, "cl_ref.pre"), shell = True)	
	os.chdir ("output")
	os.system ("cp {0} {1}".format (clfile, clout))

	os.system ("rm *")
	os.chdir ("../")
	os.chdir ("../")
	
	return clout

def run_21cm (params_list, params_input, tag, Yp_BBN = True):
	
	
	#Delete unnecessary output file in CLASS
	#os.system ("rm class_syn/output/*")
	#os.system ("rm class_new/output/*")

	# Run CLASS (syncronous gauge)
	os.chdir (path_CLASS_syn)
	oldfile = "params_prac_.ini"
	newfile = "params_prac2_.ini"
	clfile = "params_prac_cl.dat"
	clout = path_data + "/params_prac_cl_" + tag + ".dat"
	setparamsfile_CLASS (params_list, oldfile, newfile, Yp_BBN)
	call ('./class {0}'.format (newfile), shell = True)	
	
	os.chdir ("output")
	outfile_syn = path_data + "/delta_syn_{0}.dat".format (tag)
	outfile_HYREC = path_data + "/output_prac_{0}.dat".format (tag)
	os.system ("cp delta.dat {0}".format (outfile_syn))
	os.system ("cp HyRec_output.dat {0}".format (outfile_HYREC))
	os.system ("cp {0} {1}".format (clfile, clout))
	if tag == "0":
		if Yp_BBN == True:
			cl_out = path_result_Yp_BBN + "/cl_" + tag + ".dat"
		else:
			cl_out = path_result_Yp + "/cl_" + tag + ".dat"
		os.system ("cp params_prac_cl.dat {1}".format(tag, cl_out))

	os.system ("rm *")
	os.chdir ("../")
	os.chdir ("../")
	
	# Run CLASS (newtonian gauge)
	os.chdir (path_CLASS_new)
	oldfile = "params_prac_.ini"
	newfile = "params_prac2_.ini"
	setparamsfile_CLASS (params_list, oldfile, newfile, Yp_BBN)
	call ('./class  {0}'.format (newfile), shell = True)	
	os.chdir ("output")
	outfile_new = path_data + "/delta_new_{0}.dat".format (tag)
	os.system ("cp delta.dat {0}".format (outfile_new))
	os.system ("rm *")
	os.chdir ("../")
	os.chdir ("../")
	
	'''
	os.chdir (path_HYREC)
	# Run HYREC	
	oldfile = "input_prac.dat"
	newfile = "input_prac2.dat"
	outfile_HYREC = "output_prac_{0}.dat".format (tag)
	setparamsfile_HYREC (params_list, oldfile, newfile)
	os.system ("gcc -lm -O3 hyrectools.c helium.c hydrogen.c history.c hyrec.c -o hyrec")
	call ("./hyrec < {0} > {1}".format (newfile, outfile_HYREC), shell = True)	
	os.chdir ("../")
	'''

	# Calculate C_l
	file_names = np.array ([params_input, outfile_syn, outfile_new, outfile_HYREC])
	np.savetxt (path_data + '/file_names_{0}.txt'.format (tag), file_names, fmt="%s")
		
# Load cosmological parameters
#params_list = np.loadtxt (params_input)[0:,]
#tag = "baryon"
#run_21cm (params_list, tag)
#run (params_list, '0')
