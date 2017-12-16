

import numpy as np
import os
from subprocess import call
from path import *
from updateparams import *


def run_21cm (comm, rank, size, params_list, params_input, tag, Yp_BBN = True):
	
	# Delete unnecessary output file in CLASS
	os.system ("rm class_syn/output/*")
	os.system ("rm class_new/output/*")
	
	# Run CLASS
	print ('-> (Rank:{0}) Start running CLASS'.format(rank))
	if size == 1:
		run_CLASS_syn (params_input, params_list, Yp_BBN, tag, rank)
		run_CLASS_new (params_input, params_list, Yp_BBN, tag, rank)
	else:
		if rank == 0:
			run_CLASS_syn (params_input, params_list, Yp_BBN, tag, rank)
		elif rank == 1:
			run_CLASS_new (params_input, params_list, Yp_BBN, tag, rank)
	
	comm.Barrier ()
	outfile_syn = path_data + "/delta_syn_{0}.dat".format (tag)
	outfile_new = path_data + "/delta_new_{0}.dat".format (tag)
	outfile_HYREC = path_data + "/output_prac_{0}.dat".format (tag)
	
	# Save file names to load data in cl_21.py
	file_names = np.array ([params_input, outfile_syn, outfile_new, outfile_HYREC])
	np.savetxt (path_data + '/file_names_{0}.txt'.format (tag), file_names, fmt="%s")

def run_CLASS_syn (params_input, params_list, Yp_BBN, tag, rank):
	# Run CLASS (syncronous gauge)
	
	os.chdir (path_CLASS_syn)
	oldfile = "params_prac_.ini"
	newfile = "params_prac2_.ini"
	clfile = "params_prac_cl.dat"
	clout = path_data + "/params_prac_cl_" + tag + ".dat"
	
	print ('-> (Rank:{0}) Updating CLASS ini file (syncronous gauge) with parameter file:'.format(rank))
	print ('-> (Rank:{0}) {1}'.format(rank, params_input))
	setparamsfile_CLASS (params_input, params_list, oldfile, newfile, Yp_BBN)
	print ('-> (Rank:{0}) Running CLASS in syncronous gauge'.format(rank))
	call ('./class {0} > log.txt'.format (newfile), shell = True)	
	
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
	print ('   (Rank:{0}) Done'.format(rank))

def run_CLASS_new (params_input, params_list, Yp_BBN, tag, rank):
	# Run CLASS (newtonian gauge)
	
	os.chdir (path_CLASS_new)
	oldfile = "params_prac_.ini"
	newfile = "params_prac2_.ini"
	
	print ('-> (Rank:{0}) Updating CLASS ini file (newtonian gauge) with parameter file:'.format(rank))
	print ('-> (Rank:{0}) {1}'.format(rank, params_input))
	setparamsfile_CLASS (params_input, params_list, oldfile, newfile, Yp_BBN)
	print ('-> (Rank:{0}) Running CLASS in newtonian gauge'.format(rank))
	call ('./class {0} > log.txt'.format (newfile), shell = True)	
	os.chdir ("output")
	outfile_new = path_data + "/delta_new_{0}.dat".format (tag)
	os.system ("cp delta.dat {0}".format (outfile_new))
	os.system ("rm *")
	os.chdir ("../")
	os.chdir ("../")
	print ('   (Rank:{0}) Done'.format(rank))

