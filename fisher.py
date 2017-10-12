from run import *


stepsize = [0.01, 0.005, 0.015]
params_list = np.loadtxt (params_input)[0:,]
fisher_params = ["h"]


for param in fisher_params:
	tag = param + "_0"
	run (params_list, tag)

	for i in range(len(stepsize)):
		params_list_copy = params_list.copy ()
		params_list_copy[0] *= (1 - stepsize[i])
		tag = param + "_{}1".format(i)
		run (params_list_copy, tag)
	
		params_list_copy = params_list.copy ()
		params_list_copy[0] *= (1 + stepsize[i])
		tag = param + "_{}2".format(i)
		run (params_list_copy, tag)

