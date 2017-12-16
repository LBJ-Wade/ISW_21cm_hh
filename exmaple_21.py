import numpy as np
from run import *
#from cl_21_modified import *
from cl_21 import *
z_m_list = [30, 50,75,100,125,150,175,200]
w_list = [3.20283, 2.4956, 2.0427, 1.7706, 1.5841, 1.4459, 1.3383, 1.2513]
w_list = [3.1891, 2.4779, 2.0211, 1.7456, 1.5560, 1.4151, 1.3050, 1.2156]      # 4sigma 800 Mpc
infile = ['params_Yp_BBN/m0.dat','params_Yp_BBN/m1.dat','params_Yp_BBN/m2.dat']

infile = ['err_bar_plot/m1.dat','err_bar_plot/m2.dat']
tag_list = ['m1','m2']

infile = ['err_bar_plot/m0.dat', 'err_bar_plot/Neff1.dat']
tag_list = ['m0', 'Neff1']

infile = ['err_bar_plot/Neff1.dat']
tag_list = ['Neff1']

infile = ['err_bar_plot/Neff2.dat']
tag_list = ['Neff2']
for i in range(len(tag_list)):
	tag = tag_list[i]
	params_list = np.loadtxt (infile[i])[0:,]
	params_input = infile[i]
	
	run_21cm (params_list, params_input, tag, Yp_BBN=True)
	Cl = set_cl_21 (tag, Yp_BBN=True)
	plt.show()
	l_list = Cl.l_list
	cl21T_list = []
	cl21_list = []
	for j in range(len(z_m_list)):
		new_l_list = np.arange(2,5001)
		cl_zm = Cl.cl21T (z_m_list[j], w_list[j])
		cl_zm = np.interp (new_l_list, l_list, cl_zm)
		cl21T_list.append (cl_zm)
		zz = [z_m_list[j],z_m_list[j]]
		ww = [w_list[j], w_list[j]]
		cl_zmzm = Cl.cl21 (zz, ww)
		cl_zmzm = np.interp (new_l_list, l_list, cl_zmzm)
		cl21_list.append (cl_zmzm)
	out_zm = 'err_bar_plot/cl21_{}.txt'.format(tag)
	data = np.column_stack((new_l_list, cl21_list[0], cl21_list[1], cl21_list[2], cl21_list[3], cl21_list[4], cl21_list[5], cl21_list[6], cl21_list[7]))
	np.savetxt(out_zm, data, fmt = '%1.6e')
	out_zm = 'err_bar_plot/cl21T_{}.txt'.format(tag)
	data = np.column_stack((new_l_list, cl21T_list[0], cl21T_list[1], cl21T_list[2], cl21T_list[3], cl21T_list[4], cl21T_list[5], cl21T_list[6], cl21T_list[7]))
	np.savetxt(out_zm, data, fmt = '%1.6e')
