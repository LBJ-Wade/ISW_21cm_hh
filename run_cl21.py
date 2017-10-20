from cl_21 import *

infile = 'data/file_names.txt'
file_names = np.genfromtxt(infile, dtype="str")[0:]

params_input = file_names[0]
infile_syn = file_names[1]
infile_new = file_names[2]
infile_HYREC = file_names[3]
outfile = file_names[4]


params_list = np.loadtxt (params_input)[0:,]

z_m_list = [30, 50,75,100,125,150,175,200]
w_list = [3.31685, 2.47355 ,1.93014, 1.60434, 1.38246, 1.21989, 1.09467, 1]

aa = 2.7255*10**6
if not infile_new == 'None':
	Cl = cl_21 (params_list, infile_HYREC, infile_syn, infile_new)
	Cl.test_run ()
	Cl.c_z ()

	cl_list = []
	for i in [0,1]:#range(len(z_m_list)):
		z_m = z_m_list[i]
		w = w_list[i]
		cl  = Cl.cl21T (z_m, w)
		cl_list.append (cl)
	print (cl)
	l_list = Cl.l_list
	data = np.column_stack ((l_list, cl_list[0], cl_list[1]))
	np.savetxt(outfile, data, fmt = '%1.6e')

else:
	Cl = cl_21 (params_list, infile_HYREC, infile_syn)

	Cl.test_run ()
	Cl.c_z ()
	
	cl_list = []
	for i in range(2):
		for j in range(2):
			z_m = [z_m_list[i], z_m_list[j]]
			w = [w_list[i], w_list[j]]
			cl = Cl.cl21 (z_m, w)
			cl_list.append (cl)
	l_list = Cl.l_list
	data = np.column_stack ((l_list, cl_list[0], cl_list[1], cl_list[2], cl_list[3]))
	np.savetxt(outfile, data, fmt = '%1.6e')


#plt.plot (l_list, np.sqrt(aa*cl*l_list*(l_list+1)/(2*np.pi)))
#plt.xscale ('log')
#plt.yscale ('log')

#plt.show()

