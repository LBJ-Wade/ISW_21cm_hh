from cl_21 import *

infile_syn = sys.argv[1]
infile_new = sys.argv[2]
infile_21 = sys.argv[3]
outfile = sys.argv[4]

z_m_list = [30, 50,75,100,125,150,175,200]
w_list = [3.31685, 2.47355 ,1.93014, 1.60434, 1.38246, 1.21989, 1.09467, 1]

Cl = cl_21 (infile_syn, infile_new, infile_21)
aa = 2.7255*10**6
cl_list = []
for i in [0]:#range(len(z_m_list)):
	z_m = z_m_list[i]
	w = w_list[i]

	cl  = Cl.cl21T (z_m, w)
	cl_list.append (cl)

l_list = Cl.l_list
#data = np.column_stack ((l_list, cl_list[0], cl_list[1], cl_list[2], cl_list[3], cl_list[4], cl_list[5], cl_list[6], cl_list[7]))
data = np.column_stack ((l_list, cl_list[0]))
np.savetxt(outfile, data, fmt = '%1.6e')

#plt.plot (l_list, np.sqrt(aa*cl*l_list*(l_list+1)/(2*np.pi)), label = 'limber2')
#plt.plot (l_list, np.sqrt(aa*cl2*l_list*(l_list+1)/(2*np.pi)), label = 'limber')
#plt.legend ()
#plt.xscale ('log')
#plt.yscale ('log')

#plt.show()

