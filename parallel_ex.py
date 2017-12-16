from cl_21 import *
from path import *
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank ()
size = comm.Get_size ()

Cl = set_cl_21(comm, rank, size,'0',Yp_BBN = True)
l_list = np.arange(2,5001)
cl21T_list = []
#z_m_list = [30, 50,75,100,125,150,175,200]
#w_list = [3.1891, 2.4779, 2.0211, 1.7456, 1.5560, 1.4151, 1.3050, 1.2156]

z_m_list = [30]
w_list = [3.1891]

for i in range(len(z_m_list)):
	a = time.time()

	# Cross-correlation
	cl_zm = Cl.cl21T (z_m_list[i], w_list[i])
	b = time.time ()	
	# Auto-correlation
	cl_zm = Cl.cl21 ([z_m_list[i],z_m_list[i]], [w_list[i],w_list[i]])
	c = time. time()

	cl21T_list.append (cl_zm)
print 'Cross-correlation running time', b-a
print 'Auto-correlation running time', c-b

if rank == 0:
	out_zm = path_result_Yp_BBN + '/cl21_parallel{}.txt'.format(rank)
	data = np.column_stack ((Cl.l_list, cl21T_list[0]))
	np.savetxt (out_zm, data,  fmt = '%1.6e')

