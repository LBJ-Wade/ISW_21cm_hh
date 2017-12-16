from fisher import *
from mpi4py import MPI

print ('\n')

comm = MPI.COMM_WORLD
rank = comm.Get_rank ()
size = comm.Get_size ()

# Fisher 21T
F = fisher (comm, rank, size, Yp_BBN = True)
F.cl21T_deriv_vec ()
F.cov_matrix ()
if rank == 0:
	fisher_matrix = F.fisher_analysis ()

	data = np.column_stack((fisher_matrix[0],fisher_matrix[1]))
	np.savetxt('result/fisher_matrix_cl21T_BBN.txt', data, fmt = '%1.20e')

	inv_fisher = inv(fisher_matrix)

	sigma = []
	for i in range(len(fisher_matrix)):
		sigma.append (np.sqrt(inv_fisher[i,i]))
	sigma = np.array(sigma)
	data = np.column_stack((sigma))
	np.savetxt('result/sigma_cl21T_BBN.txt', data, fmt = '%1.6e')
