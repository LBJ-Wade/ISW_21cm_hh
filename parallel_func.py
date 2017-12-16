from mpi4py import MPI
import numpy as np

def job_alloc (comm, rank, size, l_list):
	total_iteration = len (l_list)
	sub_iteration = int (total_iteration) / size
	
	if size == 1:
		return l_list
	else:
		if not rank == size-1:
			return l_list[sub_iteration*rank:sub_iteration*(rank+1)]
		else:
			return l_list[sub_iteration*rank:]

def combine_result (comm, rank, size, result):
	rank_list = range(size)
	if not rank == 0:
		comm.send (result, dest = 0)
	else:
		for i in rank_list[1:]:
			data = comm.recv (source = i)
			result = np.array (list(result) + list(data))
	return result


