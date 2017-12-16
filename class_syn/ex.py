from subprocess import call
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank ()
size = comm.Get_size ()

infile = "params_prac2_.ini"

if rank == 0:
	os.chdir ('class_syn')
	call ('./class {0} > log.txt'.format(infile), shell = True)
elif rank == 1:

	os.chdir ('class_new')
	call ('./class {0} > log.txt'.format(infile), shell = True)
