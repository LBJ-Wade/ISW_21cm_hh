import numpy as np
infile1 = "fisher_matrix_BBN_3000_fsky0.4.txt"
infile2 = "fisher_matrix_cl21T_BBN_3000_fsky0.4.txt"
#infile1 = "fisher_matrix_BBN_5000_fsky0.7.txt"
#infile2 = "fisher_matrix_cl21T_BBN_5000_fsky0.7.txt"
infile_planck = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/prior_planck_logAs.txt"
prior_planck = np.zeros([8,8])
for i in range(6):
	a = np.loadtxt (infile_planck)[0:,i]
	prior_planck[i][:6] = a.copy()
F = np.zeros([8,8])
for i in range(8):
	F[i,:] = np.loadtxt(infile1)[0:,i]+np.loadtxt(infile2)[0:,i]

F = F + prior_planck

data = np.column_stack((F[0],F[1],F[2],F[3],F[4],F[5],F[6],F[7]))
#np.savetxt('Planck_S4_BBN_3000_fsky0.4.txt', data, fmt = '%f')
np.savetxt('Planck_S4_21_BBN_5000_fsky0.7.txt', data, fmt = '%f')

