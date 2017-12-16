from fisher import *
# Fisher 21T
F = fisher (Yp_BBN = True)
F.cl21T_deriv_vec ()
F.cov_matrix ()
fisher_matrix = F.fisher_analysis ()

data = np.column_stack((fisher_matrix[0],fisher_matrix[1],fisher_matrix[2],fisher_matrix[3],fisher_matrix[4],fisher_matrix[5],fisher_matrix[6],fisher_matrix[7]))
#data = np.column_stack((fisher_matrix[0]))
np.savetxt('result/fisher_matrix_cl21T_BBN.txt', data, fmt = '%1.20e')

inv_fisher = inv(fisher_matrix)

sigma = []
for i in range(len(fisher_matrix)):
	sigma.append (np.sqrt(inv_fisher[i,i]))
sigma = np.array(sigma)
data = np.column_stack((sigma))
np.savetxt('result/sigma_cl21T_BBN.txt', data, fmt = '%1.6e')
plt.show()
