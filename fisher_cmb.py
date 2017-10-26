from fisher import *
# Fisher CMB
F = prior_cmb ()
F.cmb_deriv_vec ()
F.cov_matrix ()
fisher_matrix = F.cmb_fisher_analysis ()

data = np.column_stack((fisher_matrix[0],fisher_matrix[1],fisher_matrix[2],fisher_matrix[3],fisher_matrix[4],fisher_matrix[5],fisher_matrix[6],fisher_matrix[7]))
np.savetxt('fisher_matrix.txt', data, fmt = '%1.6e')

inv_fisher = inv(fisher_matrix)

sigma = []
for i in range(len(fisher_matrix)):
	sigma.append (np.sqrt(inv_fisher[i,i]))
	sigma = np.array(sigma)
data = np.column_stack((sigma))
np.savetxt('sigma.txt', data, fmt = '%1.6e')
