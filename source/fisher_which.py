from fisher import *

def run (cmb_or_21, Yp_BBN):
	
	if cmb_or_21 == '21':
			
		F = fisher (Yp_BBN)
		F.cl21T_deriv_vec ()
		F.cov_matrix ()
		fisher_matrix = F.fisher_analysis ()
		
		if Yp_BBN == True:
			data = np.column_stack ((fisher_matrix[0], fisher_matrix[1], fisher_matrix[2], \
			                         fisher_matrix[3], fisher_matrix[4], fisher_matrix[5], \
									 fisher_matrix[6], fisher_matrix[7]))
			np.savetxt('result/fisher_matrix_cl21T_BBN.txt', data, fmt = '%1.20e')

		elif Yp_BBN == False:
			data = np.column_stack ((fisher_matrix[0], fisher_matrix[1], fisher_matrix[2], \
			                         fisher_matrix[3], fisher_matrix[4], fisher_matrix[5], \
									 fisher_matrix[6], fisher_matrix[7], fisher_matrix[8]))
			np.savetxt('result/fisher_matrix_cl21T.txt', data, fmt = '%f')
		
		else:
			print 'Keyword Error in run_fisher.py'
			return None

		inv_fisher = inv(fisher_matrix)
		sigma = []
		for i in range(len(fisher_matrix)):
			sigma.append (np.sqrt(inv_fisher[i,i]))
		sigma = np.array(sigma)
		data = np.column_stack((sigma))
		
		if Yp_BBN == True:	
			np.savetxt('result/sigma_cl21T_BBN.txt', data, fmt = '%1.6e')
		elif Yp_BBN == False:
			np.savetxt('result/sigma_cl21T.txt', data, fmt = '%1.6e')

	elif cmb_or_21 == 'cmb':
	
		F = prior_cmb (Yp_BBN)
		F.cmb_deriv_vec ()
		F.cov_matrix ()
		fisher_matrix = F.cmb_fisher_analysis ()
	
		if Yp_BBN == True:
			data = np.column_stack ((fisher_matrix[0], fisher_matrix[1], fisher_matrix[2], \
			                         fisher_matrix[3], fisher_matrix[4], fisher_matrix[5], \
									 fisher_matrix[6], fisher_matrix[7]))
			np.savetxt('result/fisher_matrix_BBN.txt', data, fmt = '%f')
	
		elif Yp_BBN == False:
			data = np.column_stack ((fisher_matrix[0], fisher_matrix[1], fisher_matrix[2], \
			                         fisher_matrix[3], fisher_matrix[4], fisher_matrix[5], \
									 fisher_matrix[6], fisher_matrix[7], fisher_matrix[8]))
			np.savetxt('result/fisher_matrix.txt', data, fmt = '%f')

		else:
			print 'Keyword Error in run_fisher.py'
			return None
		
		inv_fisher = inv(fisher_matrix)
		sigma = []
		for i in range(len(fisher_matrix)):
			sigma.append (np.sqrt(inv_fisher[i,i]))
		sigma = np.array(sigma)
		data = np.column_stack((sigma))

		if Yp_BBN == True:
			np.savetxt('result/sigma_BBN.txt', data, fmt = '%1.6e')
		elif Yp_BBN == False:
			np.savetxt('result/sigma.txt', data, fmt = '%1.6e')

	else:
		print 'Keyword Error in run_fisher.py'
		return None

