from run import *
from numpy.linalg import inv

class fisher (object): 
	def __init__ (self):
		
		self.n = 2				# Number of window functions
		self.stepsize = [0.01, 0.005, 0.015]
		self.params_list = np.loadtxt (params_input)[0:,]
		self.fisher_params = ["h","b"]		# The order of parameters should be the same as in params_list

		self.cov = None
		self.deriv_cov_dic = {}

	def cov_matrix (self):
		""" Construct covariance matrix """
		
		#cov = np.zeros ([self.n, self.n])
		cov = {}
		#tel = {'jack': 4098, 'sape': 4139}
		tag = "_0"
		run (self.params_list, tag)
		
		outfile_cl21 = path_result + "/cl21T_{}.txt".format (tag)
		
		for i in range(self.n):
			if i == 0:
				self.l_list = np.loadtxt(outfile_cl21)[0:,i]
			cl = np.loadtxt (outfile_cl21)[0:,i+1]
			#cov[i,i] = cl
			cov["{0}{1}".format(i,i)] = cl
		self.cov = cov
	
	def deriv_cov_matrix (self):
		""" Construct derivative of covariance matrix """

		outfile1_list = []
		outfile2_list = []

		for j in range(len(self.fisher_params)):
			for i in [0]:#range(len(self.stepsize)):
				param = self.fisher_params[j]

				params_list_copy = self.params_list.copy ()
				params_list_copy[j] *= (1 - self.stepsize[i])
				tag = param + "_{}1".format(i)
				run (params_list_copy, tag)
				outfile1 = path_result + "/cl21T_{}.txt".format (tag)
				outfile1_list.append (outfile1)
				
				params_list_copy = self.params_list.copy ()
				params_list_copy[j] *= (1 + self.stepsize[i])
				tag = param + "_{}2".format(i)
				run (params_list_copy, tag)
				outfile2 = path_result + "/cl21T_{}.txt".format (tag)
				outfile2_list.append (outfile2)

				#deriv_cov = np.zeros([self.n, self.n])
				deriv_cov = {}
				for k in range(self.n):
					cl1 = np.loadtxt (outfile1)[0:,i+1]
					cl2 = np.loadtxt (outfile2)[0:,i+1]
					deriv_cl = (cl2-cl1)/ (self.params_list[j]*self.stepsize[0]*2)
					deriv_cov['{0}{1}'.format(k,k)] = deriv_cl
			self.deriv_cov_dic['{0}'.format(j)] = deriv_cov
		# For conv test
		self.outfile1_list = outfile1_list	
		self.outfile2_list = outfile2_list	
	def fisher_analysis (self):
		""" Fisher Analysis """
		
		sigma = []

		fisher_matrix = np.zeros ([len(self.fisher_params), len(self.fisher_params)])
		for i in range(len(self.fisher_params)):
			for j in range(len(self.fisher_params)):
				deriv_cov1 = self.deriv_cov_dic['{0}'.format(i)]
				deriv_cov2 = self.deriv_cov_dic['{0}'.format(j)]
				
				
				for l in range(len(self.l_list)):
					print (i,j,l)

					deriv_cov1_l = np.zeros([self.n, self.n])
					deriv_cov2_l = np.zeros([self.n, self.n])
					cov = np.zeros([self.n, self.n])
					for ii in range(self.n):
						for jj in range(self.n):
							if jj == ii:		#Diagonal for now
								deriv_cov1_l[ii,jj] = deriv_cov1['{0}{1}'.format(ii,jj)] [l]
								deriv_cov2_l[ii,jj] = deriv_cov2['{0}{1}'.format(ii,jj)] [l]
								cov[ii,jj] = self.cov['{0}{1}'.format(ii,jj)] [l]
					inv_cov = inv (cov)
					fisher_matrix [i,j] += 1/2 * np.trace(np.dot (deriv_cov1_l, np.dot (inv_cov , np.dot (deriv_cov2_l ,inv_cov))))
		inv_fisher_matrix = inv (fisher_matrix)
		for i in range(len(self.fisher_params)):
			sigma.append (inv_fisher_matrix[i,i])
		self.sigma = np.array (sigma)

		return self.sigma

	def convergence_test (self):
		return None	

	
F = fisher ()
F.cov_matrix ()
F.deriv_cov_matrix ()
sigma = F.fisher_analysis ()
print (F.fisher_params)
print (sigma)
