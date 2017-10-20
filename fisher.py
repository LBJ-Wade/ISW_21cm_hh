from run import *
from numpy.linalg import inv
from cl_21 import *

class fisher (object): 
	def __init__ (self):
		
		self.n = 2				# Number of window functions
		self.stepsize = [0.01, 0.005, 0.015]
		self.params_list = np.loadtxt (params_input)[0:,]
		self.fisher_params = ["h","b"]		# The order of parameters should be the same as in params_list
		
		self.z_m_list = [30, 50]
		self.w_list = [3.31685, 2.47355]

		self.deriv_vec = {}
		self.cov = {}

	def cl21T_vec (self):
		""" Construct covariance matrix """
		
		
		for j in range(len(self.fisher_params)):
			for i in [0]:#range(len(self.stepsize)):
				param = self.fisher_params[j]

				params_list_copy = self.params_list.copy ()
				params_list_copy[j] *= (1 - self.stepsize[i])
				tag = param + "_{}1".format(i)
				run (params_list_copy, tag)
				Cl1 = set_cl_21 (tag)
				#outfile1 = path_result + "/cl21T_{}.txt".format (tag)
				#outfile1_list.append (outfile1)
				
				params_list_copy = self.params_list.copy ()
				params_list_copy[j] *= (1 + self.stepsize[i])
				tag = param + "_{}2".format(i)
				run (params_list_copy, tag)
				Cl2 = set_cl_21 (tag)
				#outfile2 = path_result + "/cl21T_{}.txt".format (tag)
				#outfile2_list.append (outfile2)
				dev_cl = {}
				for k in range(len(self.z_m_list)):
					cl1_zm = Cl1.cl21T (self.z_m_list[k], self.w_list[k])
					cl2_zm = Cl2.cl21T (self.z_m_list[k], self.w_list[k])
					dev_cl_zm = (cl2_zm-cl1_zm)/(2*self.params_list[j]*self.stepsize[i])
					dev_cl['{0}'.format (self.z_m_list[k])] = dev_cl_zm
			self.deriv_vec[param] = dev_cl


	def cov_matrix (self):
		""" Construct derivative of covariance matrix """

		cl21T = {}
		cl21 = {}
		tag = "0"
		run (self.params_list, tag)
		Cl = set_cl_21 (tag)
		for i in range(len(self.z_m_list)):
			cl_zm = Cl.cl21T (self.z_m_list[i], self.w_list[i])
			cl21T["{0}".format(i,i)] = cl_zm
		for i in range(len(self.z_m_list)):
			for j in range(len(self.z_m_list)):
				if not j < i:
					zm = [self.z_m_list[i], self.z_m_list[j]]
					w = [self.w_list[i], self.w_list[j]]
					cl_zmzl = Cl.cl21 (zm, w)
					if i == j:
						cl21["{0}{1}".format(i,j)] = cl_zmzl
					else:
						cl21["{0}{1}".format(i,j)] = cl_zmzl
						cl21["{0}{1}".format(j,i)] = cl_zmzl
		
		cl_out = path_result + "/cl_" + tag + ".dat"
		l = np.loadtxt(cl_out)[0:,0]
		aa = 2.7255**2 * 2*np.pi / (l*(l+1))
		cl = np.loadtxt(cl_out)[0:,1]*aa
		self.l_list = Cl.l_list
		clTT = np.interp (Cl.l_list, l, cl)
		for i in range(len(self.z_m_list)):
			for j in range(len(self.z_m_list)):
				element_ij = (cl21["{0}{1}".format(i,j)] * clTT + cl21T["{0}".format(i)]*cl21T["{0}".format(j)]) / (2*Cl.l_list+1)
				self.cov["{0}{1}".format(i,j)] = element_ij
		
	def fisher_analysis (self):
		F = np.zeros([len(self.fisher_params),len(self.fisher_params)])
		for m in range(len(self.fisher_params)):
			for n in range(len(self.fisher_params)):
				F_mn = 0
				param_m = self.fisher_params[m]
				param_n = self.fisher_params[n]
				for l in range(len(self.l_list)):
					vec_m = np.zeros(len(self.z_m_list))
					vec_n = np.zeros(len(self.z_m_list))
					inv_cov = np.zeros([len(self.z_m_list), len(self.z_m_list)])
					for i in range(len(self.z_m_list)):
						for j in range(len(self.z_m_list)):
							vec_m[i] = self.deriv_vec[param_m]["{0}".format(self.z_m_list[i])][l]
							vec_n[j] = self.deriv_vec[param_n]["{0}".format(self.z_m_list[j])][l]
							inv_cov[i,j] = 1 / self.cov["{0}{1}".format(i,j)][l]
							F_mn += vec_m[i]*inv_cov[i,j]*vec_n[j]
				F[m,n] = F_mn
		return F


	def convergence_test (self):
		return None	

	
F = fisher ()
F.cl21T_vec ()
F.cov_matrix ()
fisher_matrix = F.fisher_analysis ()
print (fisher_matrix)
print (inv(fisher_matrix))
print (np.dot (fisher_matrix, inv(fisher_matrix)))
