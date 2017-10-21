from run import *
from numpy.linalg import inv
from cl_21 import *

class prior_cmb (object):		# Need to add cl^dd
	def __init__ (self):
		self.stepsize = [0.01, 0.01]
		self.params_list = np.loadtxt (params_input)[0:,]
		self.fisher_params = ["h","b"]
		self.deriv_vec = {}
		self.cov = {}
		self.l_list = None
	def cmb_deriv_vec (self):
		
		for j in range(len(self.fisher_params)):
			param = self.fisher_params[j]
			stepsize = self.stepsize[j]

			params_list_copy = self.params_list.copy ()
			params_list_copy[j] -= stepsize 
			tag = param + "1"
			outfile1 = run_cmb (params_list_copy, tag)

			params_list_copy = self.params_list.copy ()
			params_list_copy[j] += stepsize 
			tag = param + "2"
			outfile2 = run_cmb (params_list_copy, tag)
			
			dev_cl = {}
			l = np.loadtxt (outfile1)[0:,0]
			clTT1 = np.loadtxt (outfile1)[0:,1] / (l*(l+1)/(2*np.pi)) 
			clTT2 = np.loadtxt (outfile2)[0:,1] / (l*(l+1)/(2*np.pi)) 
			clTE1 = np.loadtxt (outfile1)[0:,3] / (l*(l+1)/(2*np.pi)) 
			clTE2 = np.loadtxt (outfile2)[0:,3] / (l*(l+1)/(2*np.pi)) 
			clEE1 = np.loadtxt (outfile1)[0:,2] / (l*(l+1)/(2*np.pi)) 
			clEE2 = np.loadtxt (outfile2)[0:,2] / (l*(l+1)/(2*np.pi)) 
			dev_clTT = (clTT2-clTT1)/(2*self.params_list[j]*stepsize)
			dev_clTE = (clTE2-clTE1)/(2*self.params_list[j]*stepsize)
			dev_clEE = (clEE2-clEE1)/(2*self.params_list[j]*stepsize)
			dev_cl['TT'] = dev_clTT
			dev_cl['TE'] = dev_clTE
			dev_cl['EE'] = dev_clEE
			self.deriv_vec[param] = dev_cl

	def cov_matrix (self):
		""" Construct covariance matrix """

		cov = {}
		tag = "0"
		outfile = run_cmb (self.params_list, tag)
		l = np.loadtxt (outfile)[0:,0]
		self.l_list = l
		clTT = np.loadtxt (outfile)[0:,1] / (l*(l+1)/(2*np.pi)) 
		clTE = np.loadtxt (outfile)[0:,3] / (l*(l+1)/(2*np.pi)) 
		clEE = np.loadtxt (outfile)[0:,2] / (l*(l+1)/(2*np.pi)) 
		cov['11'] = 2 * clTT**2 / (2*l+1)
		
		cov['12'] = 2 * clTT*clTE / (2*l+1)
		cov['21'] = 2 * clTT*clTE / (2*l+1)
		cov['13'] = 2 * clTE**2	/ (2*l+1)
		cov['31'] = 2 * clTE**2	/ (2*l+1)
		cov['22'] = (clTT*clEE + clTE**2) / (2*l+1)
		cov['23'] = 2 * clTE*clEE / (2*l+1)
		cov['32'] = 2 * clTE*clEE / (2*l+1)
		cov['33'] = 2 * clEE**2 / (2*l+1)
		
		self.cov = cov
		
	def cmb_fisher_analysis (self):
		""" Do fisher analysis with results from deriv_vec (self) and cov_matrix (self) """
		
		F = np.zeros([len(self.fisher_params),len(self.fisher_params)])
		for m in range(len(self.fisher_params)):
			for n in range(len(self.fisher_params)):
				F_mn = 0
				param_m = self.fisher_params[m]
				param_n = self.fisher_params[n]
				vec_m = np.zeros(3)
				vec_n = np.zeros(3)
				inv_cov = np.zeros([3,3])
					
				for l in range(len(self.l_list)):
					for i in range(3):
						for j in range(3):
							if i == 0:
								vec_m[i] = self.deriv_vec[param_m]["TT"][l]
							elif i == 1:
								vec_m[i] = self.deriv_vec[param_m]["TE"][l]
							elif i == 2:
								vec_m[i] = self.deriv_vec[param_m]["EE"][l]
							if j == 0:
								vec_n[j] = self.deriv_vec[param_n]["TT"][l]
							elif j == 1:
								vec_n[j] = self.deriv_vec[param_n]["TE"][l]
							elif j == 2:
								vec_n[j] = self.deriv_vec[param_n]["EE"][l]
							
							inv_cov[i,j] = 1 / self.cov["{0}{1}".format(i+1,j+1)][l]
							F_mn += vec_m[i]*inv_cov[i,j]*vec_n[j]
				F[m,n] = F_mn
		return F




class fisher (object): 
	def __init__ (self):

		self.stepsize = [0.01, 0.005, 0.015]
		self.params_list = np.loadtxt (params_input)[0:,]
		self.fisher_params = ["h","b"]		# The order of parameters should be the same as in params_list
		
		self.z_m_list = [30, 50]
		self.w_list = [3.31685, 2.47355]

		self.deriv_vec = {}
		self.cov = {}

		self.l_list = None
	def cl21T_deriv_vec (self):
		""" Construct vector of derivative of cl21T w.r.t each parameter"""
		
		for j in range(len(self.fisher_params)):
			for i in [0]:#range(len(self.stepsize)):
				param = self.fisher_params[j]

				params_list_copy = self.params_list.copy ()
				params_list_copy[j] *= (1 - self.stepsize[i])
				tag = param + "_{}1".format(i)
				run_21cm (params_list_copy, tag)
				Cl1 = set_cl_21 (tag)
				if j == 0:
					self.l_list = Cl1.l_list

				params_list_copy = self.params_list.copy ()
				params_list_copy[j] *= (1 + self.stepsize[i])
				tag = param + "_{}2".format(i)
				run_21cm (params_list_copy, tag)
				Cl2 = set_cl_21 (tag)
				
				dev_cl = {}
				for k in range(len(self.z_m_list)):
					cl1_zm = Cl1.cl21T (self.z_m_list[k], self.w_list[k])
					cl2_zm = Cl2.cl21T (self.z_m_list[k], self.w_list[k])
					dev_cl_zm = (cl2_zm-cl1_zm)/(2*self.params_list[j]*self.stepsize[i])
					dev_cl['{0}'.format (self.z_m_list[k])] = dev_cl_zm
			self.deriv_vec[param] = dev_cl

	def cov_matrix (self):
		""" Construct covariance matrix """

		cl21T = {}
		cl21 = {}
		tag = "0"
		run_21cm (self.params_list, tag)
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
		clTT = np.interp (self.l_list, l, cl)
		for i in range(len(self.z_m_list)):
			for j in range(len(self.z_m_list)):
				element_ij = (cl21["{0}{1}".format(i,j)] * clTT + cl21T["{0}".format(i)]*cl21T["{0}".format(j)]) / (2*self.l_list+1)
				self.cov["{0}{1}".format(i,j)] = element_ij
		
	def fisher_analysis (self):
		""" Do fisher analysis with results from deriv_vec (self) and cov_matrix (self) """
		
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
		""" Convergence test for the derivative of cl21T """
		
		if len(self.deriv_vec) == 0:
			print ("ERROR : Do cl21T_deriv_vec (self) first")
			return None

		deriv_cl_dic = {}
		for j in range(len(self.fisher_params)):
			param = self.fisher_params[j]
			
			dev_cl = {}
			for i in [1,2]:
				params_list_copy = self.params_list.copy ()
				params_list_copy[j] *= (1 - self.stepsize[i])
				tag = param + "_{}1".format(i)
				run_21cm (params_list_copy, tag)
				Cl1 = set_cl_21 (tag)
				
				params_list_copy = self.params_list.copy ()
				params_list_copy[j] *= (1 + self.stepsize[i])
				tag = param + "_{}2".format(i)
				run_21cm (params_list_copy, tag)
				Cl2 = set_cl_21 (tag)
				
				for k in range(len(self.z_m_list)):
					cl1_zm = Cl1.cl21T (self.z_m_list[k], self.w_list[k])
					cl2_zm = Cl2.cl21T (self.z_m_list[k], self.w_list[k])
					dev_cl_zm = (cl2_zm-cl1_zm)/(2*self.params_list[j]*self.stepsize[i])
					dev_cl['{0}_{1}'.format (self.z_m_list[k],i)] = dev_cl_zm
			
			deriv_cl_dic[param] = dev_cl
	
		print (deriv_cl_dic)
		n = 1
		for j in range(len(self.fisher_params)):
			param = self.fisher_params[j]
			for k in range(len(self.z_m_list)):
				dcl0 = self.deriv_vec[param]['{0}'.format (self.z_m_list[k])]
				dcl1 = deriv_cl_dic[param]['{0}_{1}'.format (self.z_m_list[k],1)]
				dcl2 = deriv_cl_dic[param]['{0}_{1}'.format (self.z_m_list[k],2)]
				conv_test_b1 = np.log10(abs((dcl0-dcl1)/dcl0))
				conv_test_b2 = np.log10(abs((dcl0-dcl2)/dcl0))
				
				plt.figure(n)
				plt.plot (self.l_list, conv_test_b1, label = r'stepsize$\times$0.5')
				plt.plot (self.l_list, conv_test_b2, label = r'stepsize$\times$1.5')
				plt.savefig (path_result + "/conv/conv_{0}_{1}.pdf".format (param, self.z_m_list[k]))
				n += 1
		
		print ("Convergence test done")
		return None	

# Convergence Test	
"""
F = fisher ()
F.cl21T_deriv_vec ()
F.convergence_test ()
"""
# Fisher Analysis
"""
F = fisher ()
F.cl21T_deriv_vec ()
F.cov_matrix ()
fisher_matrix = F.fisher_analysis ()
print (fisher_matrix)
print (inv(fisher_matrix))
print (np.dot (fisher_matrix, inv(fisher_matrix)))
"""
F = prior_cmb ()
F.cmb_deriv_vec ()
F.cov_matrix ()
fisher_matrix = F.cmb_fisher_analysis ()
print (fisher_matrix)
print (inv(fisher_matrix))
print (np.dot (fisher_matrix, inv(fisher_matrix)))
