from run import *
from numpy.linalg import inv
from cl_21 import *

class prior_cmb (object):		# Need to add cl^dd
	def __init__ (self):
		self.stepsize = [0.0030, 8.0e-4, 5.0e-5, 0.02, 0.1e-9, 0.01, 0.02/3, [0.7117357, 0.721146] ]
		self.params_list = np.loadtxt (params_input)[0:,]
		self.fisher_params = ['c','b','theta','tau', 'A_s','n_s','m_nu','Neff']
		self.deriv_vec = {}
		self.cov = {}
		self.l_list = None
	def cmb_deriv_vec (self):
		print ('Start cmb_deriv_vec')

		for j in range(len(self.fisher_params)):
			param = self.fisher_params[j]
			stepsize = self.stepsize[j]

			"""
			params_list_copy = self.params_list.copy ()
			if param == 'Neff':
				params_list_copy[j] = stepsize[0]
			else:
				params_list_copy[j] -= stepsize 
			tag = param + "1"
			outfile1 = run_cmb (params_list_copy, tag)

			params_list_copy = self.params_list.copy ()
			if param == 'Neff':
				params_list_copy[j] = stepsize[1]
			else:
				params_list_copy[j] += stepsize 
			tag = param + "2"
			outfile2 = run_cmb (params_list_copy, tag)
			"""
			outfile1 = path_data + "/cl_" + param+"1" + ".dat"
			outfile2 = path_data + "/cl_" + param+"2" + ".dat"
			dev_cl = {}
			l = np.loadtxt (outfile1)[28:,0]
			clTT1 = np.loadtxt (outfile1)[28:,1] / (l*(l+1)/(2*np.pi)) 
			clTT1[2970:] = 0
			clTT2 = np.loadtxt (outfile2)[28:,1] / (l*(l+1)/(2*np.pi)) 
			clTT2[2970:] = 0
			clTE1 = np.loadtxt (outfile1)[28:,4] / (l*(l+1)/(2*np.pi)) 
			clTE2 = np.loadtxt (outfile2)[28:,4] / (l*(l+1)/(2*np.pi)) 
			clEE1 = np.loadtxt (outfile1)[28:,2] / (l*(l+1)/(2*np.pi)) 
			clEE2 = np.loadtxt (outfile2)[28:,2] / (l*(l+1)/(2*np.pi)) 
			cldd1 = np.loadtxt (outfile1)[28:,5] / (l*(l+1)/(2*np.pi))
			cldd2 = np.loadtxt (outfile2)[28:,5] / (l*(l+1)/(2*np.pi))
			
			if param == 'Neff':
				dev_clTT = (clTT2-clTT1)/(2*0.08)
				dev_clTE = (clTE2-clTE1)/(2*0.08)
				dev_clEE = (clEE2-clEE1)/(2*0.08)
				dev_cldd = (cldd2-cldd1)/(2*0.08)
			else param == 'm_nu':
				dev_clTT = (clTT2-clTT1)/(2*stepsize*3)
				dev_clTE = (clTE2-clTE1)/(2*stepsize*3)
				dev_clEE = (clEE2-clEE1)/(2*stepsize*3)
				dev_cldd = (cldd2-cldd1)/(2*stepsize*3)
			else:	
				dev_clTT = (clTT2-clTT1)/(2*stepsize)
				dev_clTE = (clTE2-clTE1)/(2*stepsize)
				dev_clEE = (clEE2-clEE1)/(2*stepsize)
				dev_cldd = (cldd2-cldd1)/(2*stepsize)
			dev_cl['TT'] = dev_clTT
			dev_cl['TE'] = dev_clTE
			dev_cl['EE'] = dev_clEE
			dev_cl['dd'] = dev_cldd
			self.deriv_vec[param] = dev_cl

	def cov_matrix (self):
		""" Construct covariance matrix """
		print ('Start cov_matrix')
		
		f_sky = 0.7	
		cov = {}
		tag = "0"
		#outfile = run_cmb (self.params_list, tag)
		outfile = path_data + "/cl_" + tag + ".dat"
		l = np.loadtxt (outfile)[28:,0]
		clTT_N = (2*0.000290888)**2 *np.e**(l*(l+1)*0.000290888**2/(8*np.log(2)))
		clEE_N = 2*clTT_N
		self.l_list = l
		clTT = np.loadtxt (outfile)[28:,1] / (l*(l+1)/(2*np.pi)) 
		clTT += clTT_N
		clTE = np.loadtxt (outfile)[28:,4] / (l*(l+1)/(2*np.pi)) 
		clEE = np.loadtxt (outfile)[28:,2] / (l*(l+1)/(2*np.pi)) 
		clEE += clEE_N
		cldd = np.loadtxt (outfile)[28:,5] / (l*(l+1)/(2*np.pi))
		clTd = np.loadtxt (outfile)[28:,6] / (l*(l+1)/(2*np.pi))
		clEd = np.loadtxt (outfile)[28:,7] / (l*(l+1)/(2*np.pi))
		
		cov['11'] = f_sky * 2 * clTT**2 / (2*l+1)
		cov['22'] = f_sky * (clTT*clEE + clTE**2) / (2*l+1)
		cov['33'] = f_sky * 2 * clEE**2 / (2*l+1)
		cov['44'] = f_sky * 2 * cldd**2 / (2*l+1)

		cov['12'] = f_sky * 2 * clTT*clTE / (2*l+1); cov['21'] = f_sky * 2 * clTT*clTE / (2*l+1)
		cov['13'] = f_sky * 2 * clTE**2	/ (2*l+1); cov['31'] = f_sky * 2 * clTE**2	/ (2*l+1)
		cov['14'] = f_sky * 2 * clTd**2 / (2*l+1); cov['41'] = f_sky * 2 * clTd**2 / (2*l+1)
		
		cov['23'] = f_sky * 2 * clTE*clEE / (2*l+1); cov['32'] = f_sky * 2 * clTE*clEE / (2*l+1)
		cov['24'] = f_sky * 2 * clTd*clEd / (2*l+1); cov['42'] = f_sky * 2 * clTd*clEd / (2*l+1)
		
		cov['34'] = f_sky * 2 * clEd**2 / (2*l+1); cov['43'] = f_sky * 2 * clEd**2 / (2*l+1)

		self.cov = cov
		
	def cmb_fisher_analysis (self):
		""" Do fisher analysis with results from deriv_vec (self) and cov_matrix (self) """
		print ('Start cmb_fisher_analysis')
			
		F = np.zeros([len(self.fisher_params),len(self.fisher_params)])
		for m in range(len(self.fisher_params)):
			for n in range(len(self.fisher_params)):
				F_mn = 0
				param_m = self.fisher_params[m]
				param_n = self.fisher_params[n]
					
				#for l in range(len(self.l_list)):
				for l in range(2971):
					vec_m = np.zeros(4)
					vec_n = np.zeros(4)
					inv_cov = np.zeros([4,4])
					for i in range(4):
						for j in range(4):
							if i == 0:
								vec_m[i] = self.deriv_vec[param_m]["TT"][l]
							elif i == 1:
								vec_m[i] = self.deriv_vec[param_m]["TE"][l]
							elif i == 2:
								vec_m[i] = self.deriv_vec[param_m]["EE"][l]
							elif i == 3:
								vec_m[i] = self.deriv_vec[param_m]["dd"][l]
							
							if j == 0:
								vec_n[j] = self.deriv_vec[param_n]["TT"][l]
							elif j == 1:
								vec_n[j] = self.deriv_vec[param_n]["TE"][l]
							elif j == 2:
								vec_n[j] = self.deriv_vec[param_n]["EE"][l]
							elif j == 3:
								vec_n[j] = self.deriv_vec[param_n]["dd"][l]
							
							inv_cov[i,j] = self.cov["{0}{1}".format(i+1,j+1)][l]
					inv_cov = inv(inv_cov)
					F_mn += np.dot(vec_m, np.dot(inv_cov, vec_n))
				for l in range(2971,len(self.l_list)):
					vec_m = np.zeros(3)
					vec_n = np.zeros(3)
					inv_cov = np.zeros([3,3])
					
					for i in range(3):
						for j in range(3):
							if i == 0:
								vec_m[i] = self.deriv_vec[param_m]["TE"][l]
							elif i == 1:
								vec_m[i] = self.deriv_vec[param_m]["EE"][l]
							elif i == 2:
								vec_m[i] = self.deriv_vec[param_m]["dd"][l]
							
							if j == 0:
								vec_n[j] = self.deriv_vec[param_n]["TE"][l]
							elif j == 1:
								vec_n[j] = self.deriv_vec[param_n]["EE"][l]
							elif j == 2:
								vec_n[j] = self.deriv_vec[param_n]["dd"][l]
				
							inv_cov[i,j] = self.cov["{0}{1}".format(i+2,j+2)][l]
					inv_cov = inv(inv_cov)
					F_mn += np.dot(vec_m, np.dot(inv_cov, vec_n))
				F[m,n] = F_mn
		return F




class fisher (object): 
	def __init__ (self):
		
		self.params_list = np.loadtxt (params_input)[0:,]
		
		self.stepsize = [0.0030, 8.0e-4, 5.0e-5, 0.02, 0.1e-9, 0.01, 0.02/3, [0.7117357, 0.721146] ]
		self.fisher_params = ['c','b','theta','tau', 'A_s','n_s','m_nu','Neff']
		
		#self.z_m_list = [30, 50]
		self.z_m_list = [30, 50,75,100,125,150,175,200]
		#self.w_list = [3.31685, 2.47355]
		self.w_list = [3.31685, 2.47355 ,1.93014, 1.60434, 1.38246, 1.21989, 1.09467, 1]
		self.deriv_vec = {}
		self.cov = {}

		self.l_list = np.arange(2,5001)
	def cl21T_deriv_vec (self):
		""" Construct vector of derivative of cl21T w.r.t each parameter"""
		
		for j in range(len(self.fisher_params)):
			
			param = self.fisher_params[j]
			stepsize = self.stepsize[j]
			params_list_copy = self.params_list.copy ()
			if param == 'Neff':
				params_list_copy[j] = stepsize[0]
			else:
				params_list_copy[j] -= stepsize
			tag = param + "_01"
			run_21cm (params_list_copy, tag)
			Cl1 = set_cl_21 (tag)
			if j == 0:
				l_list = Cl1.l_list

			params_list_copy = self.params_list.copy ()
			if param == 'Neff':
				params_list_copy[j] = stepsize[1]
			else:
				params_list_copy[j] += stepsize 
			tag = param + "_02"
			run_21cm (params_list_copy, tag)
			Cl2 = set_cl_21 (tag)
			
			
			dev_cl = {}
			dev_list = []
			for k in range(len(self.z_m_list)):
				cl1_zm = Cl1.cl21T (self.z_m_list[k], self.w_list[k])
				cl1_zm = np.interp (self.l_list, l_list, cl1_zm)
				cl2_zm = Cl2.cl21T (self.z_m_list[k], self.w_list[k])
				cl2_zm = np.interp (self.l_list, l_list, cl2_zm)
				if param == 'Neff':
					dev_cl_zm = (cl2_zm-cl1_zm)/(2*0.08)
				else:
					dev_cl_zm = (cl2_zm-cl1_zm)/(2*stepsize)
				dev_cl['{0}'.format (self.z_m_list[k])] = dev_cl_zm
				dev_list.append (dev_cl_zm)
			out_zm = path_result + '/cl21T_deriv_'+param+'.txt'
			data = np.column_stack((self.l_list, dev_list[0], dev_list[1], dev_list[2], dev_list[3], dev_list[4], dev_list[5], dev_list[6], dev_list[7]))
			np.savetxt(out_zm, data, fmt = '%1.6e')
			
			self.deriv_vec[param] = dev_cl
			
	def cov_matrix (self):
		""" Construct covariance matrix """

		cl21T = {}
		cl21 = {}
		tag = "0"
		run_21cm (self.params_list, tag)
		Cl = set_cl_21 (tag)
		l_list = Cl.l_list	
		cl21T_list = []
		
		for i in range(len(self.z_m_list)):
			cl_zm = Cl.cl21T (self.z_m_list[i], self.w_list[i])
			cl_zm = np.interp (self.l_list, l_list, cl_zm)
			cl21T["{0}".format(i,i)] = cl_zm
			cl21T_list.append (cl_zm)
		out_zm = path_result + '/cl21T_0.txt'
		data = np.column_stack((self.l_list, cl21T_list[0], cl21T_list[1], cl21T_list[2], cl21T_list[3], cl21T_list[4], cl21T_list[5], cl21T_list[6], cl21T_list[7]))
		np.savetxt(out_zm, data, fmt = '%1.6e')
		
		for i in range(len(self.z_m_list)):
			for j in range(len(self.z_m_list)):
				if not j < i:
					zm = [self.z_m_list[i], self.z_m_list[j]]
					w = [self.w_list[i], self.w_list[j]]
					cl_zmzl = Cl.cl21 (zm, w)
					cl_zmzl = np.interp (self.l_list, l_list, cl_zmzl)
					if i == j:
						cl21["{0}{1}".format(i,j)] = cl_zmzl
					else:
						cl21["{0}{1}".format(i,j)] = cl_zmzl
						cl21["{0}{1}".format(j,i)] = cl_zmzl
					out_zm = path_result + '/cl21zmzl_{0}{1}.txt'.format(i,j)
					data = np.column_stack((self.l_list, cl_zmzl))
					np.savetxt(out_zm, data, fmt = '%1.6e')
		
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
							inv_cov[i,j] = self.cov["{0}{1}".format(i,j)][l]
					inv_cov = inv(inv_cov)
					F_mn += np.dot (vec_m, np.dot (inv_cov ,vec_n))
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
'''
# Fisher Analysis
F = fisher ()
F.cl21T_deriv_vec ()
F.cov_matrix ()
fisher_matrix = F.fisher_analysis ()
inv_fisher = inv(fisher_matrix)
print (fisher_matrix)
print (inv(fisher_matrix))
print (np.dot (fisher_matrix, inv(fisher_matrix)))

'''
# Fisher CMB
F = prior_cmb ()
F.cmb_deriv_vec ()
F.cov_matrix ()
fisher_matrix = F.cmb_fisher_analysis ()
inv_fisher = inv(fisher_matrix)
#print (fisher_matrix)
#print (inv(fisher_matrix))
#print (np.dot (fisher_matrix, inv(fisher_matrix)))

sigma = []
for i in range(len(fisher_matrix)):
	sigma.append (inv_fisher[i,i])
sigma = np.array(sigma)
data = np.column_stack((sigma))

#np.savetxt('sigma_cl21T.txt', data, fmt = '%1.6e')
#data = np.column_stack((fisher_matrix[0],fisher_matrix[1],fisher_matrix[2],fisher_matrix[3],fisher_matrix[4],fisher_matrix[5],fisher_matrix[6],fisher_matrix[7]))
#np.savetxt('fisher_matrix_cl21T.txt', data, fmt = '%1.6e')
np.savetxt('sigma.txt', data, fmt = '%1.6e')
data = np.column_stack((fisher_matrix[0],fisher_matrix[1],fisher_matrix[2],fisher_matrix[3],fisher_matrix[4],fisher_matrix[5],fisher_matrix[6]))
np.savetxt('fisher_matrix.txt', data, fmt = '%1.6e')

