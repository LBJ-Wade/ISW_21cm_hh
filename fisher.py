from run import *
from numpy.linalg import inv
from cl_21 import *

class fisher (object): 
	def __init__ (self, comm, rank, size, Yp_BBN = True):
		self.comm = comm
		self.rank = rank
		self.size = size
		
		# If Yp_BBN = True, the helium fraction Yp is fixed by BBN constraint.
		self.Yp_BBN = Yp_BBN	
		if self.Yp_BBN == True:
			self.params_path = path_params_Yp_BBN
			self.params_input = self.params_path + "/params_nonu.dat"
			self.params_list = np.loadtxt (self.params_input)[0:,]
			self.fisher_params = ['c','b','theta','tau', 'A_s','n_s','m_nu','Neff']
			self.stepsize = [0.0030, 8.0e-4, 5.0e-3, 0.02, 0.045, 0.01, 0.02, 0.08]
			
			self.path_result = path_result_Yp_BBN	
		else:
			self.params_path = path_params_Yp
			self.params_input = self.params_path + "/params_nonu.dat"
			self.params_list = np.loadtxt (self.params_input)[0:,]
			self.fisher_params = ['c','b','theta','tau', 'A_s','n_s','m_nu','Neff', 'Yp']
			self.stepsize = [0.0030, 8.0e-4, 5.0e-5, 0.02, 0.045, 0.01, 0.02, 0.08, 0.0048]
			self.path_result = path_result_Yp

		# Redshifts where the code calculate correlation functions
		self.z_m_list = [30, 50, 75, 100, 125, 150, 175, 200]
		
		# Widths of window functions.
		# Window function is defied in selection_function21cm.py.
		self.w_list = [3.1891, 2.4779, 2.0211, 1.7456, 1.5560, 1.4151, 1.3050, 1.2156] 

		self.deriv_vec = {}
		self.cov = {}

		self.l_list = None

	def cl21T_deriv_vec (self):
		""" Construct vector of partial derivative of cl21T w.r.t each parameter"""
		print ('(Rank:{0}) fisher.py : Start cl21T_deriv_vec'.format(self.rank))
		
		for j in range(len(self.fisher_params)):
			param = self.fisher_params[j]
			stepsize = self.stepsize[j]
			
			# Run CLASS using run.py to get data
			# Data is saved automatically in data directory (See run.py for detail)
			infile1 = self.params_path + "/params_" + param + "1.dat"
			params_list_copy = np.loadtxt (infile1)[0:,]
			tag1 = param + "_01"
			run_21cm (self.comm, self.rank, self.size, params_list_copy, infile1, tag1, self.Yp_BBN)
			self.comm.Barrier ()

			infile2 = self.params_path + "/params_" + param + "2.dat"
			params_list_copy = np.loadtxt (infile2)[0:,] 
			tag2 = param + "_02"
			run_21cm (self.comm, self.rank, self.size, params_list_copy, infile2, tag2, self.Yp_BBN)
			self.comm.Barrier ()

			# Construct object of class cl_21 of cl_21.py for two parameter sets.
			Cl1 = set_cl_21 (self.comm, self.rank, self.size, tag1, self.Yp_BBN)
			Cl2 = set_cl_21 (self.comm, self.rank, self.size, tag2, self.Yp_BBN)
			if j == 0:
				l_list = Cl1.l_list
				self.l_list = l_list
			
			# Calculate the cross-correlation functions and the partial derivative of it.
			# The partial derivative is saved in result_Yp_BBN directory. 
			dev_cl = {}
			dev_list = []
			for k in range(len(self.z_m_list)):
				cl1_zm = Cl1.cl21T (self.z_m_list[k], self.w_list[k])
				cl2_zm = Cl2.cl21T (self.z_m_list[k], self.w_list[k])
				if self.rank == 0:
					dev_cl_zm = (cl2_zm-cl1_zm)/(2*stepsize)
					dev_cl['{0}'.format (self.z_m_list[k])] = dev_cl_zm *2.7255
					dev_list.append (dev_cl_zm)
			
			out_zm = self.path_result + '/cl21T_deriv_'+param+'.txt'
			if self.rank == 0:
				data = np.column_stack((l_list, dev_list[0], dev_list[1], dev_list[2], dev_list[3], dev_list[4], dev_list[5], dev_list[6], dev_list[7]))
				np.savetxt(out_zm, data, fmt = '%1.6e')
			
			"""	
			out_zm = self.path_result + '/cl21T_deriv_'+param+'.txt'
			dev_cl = {}
			for k in range(len(self.z_m_list)):
				dev_cl_zm = np.loadtxt(out_zm)[0:,k+1] *2.7255
				dev_cl['{0}'.format (self.z_m_list[k])] = dev_cl_zm
			"""
			self.deriv_vec[param] = dev_cl
		if self.rank == 0:
			os.system ("rm data/*")				
	def cov_matrix (self):
		""" Construct covariance matrix """
		print ('\n')
		print ('fisher.py : Start cov_matrix')
		f_sky = 0.4
		cl21T = {}
		cl21 = {}
		tag = "0"
		
		# Run CLASS with fiducial cosmological parameters using run.py.
		run_21cm (self.comm, self.rank, self.size,self.params_list, self.params_input, tag, self.Yp_BBN)
		self.comm.Barrier ()
		
		# Construct the object of class cl_21 of cl_21.py.
		Cl = set_cl_21 (self.comm, self.rank, self.size, tag, self.Yp_BBN)
		l_list = Cl.l_list	
		cl21T_list = []
		
		# Calculate the cross-correlation functions and save them in result_Yp_BBN directory.
		for i in range(len(self.z_m_list)):
			cl_zm = Cl.cl21T (self.z_m_list[i], self.w_list[i])
			if self.rank == 0:
				cl21T["{0}".format(i,i)] = cl_zm *2.7255
				cl21T_list.append (cl_zm)
		out_zm = self.path_result + '/cl21T_0.txt'
		if self.rank == 0:
			data = np.column_stack((l_list, cl21T_list[0], cl21T_list[1], cl21T_list[2], cl21T_list[3], cl21T_list[4], cl21T_list[5], cl21T_list[6], cl21T_list[7]))
			np.savetxt(out_zm, data, fmt = '%1.6e')
		
		"""	
		out_zm = self.path_result + '/cl21T_0.txt'
		for i in range(len(self.z_m_list)):
			cl_zm = np.loadtxt (out_zm)[0:,i+1] *2.7255
			cl21T["{0}".format(i,i)] = cl_zm
		"""

		# Calculate the auto-correlation functions and save them in result_Yp_BBN directory.
		for i in range(len(self.z_m_list)):
			for j in range(len(self.z_m_list)):
				if not j < i:
					zm = [self.z_m_list[i], self.z_m_list[j]]
					w = [self.w_list[i], self.w_list[j]]
									
					cl_zmzl = Cl.cl21 (zm, w)
					if self.rank == 0:
						if i == j:
							cl21["{0}{1}".format(i,j)] = cl_zmzl
						else:
							cl21["{0}{1}".format(i,j)] = cl_zmzl
							cl21["{0}{1}".format(j,i)] = cl_zmzl
						out_zm = self.path_result + '/cl21zmzl_{0}{1}.txt'.format(i,j)
						data = np.column_stack((l_list, cl_zmzl))
						np.savetxt(out_zm, data, fmt = '%1.6e')
					
					"""	
					out_zm = self.path_result + '/cl21zmzl_{0}{1}.txt'.format(i,j)
					cl_zmzl = np.loadtxt(out_zm)[0:,1]
					if i == j:
						cl21["{0}{1}".format(i,j)] = cl_zmzl
					else:
						cl21["{0}{1}".format(i,j)] = cl_zmzl
						cl21["{0}{1}".format(j,i)] = cl_zmzl
					"""
		cl_out = self.path_result + "/cl_" + tag + ".dat"
		l = np.loadtxt(cl_out)[0:,0]
		aa = 10**-12 * 2*np.pi / (l*(l+1))
		cl = np.loadtxt(cl_out)[0:,1]*aa
		clTT = np.interp (l_list, l, cl)
		clTT_N = (2*10**-6*0.000290888)**2 *np.e**(l_list*(l_list+1)*(0.000290888**2)/(8*np.log(2)))
		clTT += clTT_N
	
		if self.rank == 0:
			for i in range(len(self.z_m_list)):
				for j in range(len(self.z_m_list)):
					element_ij = (cl21["{0}{1}".format(i,j)] * clTT + cl21T["{0}".format(i)]*cl21T["{0}".format(j)]) / (2*l_list+1)/f_sky
					self.cov["{0}{1}".format(i,j)] = element_ij
		
	def fisher_analysis (self):
		""" Do fisher analysis with results from deriv_vec (self) and cov_matrix (self) """
		print ('\n')
		print ('fisher.py : Start cl21T fisher_analysis')
		
		F = np.zeros([len(self.fisher_params),len(self.fisher_params)])
		for m in range(len(self.fisher_params)):
			for n in range(len(self.fisher_params)):
				F_mn = 0
				param_m = self.fisher_params[m]
				param_n = self.fisher_params[n]
				#for l in range(28,2999):
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

