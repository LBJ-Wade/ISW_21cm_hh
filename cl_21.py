import numpy as np
from scipy.interpolate import interp2d
from scipy.integrate import simps
import selection_function21cm as sf
import matplotlib.pyplot as plt
import scipy.special
import sys
import matplotlib.ticker

def set_cl_21 (tag):		# At the end, it would take z_m_list, w_list as arguments
	""" Construct object of class cl_21 """
	
	infile = 'data/file_names_{0}.txt'.format (tag)
	file_names = np.genfromtxt(infile, dtype="str")[0:]
	params_input = file_names[0]
	infile_syn = file_names[1]
	infile_new = file_names[2]
	infile_HYREC = file_names[3]
	outfile = file_names[4]
	
	params_list = np.loadtxt (params_input)[0:,]
	print infile_HYREC
	Cl = cl_21 (params_list, infile_HYREC, infile_syn, infile_new)
	Cl.test_run () 
	Cl.c_z ()
	return Cl
	

class cl_21 (object):
	def __init__ (self, params_list, infile_HYREC, infile_syn, infile_new = None):

		self.params_list = params_list
		self.infile_syn = infile_syn
		self.infile_new = infile_new
		self.infile_HYREC = infile_HYREC
		#self.infile_21 = infile_21

		z = np.loadtxt(self.infile_syn)[0:,0]
		print ('z data', len(z))
		k = np.loadtxt(self.infile_syn)[0:,1]
		print ('k data', len(k))
		self.zlist2 = np.array(sorted(set(z)))
		self.klist2 = np.array(sorted(set(k)))
		self.number_of_z2 = len(self.zlist2)
		self.number_of_k2 = len(self.klist2)
		self.hubble_class = np.loadtxt(self.infile_syn)[0:self.number_of_z2,4][::-1]
		print ('hubble data', len(self.hubble_class))
		self.baryon = -np.loadtxt(self.infile_syn)[0:,2]
		self.baryon_dot = -np.loadtxt(self.infile_syn)[0:,3]

		As = self.params_list[5]
		n_s = self.params_list[6]
		k_pivot = 0.05
		self.k_list = np.logspace (-4,5, 5000)
		self.P_phi = As * (self.k_list/k_pivot)**(n_s-1) * 2*np.pi**2 / self.k_list**3

		a_0 = 10**-2/4
		n = 10000
		scale_factor = np.logspace (np.log10(a_0), 0, n)
		scale_factor_reverse = scale_factor[::-1]
		self.redshift2 = 1/scale_factor_reverse - 1
		hubble_class2 = np.interp (self.redshift2, self.zlist2, self.hubble_class)

		chi_class = []
		for i in range(len(hubble_class2)):
			chi_class.append (simps (1/hubble_class2[:i+1], self.redshift2[:i+1]))
		chi_class = np.array (chi_class)
		self.chi_class = np.interp (self.zlist2, self.redshift2, chi_class)
		
		# HYREC 
		self.z_HYREC = np.loadtxt(self.infile_HYREC)[0:,0]
		self.x_HYREC = np.loadtxt(self.infile_HYREC)[0:,1]
		Tm_Tr = np.loadtxt(self.infile_HYREC)[0:,2]
		T_cmb = 2.7255
		self.Tr = T_cmb*(1+self.z_HYREC)
		self.Tm = self.Tr * Tm_Tr
		self.T_T = None
		self.T_H = None
		self.T_b = None

		self.c = 299792458
		self.Mpc_to_m = 3.0857*10**22
		self.Yp = 0.245
		self.Omega_b = self.params_list[1]
		self.eV_to_m_inv = 5076142.131979696
		self.h = self.params_list[0]
		self.rho_cr = 8.056*10**-11 * self.h**2 # eV^4
		self.mp = 938.2720813*10**6  #eV
		self.me = 0.5109989461*10**6	#eV
		self.sigma_T = 6.6524587158 * 10 **-29
		self.J_to_eV = 6.2415093433*10**18
		self.k_B = 8.613303*10**-5 # eV/K
		self.wavelength = 0.21 #m
		self.E10 = 0.068															# K
		self.A10 = 2.85*10**-15 / self.c											# /m	
		self.B10 = self.A10*(1+1/(np.e**(self.E10/self.Tr)-1))							# /m
		
		l_list = np.logspace(np.log10(2), np.log10(5000), 1000)
		#l_list = np.logspace(np.log10(2), np.log10(10), 2)
		for i in range(len(l_list)):
			l_list[i] = int(l_list[i])
		l_list = sorted(set(l_list))
		l_list[-1] += 1
		self.l_list = np.array (l_list)


	def test_run (self):
		""" Calculate coefficients of 21 cm fluctuations (linear terms) """

		zlist = self.zlist2.copy ()
		hubble_class = self.hubble_class.copy ()
	
		z = self.z_HYREC.copy ()
		x = self.x_HYREC.copy ()

		hubble = np.interp (z[::-1], zlist, hubble_class)[::-1]
		n_H = (1-self.Yp)*self.rho_cr*self.Omega_b/self.mp*(1+z)**3 		# eV^3
		n_H *= self.eV_to_m_inv**3											# /m^3
		kappa = 3.1*10**-11 *self.Tm**0.357 * np.e**(-32/self.Tm) * 10**-6 			# m^3/s
		kappa /= self.c														# m^2
		wavelength = self.wavelength	
		C10 = n_H*kappa														# /m
		C01 = 3*np.e**(-self.E10/self.Tm)*C10											# /m
		I_nu = 2*self.k_B*self.Tr/wavelength**2									# eV/m^2
		I_nu *= self.eV_to_m_inv											# /m^3
	
		hubble = hubble/self.Mpc_to_m
		Ts = self.Tr + (self.Tm-self.Tr)*C10/(C10+self.A10*self.Tm/self.E10)
	
		a = kappa*n_H
		b = self.A10*self.Tr/self.E10
		D = (1+z)**-1 * 3*self.E10/(32*np.pi)*1/4*self.Yp*wavelength**3*self.A10*n_H**2*kappa/(hubble*self.Tm*(kappa*n_H+self.A10*self.Tr/self.E10))
		C = (1+z)**-1 * 3*self.E10/(32*np.pi)*(1-x)*wavelength**3*self.A10
		T_b = -C*self.E10*kappa*n_H**2*(self.Tr-self.Tm)/(hubble*(self.E10*kappa*n_H + self.A10*self.Tr)*self.Tm)
		T_H = -C*self.E10*kappa*n_H**2*(self.E10*kappa*n_H+2*self.A10*self.Tr)*(self.Tr-self.Tm)/(hubble*(self.E10*kappa*n_H+self.A10*self.Tr)**2*self.Tm)
		T_T = (0.357*C*self.E10*kappa*n_H**2*(-89.6359*self.A10*self.Tr**2+89.6359*self.A10*self.Tr*self.Tm+2.80112*self.E10*kappa*n_H*self.Tr*self.Tm+1.80112*self.A10*self.Tr**2*self.Tm+self.A10*self.Tr*self.Tm**2))/(hubble*(self.E10*kappa*n_H+self.A10*self.Tr)**2*self.Tm**2)
		T_HH = - self.A10**2*C*self.E10*kappa*n_H**2*self.Tr**2*(self.Tr-self.Tm)/(hubble*(self.E10*kappa*n_H+self.A10*self.Tr)**3*self.Tm)
		T_HT = (0.714*C*self.E10*kappa*n_H**2*(-89.6359*self.A10**2*self.Tr**3+1.40056*self.E10**2*kappa**2*n_H**2*self.Tr*self.Tm+89.6359*self.A10**2*self.Tr**2*self.Tm+4.20168*self.A10*self.E10*kappa*n_H*self.Tr**2*self.Tm+1.80112*self.A10**2*self.Tr**3*self.Tm+self.A10**2*self.Tr**2*self.Tm**2))/(hubble*(self.E10*kappa*n_H+self.A10*self.Tr)**3*self.Tm**2)
		T_TT = (C*self.E10*kappa*n_H**2*(512*self.A10*self.E10*kappa*n_H*self.Tr**2-512*self.A10**2*self.Tr**3-512*self.A10*self.E10*kappa*n_H*self.Tr*self.Tm+512*self.A10**2*self.Tr**2*self.Tm+75.424*self.A10*self.E10*kappa*n_H*self.Tr**2*self.Tm+52.576*self.A10**2*self.Tr**3*self.Tm-43.424*self.A10*self.E10*kappa*n_H*self.Tr*self.Tm**2-self.E10**2*kappa**2*n_H**2*self.Tr*self.Tm**2-20.576*self.A10**2*self.Tr**2*self.Tm**2-1.40078*self.A10*self.E10*kappa*n_H*self.Tr**2*self.Tm**2-0.528225*self.A10**2*self.Tr**3*self.Tm**2-1.11022*10**-16*self.E10**2*kappa**2*n_H**2*self.Tm**3-0.242225*self.A10*self.E10*kappa*n_H*self.Tr*self.Tm**3-0.114776*self.A10**2*self.Tr**2*self.Tm**3))/(hubble*(self.E10*kappa*n_H+self.A10*self.Tr)**3*self.Tm**3)

		self.T_T = T_T
		self.T_H = T_H
		self.T_b = T_b
		
		
		return z, T_T, T_H, T_b

	def c_z (self):
		""" Calculate C1(z) which is defiend as T_{Tgas} = C1(z) T_{b} """

		hubble_class = self.hubble_class[::-1].copy ()

		z = self.z_HYREC.copy ()
		x = self.x_HYREC.copy ()

		new_z = np.linspace(1000,0,100000)
		x = np.interp (new_z[::-1], z[::-1], x[::-1])[::-1]
		hubble = np.interp (new_z[::-1], self.zlist2, hubble_class[::-1])[::-1]
		Tr = np.interp (new_z[::-1], z[::-1], self.Tr[::-1])[::-1]
		Tm = np.interp (new_z[::-1], z[::-1], self.Tm[::-1])[::-1]
		a_r = 4*5.670373 * 10**-8 * self.J_to_eV	#eV m^-2 s^-1 K^-4
		x_He = self.Yp/(4*(1-self.Yp))
		gamma = (8*self.sigma_T*a_r*Tr**4)/(3*(1+x_He+x)*self.me) *x	# s^-1
		gamma /= self.c		# m^-1
		hubble /= self.Mpc_to_m
		hubble_class /= self.Mpc_to_m

		T21 = []
		redshift_distortion = []
		redshift = []
		wavenumber = []
		zz = np.linspace(400,0,1000)
		T_T = np.interp(zz[::-1], self.z_HYREC[::-1], self.T_T[::-1])[::-1]
		T_H = np.interp(zz[::-1], self.z_HYREC[::-1], self.T_H[::-1])[::-1]
		T_b = np.interp(zz[::-1], self.z_HYREC[::-1], self.T_b[::-1])[::-1]
		
		for i in [self.number_of_k2-1]:#range(number_of_k):#[20, 28, 184, 443, 551, 584,594]:#range(20,number_of_k):
			kk = self.klist2[::-1][i]
			print (i, kk)
			b = self.baryon[i*self.number_of_z2:(i+1)*self.number_of_z2]
			b = np.interp (new_z[::-1], self.zlist2, b[::-1])[::-1]
			b_dot = self.baryon_dot[i*self.number_of_z2:(i+1)*self.number_of_z2]*(-hubble_class)
			b_dot = np.interp (new_z[::-1], self.zlist2, b_dot[::-1])[::-1] 
			
			C_list = [0]
			dCdz_list = []
			C = 0
			for j in range(len(new_z)-1):
				dz = new_z[j]-new_z[j+1]
				dCdz = 1/((1+new_z[j])*hubble[j]*b[j]) * (-(1+new_z[j])*2/3*b_dot[j] + ((1+new_z[j])*b_dot[j]+Tr[j]/Tm[j]*gamma[j]*b[j])*C)
				C -= dCdz*dz
				dCdz_list.append (dCdz)
				C_list.append (C)
			C_list = np.array (C_list)
			C1 = np.interp(zz[::-1], new_z[::-1], C_list[::-1])[::-1]
		
			transfer_21 = (T_H + T_T*C1)
			distortion = T_b
			T21 += list(transfer_21)
			redshift_distortion += list(distortion)
			redshift += list(zz)
			wavenumber += list(np.ones(len(zz))*kk)
		
		self.T21 = np.array (T21)
		self.redshift_distortion = np.array (redshift_distortion)
		self.zlist = np.array (redshift)
		wavenumber = np.array (wavenumber)
		self.zz = zz	
		#data = np.column_stack((redshift, wavenumber, T21, redshift_distortion, hubble_list, baryon))
		#np.savetxt('transfer_21_nonu_Neff.txt', data, fmt = '%1.6e')
		#np.savetxt(outfile, data, fmt = '%1.6e')
		
			
	def cl21T (self, z_m, w):
		""" Calculate cross-correlation functions of ISW and 21 cm """

		z, sel, _ = sf.run_sel (w, z_m)
		
		chi_class_local = np.interp (z, self.zlist2, self.chi_class)
		hubble_local = np.interp (z, self.zlist2, self.hubble_class)
		dphidz = np.loadtxt(self.infile_new)[0:,5]
		T_dphidz = []
		T_baryon = []
		for i in range(self.number_of_k2):
			p = dphidz[self.number_of_z2*i:self.number_of_z2*(i+1)][::-1]
			T_dphidz.append (p)
			bb = self.baryon[self.number_of_z2*i:self.number_of_z2*(i+1)][::-1]
			T_baryon.append (bb)
		T_dphidz = interp2d (self.zlist2, self.klist2, T_dphidz[::-1], kind = 'quintic')
		T_baryon = interp2d (self.zlist2, self.klist2, T_baryon[::-1], kind = 'quintic')

		delta_21 = self.T21[::-1].copy ()
		zz = self.zz[::-1].copy ()
		cl_list = []
		for l in self.l_list:
			print (l)
			kk = (l+1/2)/chi_class_local
			P_phi_local = np.interp (kk, self.k_list, self.P_phi)
			
			transfer_21 = []
			transfer_dphidz = []
			for j in range(len(kk)):
				#T = delta_21 (z[j], kk[j])[0]
				T = np.interp (z[j], zz, delta_21)
				bb = T_baryon (z[j], kk[j])[0]
				transfer_21.append (T*bb)
				p = T_dphidz (z[j], kk[j])[0]
				transfer_dphidz.append (p)
			transfer_21 = np.array (transfer_21)
			transfer_dphidz = np.array (transfer_dphidz)
			"""	
			plt.figure(1)
			plt.plot (z,P_phi_local)
			plt.figure(2)
			plt.plot (z,transfer_21)
			plt.figure(3)
			plt.plot (z,transfer_dphidz)
			plt.figure(4)
			plt.plot (z,hubble_local)
			plt.figure(5)
			plt.plot (z,chi_class_local)
			plt.show()
			"""
			
			integrand = -2 * P_phi_local * sel * transfer_21 * transfer_dphidz * hubble_local / chi_class_local**2
			cl = simps (integrand, z)
			cl_list.append (cl)
	
		cl_list = np.array (cl_list)
		
		return cl_list



	def cl21 (self, z_m, w):
		""" Calculate 21 cm auto-correlation functions """
		
		z = np.linspace(10,400,1000)
		z1, sel1, _ = sf.run_sel (w[0], z_m[0])
		z2, sel2, _ = sf.run_sel (w[1], z_m[1])
		sel1 = np.interp (z, z1, sel1)
		sel2 = np.interp (z, z2, sel2)

		chi_class_local = np.interp (z, self.zlist2, self.chi_class)
		hubble_local = np.interp (z, self.zlist2, self.hubble_class)
		
		T_baryon = []
		for i in range(self.number_of_k2):
			bb = self.baryon[self.number_of_z2*i:self.number_of_z2*(i+1)][::-1]
			T_baryon.append (bb)
		T_baryon = interp2d (self.zlist2, self.klist2, T_baryon[::-1], kind = 'quintic')
		
		"""
		delta_21 =[]
		distortion = []
		for i in range(number_of_k):
			bb = self.baryon[number_of_z*i:number_of_z*(i+1)][::-1]

			d = self.T21[number_of_z*(number_of_k-1):number_of_z*number_of_k][::-1]
			d = d*bb
			delta_21.append (d)
			d = self.redshift_distortion[number_of_z*i:number_of_z*(i+1)][::-1]
			d = d*bb
			distortion.append (d)
		delta_21 = interp2d (zlist, klist, delta_21[::-1], kind = 'quintic')
		distortion = interp2d (zlist, klist, distortion[::-1], kind = 'quintic')
		"""
		delta_21 = self.T21[::-1].copy ()
		zz = self.zz[::-1].copy ()

		cl_list = []
		for l in self.l_list:
			print (l)
			kk = (l+1/2)/chi_class_local
			P_phi_local = np.interp (kk, self.k_list, self.P_phi)
			
			transfer_21 = []
			for j in range(len(kk)):
				T = np.interp (z[j], zz, delta_21)
				bb = T_baryon (z[j], kk[j])[0]
				transfer_21.append (T*bb)
			transfer_21 = np.array (transfer_21)
	
			integrand = 2*np.pi**2/l**3 * kk**3*P_phi_local/(2*np.pi**2) * sel1 * sel2 * transfer_21**2 * chi_class_local * hubble_local
			cl = simps (integrand, z)
			cl_list.append (cl)
	
		cl_list = np.array (cl_list)
		
		return cl_list


	def cl21_sharp (self, z_m, w):
		
		z, sel, _ = sf.run_sel (w, z_m)

		chi_class_local = np.interp (z, self.redshift2, self.chi_class)
		hubble_local = np.interp (z, self.redshift2, self.hubble_class)
		
		dphidz = np.loadtxt(self.infile_new)[0:,5]
		T_baryon = []
		for i in range(self.number_of_k2):
			bb = self.baryon[self.number_of_z2*i:self.number_of_z2*(i+1)][::-1]
			T_baryon.append (bb)
		T_baryon = interp2d (self.zlist2, self.klist2, T_baryon[::-1], kind = 'quintic')
		
		"""
		delta_21 =[]
		distortion = []
		for i in range(number_of_k):
			bb = self.baryon[number_of_z*i:number_of_z*(i+1)][::-1]
			d = self.T21[number_of_z*(number_of_k-1):number_of_z*number_of_k][::-1]
			d = d*bb
			delta_21.append (d)
			d = self.redshift_distortion[number_of_z*i:number_of_z*(i+1)][::-1]
			d = d*bb
			distortion.append (d)
		delta_21 = interp2d (zlist, klist, delta_21[::-1], kind = 'quintic')
		distortion = interp2d (zlist, klist, distortion[::-1], kind = 'quintic')
		"""
		
		delta_21 = self.T21[0:self.number_of_z][::-1]
		distortion = self.redshift_distortion[0:self.number_of_z][::-1]

		chi_class_z_m = np.interp (z_m, z, chi_class_local) 
		cl_list = []
		for l in l_list:
			print (l)
			jl = scipy.special.spherical_jn (int(l), self.k_list*chi_class_z_m)
			jl_2 = scipy.special.spherical_jn (int(l), self.k_list*chi_class_z_m, 2) 
			
			transfer_21 = []
			distor = []
			
			transfer_21_z_m = np.interp (z_m, self.zlist, delta_21)
			distor_z_m = np.interp (z_m, self.zlist, distortion)
			
			for j in range(len(self.k_list)):
				bb = T_baryon (z_m, k_list[j])[0]
				transfer_21.append (transfer_21_z_m*bb)
				distor.append (distor_z_m*bb)
				
			transfer_21 = np.array (transfer_21)
			distor = np.array (distor)
			
			integrand = 2/np.pi * self.P_phi * self.k_list**2 * (transfer_21**2*jl**2 + 2*transfer_21*distor*jl*jl_2 + distor**2*jl_2**2)
			cl = simps (integrand, self.k_list)
			cl_list.append (cl)
	
		cl_list = np.array (cl_list)
		print (cl_list)
		return cl_list



	def cl21_exact (self, z_m, w):
		
		z, sel, _ = sf.run_sel (w, z_m)

		chi_class_local = np.interp (z, self.redshift2, self.chi_class)
		hubble_local = np.interp (z, self.redshift2, self.hubble_class)
		
		T21 = self.T21[0:self.number_of_z][::-1]
		redshift_distortion = self.redshift_distortion[0:self.number_of_z][::-1]
		
		delta_21 =[]
		distortion = []
		for i in range(self.number_of_k2):
			bb = self.baryon[self.number_of_z2*i:self.number_of_z2*(i+1)][::-1]
			delta_21.append (T21*bb)
			distortion.append (redshift_distortion*bb)
		delta_21 = interp2d (self.zlist, self.klist, delta_21[::-1], kind = 'quintic')
		distortion = interp2d (self.zlist, self.klist, distortion[::-1], kind = 'quintic')
	
		"""
		delta_21 =[]
		distortion = []
		for i in range(number_of_k):
			bb = baryon[number_of_z*i:number_of_z*(i+1)][::-1]
				
			#d = T21[number_of_z*i:number_of_z*(i+1)][::-1]
			d = T21[number_of_z*(number_of_k-1):number_of_z*number_of_k][::-1]
			d = d*bb
			delta_21.append (d)
			d = redshift_distortion[number_of_z*i:number_of_z*(i+1)][::-1]
			d = d*bb
			distortion.append (d)
		delta_21 = interp2d (zlist, klist, delta_21[::-1], kind = 'quintic')
		distortion = interp2d (zlist, klist, distortion[::-1], kind = 'quintic')
		"""
		
		cl_list = []
		for l in l_list:
			print (l)
			alpha_list = []
			beta_list = []
			
			#P0_list = []
			#P0v_list = []
			#Pv_list = []
			for j in range(len(self.k_list)):
				jl = scipy.special.spherical_jn (int(l), self.k_list[j]*chi_class_local)
				jl_2 = scipy.special.spherical_jn (int(l), self.k_list[j]*chi_class_local, 2) 
				
				transfer_21 = delta_21 (z, self.k_list[j])[0]
				distor = distortion (z, self.k_list[j])[0]
				alpha = simps (jl*sel*transfer_21, z)
				beta = simps (jl_2*sel*distor, z)
				
				#transfer_21 = delta_21 (z_m, k_list[j])[0]
				#distor = distortion (z_m, k_list[j])[0]
				#alpha = simps (jl*sel, redshift)
				#beta = simps (jl_2*sel, redshift)
				alpha_list.append (alpha)
				beta_list.append (beta)
			
				#P0 = P_phi[j]*transfer_21**2
				#P0v = P_phi[j]*transfer_21*distor
				#Pv = P_phi[j]*distor**2
				#P0_list.append (P0)
				#P0v_list.append (P0v)
				#Pv_list.append (Pv)
	
			alpha_list = np.array (alpha_list)
			beta_list = np.array (beta_list)
			#P0_list = np.array (P0_list)
			#P0v_list = np.array (P0v_list)
			#Pv_list = np.array (Pv_list)
	
			integrand = 2/np.pi * self.P_phi * self.k_list**2 * (alpha_list**2 + 2*alpha_list*beta_list + beta_list**2)
			#integrand1 = 2/np.pi * k_list**2 * (P0_list*alpha_list**2 + P0v_list*2*alpha_list*beta_list + Pv_list*beta_list**2)
			#integrand1 = 2/np.pi * P_phi * k_list**2 * alpha_list**2
			#integrand2 = 2/np.pi * P_phi * k_list**2 * 2*alpha_list*beta_list 
			#integrand3 = 2/np.pi * P_phi * k_list**2 * beta_list**2
			
			cl = simps (integrand, self.k_list)
			#cl_ab = simps (integrand2, k_list)
			#cl_bb = simps (integrand3, k_list)
			
			cl_list.append (cl)
			#cl_list_ab.append (cl_ab)
			#cl_list_bb.append (cl_bb)
	
		cl_list = np.array (cl_list)
		#cl_list_ab = np.array (cl_list_ab)
		#cl_list_bb = np.array (cl_list_bb)
		print (cl_list)
		
		return cl_list

"""
#cl, cl2 = cl21 (z_m)
cl, cl2 = cl21T (z_m)
data = np.column_stack((l_list, cl2, cl))
#np.savetxt('cl_21cm_l10^7_6MHz_fix.txt', data, fmt = '%1.6e')
np.savetxt(outfile, data, fmt = '%1.6e')
#np.savetxt('result/antony/nanoom/21T_nu2_{}.txt'.format(z_m), data, fmt = '%1.6e')
#print ('21T_nu2_{}.txt'.format(z_m))
aa = 2.7255*10**6
plt.plot (l_list, np.sqrt(aa*cl*l_list*(l_list+1)/(2*np.pi)), label = 'limber2')
plt.plot (l_list, np.sqrt(aa*cl2*l_list*(l_list+1)/(2*np.pi)), label = 'limber')
#plt.plot (l_list, np.sqrt(aa*mono/l_list), label = 'mono')
#plt.plot (l_list, np.sqrt(aa*monovel/l_list), label = 'monovel')
#plt.plot (l_list, np.sqrt(aa*vel/l_list), label = 'vel')
#plt.plot (l_list, np.sqrt(aa*(mono+monovel+vel)/l_list), label = 'approx')
plt.legend ()
#plt.axis([2,10**7,10**-3,10])
plt.xscale ('log')
plt.yscale ('log')
"""
"""
aa = 2.7255*10**6
#cl = cl21_sharp ()
#cl, ab, bb = cl21_exact ()
cl = cl21_exact ()
#data = np.column_stack((l_list, cl,ab,bb))
data = np.column_stack((l_list, cl))
#np.savetxt('result/antony/antony_sharp_10000.txt', data, fmt = '%1.6e')
#print ('antony_sharp.txt, 10000')
np.savetxt('result/antony/antony_6MHz_5000.txt', data, fmt = '%1.6e')
print ('antony_exact_6MHz_5000.txt')
#print ('antony_exact_6MHz_100000, l=2, not saved')

#np.savetxt('cl_21cm_l10^7_exact_novel_01MHz.txt', data, fmt = '%1.6e')
#np.savetxt('cl_21cm_l10^7_sharp2.txt', data, fmt = '%1.6e')
aa = 2.7255*10**6
plt.plot (l_list, np.sqrt(aa*cl*l_list*(l_list+1)/(2*np.pi)), label = 'aa')
#plt.plot (l_list, np.sqrt(aa*ab*l_list*(l_list+1)/(2*np.pi)), label = 'ab')
#plt.plot (l_list, np.sqrt(aa*bb*l_list*(l_list+1)/(2*np.pi)), label = 'bb')
#plt.plot (l_list, np.sqrt(aa*cl*l_list*(l_list+1)/(2*np.pi)), label = 'sharp')
plt.legend ()
plt.axis([2,10**7,10**-3,10])
plt.xscale ('log')
plt.yscale ('log')

print (np.sqrt(aa*cl*l_list*(l_list+1)/(2*np.pi)))
plt.show()
"""

