import numpy as np
from scipy.interpolate import interp2d
from scipy.integrate import simps
import selection_function21cm as sf
import matplotlib.pyplot as plt
import scipy.special
import sys

class cl_21 (object):
	def __init__ (self, infile_syn, infile_new, infile_21):

		self.infile_syn = infile_syn
		self.infile_new = infile_new
		self.infile_21 = infile_21

		z = np.loadtxt(self.infile_syn)[0:,0]
		print ('z data', len(z))
		k = np.loadtxt(self.infile_syn)[0:,1]
		print ('k data', len(k))
		zlist2 = np.array(sorted(set(z)))
		klist2 = np.array(sorted(set(k)))
		number_of_z2 = len(zlist2)
		number_of_k2 = len(klist2)
		hubble_class = np.loadtxt(self.infile_syn)[0:number_of_z2,4][::-1]
		print ('hubble data', len(hubble_class))
		self.baryon = -np.loadtxt(self.infile_syn)[0:,2]

		As = (np.e**(3.062))*10**-10
		n_s = 0.968
		k_pivot = 0.05
		self.k_list = np.logspace (-4,5, 5000)
		self.P_phi = As * (self.k_list/k_pivot)**(n_s-1) * 2*np.pi**2 / self.k_list**3

		a_0 = 10**-2/4
		n = 10000
		scale_factor = np.logspace (np.log10(a_0), 0, n)
		scale_factor_reverse = scale_factor[::-1]
		redshift2 = 1/scale_factor_reverse - 1
		hubble_class = np.interp (redshift2, zlist2, hubble_class)

		chi_class = []
		for i in range(len(hubble_class)):
			chi_class.append (simps (1/hubble_class[:i+1], redshift2[:i+1]))
		chi_class = np.array (chi_class)

		#21 data
		z = np.loadtxt(self.infile_21)[0:,0]
		print ('z data', len(z))
		k = np.loadtxt(self.infile_21)[0:,1]
		print ('k data', len(k))
		zlist = sorted(set(z))
		klist = sorted(set(k))
		number_of_z = len(zlist)
		number_of_k = len(klist)
		self.T21 = np.loadtxt(self.infile_21)[0:,2]
		print ('T21 data', len(self.T21))
		self.redshift_distortion = np.loadtxt(self.infile_21)[0:,3]
		print ('redshift distortion', len(self.redshift_distortion))
		#self.baryon = np.loadtxt(self.infile_21)[0:,5]
		

		l_list = np.logspace(np.log10(2), np.log10(5000), 1000)
		l_list = np.logspace(np.log10(2), np.log10(100), 100)
		for i in range(len(l_list)):
			l_list[i] = int(l_list[i])
		l_list = sorted(set(l_list))
		#l_list[-1] += 1
		l_list = np.array (l_list)

		self.hubble_class = hubble_class
		self.chi_class = chi_class
		self.l_list = l_list
		self.zlist = zlist
		self.klist = klist
		self.zlist2 = zlist2
		self.klist2 = klist2
		self.number_of_z = number_of_z
		self.number_of_k = number_of_k
		self.number_of_z2 = number_of_z2
		self.number_of_k2 = number_of_k2
		self.redshift2 = redshift2

	def cl21T (self, z_m, w):
		z, sel, _ = sf.run_sel (w, z_m)

		chi_class_local = np.interp (z, self.redshift2, self.chi_class)
		hubble_local = np.interp (z, self.redshift2, self.hubble_class)
		
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
		"""
		delta_21 =[]
		distortion = []
		for i in range(self.number_of_k):
			#bb = self.baryon[self.number_of_z*i:self.number_of_z*(i+1)][::-1]
	
			d = self.T21[self.number_of_z*(self.number_of_k-1):self.number_of_z*self.number_of_k][::-1]
			#d = d*bb
			delta_21.append (d)
			d = self.redshift_distortion[self.number_of_z*i:self.number_of_z*(i+1)][::-1]
			#d = d*bb
			distortion.append (d)
		delta_21 = interp2d (self.zlist, self.klist, delta_21[::-1], kind = 'quintic')
		distortion = interp2d (self.zlist, self.klist, distortion[::-1], kind = 'quintic')
		"""

		delta_21 = self.T21[0:self.number_of_z][::-1]

		cl_list = []
		for l in self.l_list:
			print (l)
			kk = (l+1/2)/chi_class_local
			P_phi2 = np.interp (kk, self.k_list, self.P_phi)
			
			transfer_21 = []
			transfer_dphidz = []
			for j in range(len(kk)):
				#T = delta_21 (z[j], kk[j])[0]
				T = np.interp (z[j], self.zlist, delta_21)
				bb = T_baryon (z[j], kk[j])[0]
				transfer_21.append (T*bb)
				p = T_dphidz (z[j], kk[j])[0]
				transfer_dphidz.append (p)
			transfer_21 = np.array (transfer_21)
			transfer_dphidz = np.array (transfer_dphidz)
	
			integrand2 = -2 * P_phi2 * sel * transfer_21 * transfer_dphidz * hubble_local / chi_class_local**2
			cl = simps (integrand2, z)
			cl_list.append (cl)
	
		cl_list = np.array (cl_list)
		
		return cl_list



def cl21 (z_m):
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
	
	cl_list = []
	cl_list2 = []
	for l in l_list:
		print (l)
		kk = (l+1/2)/chi_class
		P_phi2 = np.interp (kk, k_list, P_phi)
		
		transfer_21 = []
		for j in range(len(kk)):
			T = delta_21 (redshift[j], kk[j])[0]
			transfer_21.append (T)
		transfer_21 = np.array (transfer_21)
		
		chi_class_z_m = np.interp (z_m, redshift, chi_class)
		transfer_21_z_m = np.interp (z_m, redshift, transfer_21)
		kk_z_m = (l+1/2)/chi_class_z_m
		integrand = 2*np.pi**2/l**3 * kk_z_m**3*P_phi2/(2*np.pi**2) * sel**2 * transfer_21_z_m**2 * chi_class_z_m * hubble
		cl = simps (integrand, redshift)
		cl_list.append (cl)
		
		integrand2 = 2*np.pi**2/l**3 * kk**3*P_phi2/(2*np.pi**2) * sel**2 * transfer_21**2 * chi_class * hubble
		cl2 = simps (integrand2, redshift)
		cl_list2.append (cl2)

		
	cl_list = np.array (cl_list)
	cl_list2 = np.array (cl_list2)
	
	return cl_list, cl_list2

def cl21_sharp ():
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

	z_m = 50
	chi_class2 = np.interp (z_m, redshift, chi_class) 
	cl_list = []
	for l in l_list:
		print (l)
		jl = scipy.special.spherical_jn (int(l), k_list*chi_class2)
		jl_2 = scipy.special.spherical_jn (int(l), k_list*chi_class2, 2) 
		
		transfer_21 = []
		distor = []
		for j in range(len(k_list)):
			T = delta_21 (z_m, k_list[j])[0]
			transfer_21.append (T)
			T = distortion (z_m, k_list[j])[0]
			distor.append (T)
		transfer_21 = np.array (transfer_21)
		distor = np.array (distor)
		
		integrand2 = 2/np.pi * P_phi * k_list**2 * (transfer_21**2*jl**2 + 2*transfer_21*distor*jl*jl_2 + distor**2*jl_2**2)
		cl = simps (integrand2, k_list)
		cl_list.append (cl)

	cl_list = np.array (cl_list)
	print (cl_list)
	return cl_list

def cl21_exact ():
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

	z_m = 50
	cl_list = []
	cl_list_ab = []
	cl_list_bb = []
	for l in l_list:
		print (l)
		alpha_list = []
		beta_list = []
		
		P0_list = []
		P0v_list = []
		Pv_list = []
		for j in range(len(k_list)):
			jl = scipy.special.spherical_jn (int(l), k_list[j]*chi_class)
			jl_2 = scipy.special.spherical_jn (int(l), k_list[j]*chi_class, 2) 
			
			#transfer_21 = delta_21 (redshift, k_list[j])[0]
			#distor = distortion (redshift, k_list[j])[0]
			#alpha = simps (jl*sel*transfer_21, redshift)
			#beta = simps (jl_2*sel*distor, redshift)
			
			transfer_21 = delta_21 (z_m, k_list[j])[0]
			distor = distortion (z_m, k_list[j])[0]
			alpha = simps (jl*sel, redshift)
			beta = simps (jl_2*sel, redshift)
			alpha_list.append (alpha)
			beta_list.append (beta)
		
			P0 = P_phi[j]*transfer_21**2
			P0v = P_phi[j]*transfer_21*distor
			Pv = P_phi[j]*distor**2
			P0_list.append (P0)
			P0v_list.append (P0v)
			Pv_list.append (Pv)

		alpha_list = np.array (alpha_list)
		beta_list = np.array (beta_list)
		P0_list = np.array (P0_list)
		P0v_list = np.array (P0v_list)
		Pv_list = np.array (Pv_list)

		#integrand1 = 2/np.pi * P_phi * k_list**2 * (alpha_list**2 + 2*alpha_list*beta_list + beta_list**2)
		integrand1 = 2/np.pi * k_list**2 * (P0_list*alpha_list**2 + P0v_list*2*alpha_list*beta_list + Pv_list*beta_list**2)
		#integrand1 = 2/np.pi * P_phi * k_list**2 * alpha_list**2
		#integrand2 = 2/np.pi * P_phi * k_list**2 * 2*alpha_list*beta_list 
		#integrand3 = 2/np.pi * P_phi * k_list**2 * beta_list**2
		
		cl = simps (integrand1, k_list)
		#cl_ab = simps (integrand2, k_list)
		#cl_bb = simps (integrand3, k_list)
		
		cl_list.append (cl)
		#cl_list_ab.append (cl_ab)
		#cl_list_bb.append (cl_bb)

	cl_list = np.array (cl_list)
	#cl_list_ab = np.array (cl_list_ab)
	#cl_list_bb = np.array (cl_list_bb)
	print (cl_list)
	#return cl_list, cl_list_ab, cl_list_bb
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

