import numpy as np
import matplotlib.pyplot as plt
#import data_no_nu as data
import data
import matplotlib.ticker

def test_run (infile_syn, infile_HYREC):
	infile = data.transfer_syn
	infile = infile_syn
	z = np.loadtxt(infile)[0:,0]
	print ('z data', len(z))
	zlist = sorted(set(z))
	zlist = np.array(zlist)
	rev_zlist = zlist[::-1]
	number_of_z = len(zlist)
	hubble_class = np.loadtxt(infile)[0:number_of_z,4][::-1]
	print ('hubble', len(hubble_class))

	c = 299792458
	Yp = 0.245
	Omega_b = 0.041283542878588854
	h = 0.73
	m_s = (299792458)**-1
	H_0 = 10**5 * m_s * h
	eV_to_K = 11604.5250061657
	eV_to_m_inv = 5076142.131979696
	Mpc_to_m = 3.0857*10**22
	rho_cr = 8.056*10**-11 * h**2 # eV^4
	mp = 938.2720813*10**6  #eV
	sigma_T = 6.6524587158 * 10 **-29


	infile = 'output_nanoom.dat'
	infile = 'output_nanoom_Neff.dat'
	infile = infile_HYREC
	z = np.loadtxt(infile)[0:,0]
	x = np.loadtxt(infile)[0:,1]
	Tm_Tr = np.loadtxt(infile)[0:,2]
	T_cmb = 2.7255
	Tr = T_cmb*(1+z)
	Tm = Tr * Tm_Tr
	
	hubble = np.interp (z[::-1], zlist, hubble_class)#/Mpc_to_m
	hubble = hubble[::-1]
	n_H = (1-Yp)*rho_cr*Omega_b/mp*(1+z)**3 #eV^3
	n_H *= eV_to_m_inv**3	#/m^3
	kappa = 3.1*10**-11 *Tm**0.357 * np.e**(-32/Tm) * 10**-6 # m^3/s
	kappa /= c	# m^2
	k_B = 8.613303*10**-5 # eV/K
	wavelength = 0.21 #m
	A10 = 2.85*10**-15 #/s
	A10 /= c	#/m
	E10 = 0.068	#K
	hc = 2*np.pi*197.3269788*10**6*10**-15 #eV m
	#B10 = (wavelength**3/(2*hc))*A10	#m/eV
	#B10 *= 1/(eV_to_m_inv) #m^2
	#B01 = 3*B10		#m^2
	B10 = A10*(1+1/(np.e**(E10/Tr)-1))	#/m
	B01 = 3*np.e**(-E10/Tr)*B10
	
	C10 = n_H*kappa		#/m
	C01 = 3*np.e**(-E10/Tm)*C10		#/m
	I_nu = 2*k_B*Tr/wavelength**2	#eV/m^2
	I_nu *= eV_to_m_inv		#/m^3
	gamma_list = []
	n1_n0_list = []
	Ts_list = []
	init = 7800

	hubble = hubble/Mpc_to_m
	Ts = Tr + (Tm-Tr)*C10/(C10+A10*Tm/E10)
	plt.figure(1)
	plt.plot(z+1, Tr, label = r'$T_{\mathrm{cmb}}$')
	plt.plot(z+1, Tm, label = r'$T_{\mathrm{gas}}$')
	#plt.loglog(z[init:]+1, Ts_list, label = 'spin')
	plt.plot(z+1, Ts, label = r'$T_s$')
	plt.xscale('log')
	#plt.yscale('log')
	plt.legend(bbox_to_anchor=(0.25,1))
	plt.grid()
	plt.xlabel (r'1+z')
	plt.ylabel (r'T (K)')
	#plt.savefig ('Ts.pdf', format='pdf')

	a = kappa*n_H
	b = A10*Tr/E10
	D = (1+z)**-1 * 3*E10/(32*np.pi)*1/4*Yp*wavelength**3*A10*n_H**2*kappa/(hubble*Tm*(kappa*n_H+A10*Tr/E10))
	C = (1+z)**-1 * 3*E10/(32*np.pi)*(1-x)*wavelength**3*A10
	#T_b = D* (Tm-Tr)
	#T_b = 3.1*10**-17 *C * E10*n_H**2 *(-Tr+Tm)/(hubble*(c*A10*np.e**(32/Tm)*Tr + 3.1*10**-17 *E10*n_H*Tm**0.357)*Tm**0.643)
	T_b = -C*E10*kappa*n_H**2*(Tr-Tm)/(hubble*(E10*kappa*n_H + A10*Tr)*Tm)
	#T_b = 0.00244785*C*E10*n_H**2*(-Tr+Tm)/(hubble*(A10*c*np.e**(32/Tm+32*Tm)*Tr+0.00244785*E10*n_H*Tm**0.357)*Tm**0.643)
	#T_H = D*(a+2*b)/(a+b)*(Tm-Tr)
	#T_H =  6.2*10**-17 * C * n_H**2*(c*A10*np.e**(32/Tm)*E10*Tr + 1.55*10**-17*E10**2*n_H*Tm**0.357)*(-Tr+Tm)/(hubble*(c*A10*np.e**(32/Tm)*Tr + 3.1*10**-17*E10*n_H*Tm**0.357)*Tm**0.643)
	T_H = -C*E10*kappa*n_H**2*(E10*kappa*n_H+2*A10*Tr)*(Tr-Tm)/(hubble*(E10*kappa*n_H+A10*Tr)**2*Tm)
	#T_T = D* (2*Tm_Tr)
	#T_T = 0.643*C*E10*kappa*n_H**2*(Tr+0.55521*Tm)/(hubble*(E10*kappa*n_H+A10*Tr)*Tm)
	T_T = (0.357*C*E10*kappa*n_H**2*(-89.6359*A10*Tr**2+89.6359*A10*Tr*Tm+2.80112*E10*kappa*n_H*Tr*Tm+1.80112*A10*Tr**2*Tm+A10*Tr*Tm**2))/(hubble*(E10*kappa*n_H+A10*Tr)**2*Tm**2)
	#T_HH = D* (Tm-Tr + a**2/(a+b)**2)
	T_HH = - A10**2*C*E10*kappa*n_H**2*Tr**2*(Tr-Tm)/(hubble*(E10*kappa*n_H+A10*Tr)**3*Tm)
	
	#T_HT = D* (4*Tm-(3*a+2*b)/(a+b)*Tr)
	#T_HT = C*E10*kappa*n_H**2*(E10**2*kappa**2*n_H**2*Tr + 3*A10*E10*kappa*n_H*Tr**2+1.286*A10**2*Tr**3+0.714*A10**2*Tr**2*Tm)/(hubble*(E10*kappa*n_H+A10*Tr)**3*Tm)
	T_HT = (0.714*C*E10*kappa*n_H**2*(-89.6359*A10**2*Tr**3+1.40056*E10**2*kappa**2*n_H**2*Tr*Tm+89.6359*A10**2*Tr**2*Tm+4.20168*A10*E10*kappa*n_H*Tr**2*Tm+1.80112*A10**2*Tr**3*Tm+A10**2*Tr**2*Tm**2))/(hubble*(E10*kappa*n_H+A10*Tr)**3*Tm**2)
	#T_TT = D* (2*Tm-Tr)
	#T_TT = - (C*kappa*n_H**2 *(E10**3*kappa**2*n_H**2*Tr+1.40078*A10*E10**2*kappa*n_H*Tr**2+0.528225*A10**2*E10*Tr**3+0.242225*A10*E10**2*kappa*n_H*Tr*Tm+0.114776*A10**2*E10*Tr**2*Tm))/(hubble*(E10*kappa*n_H+A10*Tr)**3*Tm) 
	#T_TT = -(C*E10*kappa*n_H**2*(A10**2*Tr**2*(Tm*(-512-11.424*Tm+0.114776*Tm**2)+Tr*(512-20.576*Tm+0.528225*Tm**2)) + E10**2 *kappa**2*n_H**2*(Tr*(2.27374*10**-13-1.77636*10**-15*Tm+1.11022*10**-16*Tm**2)+Tr*(-2.27374*10**-13+1.77636*10**-15*Tm+Tm**2))+A10*E10*kappa*n_H*Tr*(Tm*(512+11.424*Tm+0.24225*Tm**2)+Tr*(-512-43.424*Tm+1.40078*Tm**2))))/(hubble*(E10*kappa*n_H+A10*Tr)**3*Tm**3)
	
	T_TT = (C*E10*kappa*n_H**2*(512*A10*E10*kappa*n_H*Tr**2-512*A10**2*Tr**3-512*A10*E10*kappa*n_H*Tr*Tm+512*A10**2*Tr**2*Tm+75.424*A10*E10*kappa*n_H*Tr**2*Tm+52.576*A10**2*Tr**3*Tm-43.424*A10*E10*kappa*n_H*Tr*Tm**2-E10**2*kappa**2*n_H**2*Tr*Tm**2-20.576*A10**2*Tr**2*Tm**2-1.40078*A10*E10*kappa*n_H*Tr**2*Tm**2-0.528225*A10**2*Tr**3*Tm**2-1.11022*10**-16*E10**2*kappa**2*n_H**2*Tm**3-0.242225*A10*E10*kappa*n_H*Tr*Tm**3-0.114776*A10**2*Tr**2*Tm**3))/(hubble*(E10*kappa*n_H+A10*Tr)**3*Tm**3)
	plt.figure(2)
	fig = plt.subplot(1,1,1)
	plt.plot(z+1, T_b*10**3, label = r'$\bar{\mathcal{T}}_b$')
	plt.plot(z+1, T_H*10**3, label = r'$\mathcal{T}_H$')
	plt.plot(z+1, T_T*10**3, label = r'$\mathcal{T}_T$')
	plt.plot(z+1, T_HH*10**3, label = r'$\mathcal{T}_{HH}$')
	plt.plot(z+1, T_HT*10**3, label = r'$\mathcal{T}_{HT}$')
	plt.plot(z+1, T_TT*10**3, label = r'$\mathcal{T}_{TT}$')
	plt.xscale('log', basex=2)
	plt.grid()
	plt.xlabel (r'1+z')
	plt.ylabel (r'Coefficient of $T_b$ (mK)')
	plt.legend(bbox_to_anchor=(0.2,1),prop={'size':12})
	plt.axis([20,200,-130,140])
	fig.set_xticks([20,30,50,70,100,150,200])
	fig.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	#plt.savefig ('Tb_coeffi.pdf', format='pdf')

	return z, T_T, T_H, T_b
