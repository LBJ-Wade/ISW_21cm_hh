import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import test
#import data_no_nu as data
import data
import sys



infile_syn = sys.argv[1]
infile_HYREC = sys.argv[2]
outfile = sys.argv[3]

test_z, test_T_T, test_T_H, test_T_b = test.test_run (infile_syn, infile_HYREC)

#infile = data.transfer_syn
infile = infile_syn
z = np.loadtxt(infile)[0:,0]
print ('z data', len(z))
zlist = sorted(set(z))
zlist = np.array(zlist)
rev_zlist = zlist[::-1]
k = np.loadtxt(infile)[0:,1]
print ('k data', len(k))
klist = sorted(set(k))
rev_klist = klist[::-1]
number_of_z = len(zlist)
number_of_k = len(klist)
hubble_class = np.loadtxt(infile)[0:number_of_z,4]
print ('hubble', len(hubble_class))
delta_b = -np.loadtxt(infile)[0:,2]
delta_b_dot = -np.loadtxt(infile)[0:,3]
plt.figure(3)
plt.loglog(rev_zlist, hubble_class)

infile = 'output_nanoom.dat'
infile = 'output_nanoom_Neff.dat'
infile = infile_HYREC
z = np.loadtxt(infile)[0:,0]
x = np.loadtxt(infile)[0:,1]
Tm_Tr = np.loadtxt(infile)[0:,2]
T_cmb = 2.7255
Tr = T_cmb*(1+z)
Tm = Tr * Tm_Tr

#zz = np.logspace(3,0,10000)
#new_z = zz-1
new_z = np.linspace(1000,0,100000)
x = np.interp (new_z[::-1], z[::-1], x[::-1])[::-1]
hubble = np.interp (new_z[::-1], rev_zlist[::-1], hubble_class[::-1])[::-1]
Tr = np.interp (new_z[::-1], z[::-1], Tr[::-1])[::-1]
Tm = np.interp (new_z[::-1], z[::-1], Tm[::-1])[::-1]
Mpc_to_m = 3.0857*10**22
c = 299792458
sigma_T = 6.6524587158 * 10 **-29	#m^2
J_to_eV = 6.2415093433*10**18
#a_r = 7.5657*10**-16 *J_to_eV		#eV m^-3 K^-4
a_r = 4*5.670373 * 10**-8 * J_to_eV	#eV m^-2 s^-1 K^-4
Yp = 0.245
x_He = Yp/(4*(1-Yp))
me = 0.5109989461*10**6	#eV
gamma = (8*sigma_T*a_r*Tr**4)/(3*(1+x_He+x)*me) *x	# s^-1
gamma /= c		# m^-1
#gamma *= Mpc_to_m	#Mpc^-1
hubble /= Mpc_to_m
hubble_class /= Mpc_to_m
plt.figure(33)
fig = plt.subplot(1,1,1)
T21 = []
redshift_distortion = []
redshift = []
wavenumber = []
hubble_list = []
baryon = []
zz = np.linspace(400,0,1000)
hubble2 = np.interp (zz[::-1], new_z[::-1], hubble[::-1])[::-1] * Mpc_to_m
T_T = np.interp(zz[::-1], test_z[::-1], test_T_T[::-1])[::-1]
T_H = np.interp(zz[::-1], test_z[::-1], test_T_H[::-1])[::-1]
T_b = np.interp(zz[::-1], test_z[::-1], test_T_b[::-1])[::-1]
#plt.plot(test.z+1, test.T_b*1000, label = r'$\bar{T}_{\mathrm{b}}$')
	
kkkk = np.arange(59)*10
for i in range(number_of_k):#[20, 28, 184, 443, 551, 584,594]:#range(20,number_of_k):
	kk = rev_klist[i]
	print (i, kk)
	b = delta_b[i*number_of_z:(i+1)*number_of_z]
	b = np.interp (new_z[::-1], rev_zlist[::-1], b[::-1])[::-1]
	b_dot = delta_b_dot[i*number_of_z:(i+1)*number_of_z]*(-hubble_class)
	b_dot = np.interp (new_z[::-1], rev_zlist[::-1], b_dot[::-1])[::-1] #*(1+new_z)
	
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
	#plt.plot(new_z+1, C_list, label = 'k = {0:.4f}/Mpc'.format(kk))
	#plt.plot(new_z+1, C_list, label = 'k = {}/Mpc'.format(kk))
	
	
	c1 = np.interp(test_z[::-1],new_z[::-1],C_list[::-1])[::-1]
	C1 = np.interp(zz[::-1], new_z[::-1], C_list[::-1])[::-1]
	b2 = np.interp(zz[::-1], new_z[::-1], b[::-1])[::-1]
		

	transfer_21 = (T_H + T_T*C1)
	distortion = T_b
	baryon += list(b2)	
	T21 += list(transfer_21)
	redshift_distortion += list(distortion)
	redshift += list(zz)
	wavenumber += list(np.ones(len(zz))*kk)
	hubble_list += list(hubble2)

T21 = np.array (T21)
redshift_distortion = np.array (redshift_distortion)
redshift = np.array (redshift)
wavenumber = np.array (wavenumber)
hubble_list = np.array (hubble_list)

data = np.column_stack((redshift, wavenumber, T21, redshift_distortion, hubble_list, baryon))
#np.savetxt('transfer_21_nonu_Neff.txt', data, fmt = '%1.6e')
np.savetxt(outfile, data, fmt = '%1.6e')
#plt.legend(bbox_to_anchor=(0.8,1),prop={'size':12})
#plt.legend()
plt.xscale('log',basex=2)
#fig.set_xticks([20,50,100,200])
fig.set_xticks([5,10,50,100,500])
fig.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#plt.axis([20,200,-42,10])
plt.axis([2,1000,-0.2,0.7])
plt.xlabel (r'$1+z$')
#plt.ylabel (r'$T_b$, $\alpha$')
plt.ylabel (r'$C_1(z)$')
plt.grid()
#plt.savefig ('c_z_kmin1.pdf', format='pdf')
#plt.savefig ('alpha_kmin.pdf', format='pdf')

#plt.show()

