import numpy as np
import matplotlib.pyplot as plt
import time
import nanoom as nn
import matplotlib.ticker
from scipy.integrate import odeint
def maketab (xmin, xmax, Nx):
	h = (xmax - xmin)/(Nx-1.0)
	xtab = np.ones(Nx)*xmin + np.arange(Nx)*h
	return xtab	

def interpolate_rates (TR, TM_TR, logAlpha_tab, logR2p2s_tab):
	global TM_TR_MIN, DTM_TR, TR_MIN, DlogTR, NTM, NTR
	coeff1 = np.zeros(4)
	coeff2 = np.zeros(4)
	temp = np.zeros(4)
	logTR = np.log(TR)
	Alpha = np.zeros(2)
	Beta = np.zeros(2)
	R2p2s = None
	iTM = np.floor ( (TM_TR - TM_TR_MIN)/DTM_TR )
	if iTM < 1:
		iTM = 1
	if iTM > NTM-3:
		iTM = NTM-3
	frac1 = (TM_TR - TM_TR_MIN)/DTM_TR - iTM
	coeff1[0] = frac1*(frac1-1.)*(2.-frac1)/6.
	coeff1[1] = (1.+frac1)*(1.-frac1)*(2.-frac1)/2.
	coeff1[2] = (1.+frac1)*frac1*(2.-frac1)/2.
	coeff1[3] = (1.+frac1)*frac1*(frac1-1.)/6.

	iTR = np.floor ( (logTR - np.log(TR_MIN))/DlogTR )
	if iTR < 1:
		iTR = 1
	if iTR > NTR-3:
		iTR = NTR-3
	frac2 = (logTR - np.log(TR_MIN))/DlogTR - iTR
	coeff2[0] = frac2*(frac2-1.)*(2.-frac2)/6.
	coeff2[1] = (1.+frac2)*(1.-frac2)*(2.-frac2)/2.
	coeff2[2] = (1.+frac2)*frac2*(2.-frac2)/2.
	coeff2[3] = (1.+frac2)*frac2*(frac2-1.)/6.

	for l in range(2):
		for k in range(4):
			temp[k] = logAlpha_tab[l][iTM-1+k][iTR-1]*coeff2[0]\
					+ logAlpha_tab[l][iTM-1+k][iTR]*coeff2[1]\
					+ logAlpha_tab[l][iTM-1+k][iTR+1]*coeff2[2]\
					+ logAlpha_tab[l][iTM-1+k][iTR+2]*coeff2[3]
		Alpha[l] = np.exp(temp[0]*coeff1[0]+temp[1]*coeff1[1]\
					+temp[2]*coeff1[2]+temp[3]*coeff1[3])#/10**6
		Beta[l] = np.exp(logAlpha_tab[l][NTM-1][iTR-1]*coeff2[0]\
				+ logAlpha_tab[l][NTM-1][iTR]*coeff2[1]\
				+ logAlpha_tab[l][NTM-1][iTR+1]*coeff2[2]\
				+ logAlpha_tab[l][NTM-1][iTR+2]*coeff2[3])#/10**6

	factor = 3.016103031869581e21 *  TR*np.sqrt(TR) * np.exp(-3.399571517984581/TR)
	#factor = 3.0917400259648894e+27*  TR*np.sqrt(TR) * np.exp(-3.399571517984581/TR)
	Beta[0] *= factor
	Beta[1] *= factor/3.

	R2p2s = np.exp(logR2p2s_tab[iTR-1]*coeff2[0]\
			+ logR2p2s_tab[iTR]*coeff2[1]\
			+ logR2p2s_tab[iTR+1]*coeff2[2]\
			+ logR2p2s_tab[iTR+2]*coeff2[3])
	return Alpha, Beta, R2p2s	

ALPHA_FILE = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/class_syn/hyrec/Alpha_inf.dat"
RR_FILE = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/class_syn/hyrec/R_inf.dat"

NTR = 100
NTM = 40
TR_MIN = 0.004
TR_MAX = 0.4
TM_TR_MIN = 0.1
TM_TR_MAX = 1.0

logTR_tab = maketab (np.log(TR_MIN), np.log(TR_MAX), NTR)
TM_TR_tab = maketab (TM_TR_MIN, TM_TR_MAX, NTM)
DlogTR = logTR_tab[1] - logTR_tab[0]
DTM_TR = TM_TR_tab[1] - TM_TR_tab[0]


logR2p2s_tab = np.log(np.loadtxt (RR_FILE)[0:])
logAlpha_tab = np.zeros([2,NTM,NTR])
for i in range(NTR):
	Alpha1 = np.loadtxt(ALPHA_FILE)[i*NTM:(i+1)*NTM,0]
	Alpha2 = np.loadtxt(ALPHA_FILE)[i*NTM:(i+1)*NTM,1]
	logAlpha_tab[0,:,i] = np.log(Alpha1)
	logAlpha_tab[1,:,i] = np.log(Alpha2)

#Alpha, Beta, R2p2s = interpolate_rates (.004,0.1, logAlpha_tab,logR2p2s_tab)
#print (Alpha[0], Alpha[1])

infile_syn = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/data/delta_syn_0.dat"
infile_HYREC = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/data/output_prac_0.dat"


z = np.loadtxt(infile_syn)[0:,0]
print ('z data', len(z))
k = np.loadtxt(infile_syn)[0:,1]
print ('k data', len(k))
zlist2 = np.array(sorted(set(z)))
klist2 = np.array(sorted(set(k)))
number_of_z2 = len(zlist2)
number_of_k2 = len(klist2)
hubble_class = np.loadtxt(infile_syn)[0:number_of_z2,4][::-1]
print ('hubble data', len(hubble_class))
baryon = -np.loadtxt(infile_syn)[0:,2]
baryon_dot = -np.loadtxt(infile_syn)[0:,3]

a_0 = 10.**-2./4.
n = 10000
scale_factor = np.logspace (np.log10(a_0), 0, n)
scale_factor_reverse = scale_factor[::-1]
redshift2 = 1./scale_factor_reverse - 1.
hubble_class2 = np.interp (redshift2, zlist2, hubble_class)


z_HYREC = np.loadtxt(infile_HYREC)[0:,0]
x_HYREC = np.loadtxt(infile_HYREC)[0:,1]
Tm_Tr = np.loadtxt(infile_HYREC)[0:,2]
T_cmb = 2.7255
Tr = T_cmb*(1.+z_HYREC)
Tm = Tr * Tm_Tr
c = 299792458.
Mpc_to_m = 3.0857*10.**22.
h = 6.711034e-01
Omega_b = 0.0222/h**2
rho_cr = 8.056*10.**-11. * h**2.
mp = 938.2720813*10.**6.
me = 0.5109989461*10.**6.
sigma_T = 6.6524587158 * 10. **-29.
new_z = np.linspace(1000,0,20000)
k_to_eV = 1/11604.5250061657
Yp = 0.2467
x_He = Yp/(4*(1-Yp))
J_to_eV = 6.2415093433*10.**18.
eV_to_m_inv = 5076142.131979696
a_r = 4*5.670373 * 10**-8 * J_to_eV    #eV m^-2 s^-1 K^-4


As = (np.e**3.0892226166413299)*10**(-10)
n_s = 0.9655
k_pivot = 0.05
k_list = np.logspace (-4,5, 5000)
P_phi = As * (k_list/k_pivot)**(n_s-1.) * 2.*np.pi**2. / k_list**3

Hubble = np.interp (new_z[::-1], zlist2, hubble_class)[::-1]
Hubble = Hubble *c / Mpc_to_m
#Hubble = Hubble / Mpc_to_m
nH = (1-Yp)*rho_cr*Omega_b/mp*(1+new_z)**3		# eV^3
nH *= eV_to_m_inv**3		# m^-3
nH *= 1e-6
L2s1s = 8.2206

ind = 700
ind = 754
#ind = 340
ind = 600
ind = 100
kk = klist2[::-1][ind]
P_phi_kk = np.interp (kk, k_list, P_phi)
print kk, P_phi_kk
delta_b = baryon[ind*number_of_z2:(ind+1)*number_of_z2]
delta_b = np.interp (new_z[::-1], zlist2, delta_b[::-1])[::-1]

ddelta_b_dz = baryon_dot[ind*number_of_z2:(ind+1)*number_of_z2]
ddelta_b_dz = np.interp (new_z[::-1], zlist2, ddelta_b_dz[::-1])[::-1]
ddelta_b_deta = ddelta_b_dz * (-Hubble) 
theta_b = -ddelta_b_deta*(1+new_z) 

xe = np.interp (new_z[::-1], z_HYREC[::-1], x_HYREC[::-1])[::-1]
#RLya = 4.662899067555897e15 *Hubble /nH/(1.-xe)
#RLya = 4.662899067555897e15 *(Hubble + theta_b/3.) /nH/(1.-xe)
RLya = 4.662899067555897e15 *(Hubble + theta_b/3.*np.sqrt(P_phi_kk)) /nH/(1.-xe)
wavelength = 1216*10**-10


TR = np.interp (new_z[::-1], z_HYREC[::-1], Tr[::-1])[::-1]
TM = np.interp (new_z[::-1], z_HYREC[::-1], Tm[::-1])[::-1]

gamma = (8*sigma_T*a_r*TR**4)/(3*(1+x_He+xe)*me) *xe	# s^-1
#gamma /= c

TR *= k_to_eV
TM *= k_to_eV
TM_TR = TM/TR

Beta1 = []
Beta2 = []
Alpha1 = []
Alpha2 = []
for i in range(len(new_z)):
	Alpha, Beta, R2p2s = interpolate_rates (TR[i], TM_TR[i], logAlpha_tab,logR2p2s_tab)
	Beta1.append (Beta[0])
	Beta2.append (Beta[1])
	Alpha1.append (Alpha[0])
	Alpha2.append (Alpha[1])
AB = np.array(Alpha1) + np.array(Alpha2)
Beta = np.array(Beta1) + np.array(Beta2)
C = (3*RLya + L2s1s)/(3*RLya + L2s1s + 4*Beta)
dlogC_dlogRLya = -nn.derivative (np.log(RLya), np.log(C))

xe_dot = -C*AB*nH*xe**2

delta_xe = [0]
delta_TM = [0]
delta_xe_ini = 0
delta_TM_ini = 0

def model (y, t, y2):
	i = t
	return -1./(Hubble[i]*(1.+new_z[i])) * ( gamma[i]*( (TR[i]-TM[i])/TM[i]*y2*0 - TR[i]/TM[i]*y) + 2./3.*(1.+new_z[i])*ddelta_b_deta[i])

for i in range(len(new_z)-1):
	dz = new_z[i]-new_z[i+1]
	Alpha, Beta, R2p2s = interpolate_rates (TR[i], TM_TR[i], logAlpha_tab,logR2p2s_tab)
	
	
	dlogAB_dlogTM = []
	#for j in range(len(new_z)):
	for j in [i-1, i, i+1]:
		Alpha2, _, _ = interpolate_rates (TR[i], TM_TR[j], logAlpha_tab,logR2p2s_tab)
		#Alpha2, _, _ = interpolate_rates (TR[j], TM_TR[j], logAlpha_tab,logR2p2s_tab)
		dlogAB_dlogTM.append (np.log(Alpha2[0]+Alpha2[1]))
	dlogAB_dlogTM = np.array (dlogAB_dlogTM)
	#dlogAB_dlogTM = nn.derivative (TM[::-1], dlogAB_dlogTM[::-1])[::-1]
	dlogAB_dlogTM = nn.derivative ([np.log(TM[i-1]), np.log(TM[i]), np.log(TM[i+1])][::-1], dlogAB_dlogTM[::-1])[::-1]

	#ddelta_TM_dz = -1./(Hubble[i]*(1.+new_z[i])) * ( gamma[i]*( (TR[i]-TM[i])/TM[i]*delta_xe_ini*0 - TR[i]/TM[i]*delta_TM_ini) + 2./3.*(1.+new_z[i])*ddelta_b_deta[i])
	ddelta_TM_dz = odeint (model,delta_TM_ini,i, args=(delta_xe_ini,))[0] 
	
	delta_TM_ini -= ddelta_TM_dz*dz
	delta_TM.append (delta_TM_ini)

	ddelta_xe_dz = -1./(Hubble[i]*(1.+new_z[i]))*xe_dot[i]/xe[i] * (delta_xe_ini + delta_b[i] + dlogAB_dlogTM[1] * delta_TM_ini + dlogC_dlogRLya[i]*(theta_b[i]/(3.*Hubble[i])-delta_b[i]))
	print i, ddelta_xe_dz, ddelta_TM_dz
	delta_xe_ini -= ddelta_xe_dz*dz
	delta_xe.append (delta_xe_ini)
plt.figure(1)
fig = plt.subplot(1,1,1)
plt.xscale('log', basex=2)
#plt.yscale('log')
plt.plot(1+new_z, delta_TM*(1+new_z), label = r'Without $\delta_{x_e}$')


delta_xe = [0]
delta_TM = [0]
delta_xe_ini = 0
delta_TM_ini = 0
for i in range(len(new_z)-1):
	dz = new_z[i]-new_z[i+1]
	Alpha, Beta, R2p2s = interpolate_rates (TR[i], TM_TR[i], logAlpha_tab,logR2p2s_tab)
	
	
	dlogAB_dlogTM = []
	#for j in range(len(new_z)):
	for j in [i-1, i, i+1]:
		Alpha2, _, _ = interpolate_rates (TR[i], TM_TR[j], logAlpha_tab,logR2p2s_tab)
		#Alpha2, _, _ = interpolate_rates (TR[j], TM_TR[j], logAlpha_tab,logR2p2s_tab)
		dlogAB_dlogTM.append (np.log(Alpha2[0]+Alpha2[1]))
	dlogAB_dlogTM = np.array (dlogAB_dlogTM)
	#dlogAB_dlogTM = nn.derivative (TM[::-1], dlogAB_dlogTM[::-1])[::-1]
	dlogAB_dlogTM = nn.derivative ([np.log(TM[i-1]), np.log(TM[i]), np.log(TM[i+1])][::-1], dlogAB_dlogTM[::-1])[::-1]

	ddelta_TM_dz = -1./(Hubble[i]*(1.+new_z[i])) * ( gamma[i]*( (TR[i]-TM[i])/TM[i]*delta_xe_ini - TR[i]/TM[i]*delta_TM_ini) + 2./3.*(1.+new_z[i])*ddelta_b_deta[i])
	delta_TM_ini -= ddelta_TM_dz*dz
	delta_TM.append (delta_TM_ini)

	ddelta_xe_dz = -1./(Hubble[i]*(1.+new_z[i]))*xe_dot[i]/xe[i] * (delta_xe_ini + delta_b[i] + dlogAB_dlogTM[1] * delta_TM_ini + dlogC_dlogRLya[i]*(theta_b[i]/(3.*Hubble[i])-delta_b[i]))
	print i, ddelta_xe_dz, ddelta_TM_dz
	delta_xe_ini -= ddelta_xe_dz*dz
	delta_xe.append (delta_xe_ini)
plt.figure(1)
plt.xscale('log')
#plt.xscale('log', basex=2)
plt.plot(1+new_z, delta_TM*(1+new_z), label = r'With $\delta_{x_e}$')
#plt.axis([1,10**3,-0.0001,0.0008])
plt.xlabel (r'$1+z$', size = 15)
plt.ylabel (r'$\delta_{T_\mathrm{gas}} (1+z)$', size = 15)
fig.set_xticks([5,10,50,100,500])
fig.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.legend ()
#plt.savefig('deltaTgas_k1e-1.pdf')



plt.figure(2)
fig = plt.subplot(1,1,1)
#plt.xscale('log', basex=2)
plt.xscale('log')
fig.set_xticks([5,10,50,100,500])
fig.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.xlabel (r'$1+z$', size = 15)
plt.ylabel (r'$\delta_{x_e}$', size = 15)
#plt.yscale('log')
plt.plot(1+new_z, delta_xe)
#plt.savefig('deltaxe_k1e-1.pdf')
plt.show()







#t4 = Tm/kBoltz/1.e4
#alphaB_PPB = 4.309e-13*pow(t4,-0.6166)/(1.+ 0.6703*pow(t4,0.5300))
