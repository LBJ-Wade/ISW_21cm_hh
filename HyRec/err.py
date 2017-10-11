#from D import *
#from isw_transfer import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator
import data
import data_no_nu
index = [1,2,4,5]
l = 1000
l_list = np.arange(2, l+1)
number_of_l = 4998

infile = data_no_nu.cl
ll = np.loadtxt(infile)[0:number_of_l,0]
aa = ll*(ll+1)/(2*np.pi) /( 10**6 * 2.7255)
aa = 1/(10**6*2.7255**2)
cl_TT_nonu = np.loadtxt(infile)[0:number_of_l,1]/aa
cl_TE_nonu = np.loadtxt(infile)[0:number_of_l,3]/aa
cl_EE_nonu = np.loadtxt(infile)[0:number_of_l,2]/aa

infile = data.cl
cl_TT = np.loadtxt(infile)[0:number_of_l,1]/aa
cl_TE = np.loadtxt(infile)[0:number_of_l,3]/aa
cl_EE = np.loadtxt(infile)[0:number_of_l,2]/aa


cl21T_nonu = []
cl21T_nonu_Neff = []
cl21T_nu1 = []
cl21T_nu2 = []
cl21_nonu = []
cl21_nonu_Neff = []
cl21_nu1 = []
cl21_nu2 = []
n = 8

l_list = np.loadtxt('result/antony/nanoom/21T_nonu_30.txt')[0:,0]
print (len(l_list))
aaa = 10**6* l_list*(l_list+1)/(2*np.pi)

z_m_list = [30,50,75,100,125,150,175,200]
for i in range(n):
	nonu = np.loadtxt('result/antony/nanoom/21T_nonu_{}.txt'.format(z_m_list[i]))[0:,1]*aaa
	nonu_Neff = np.loadtxt('result/antony/nanoom/21T_nonu_Neff_{}.txt'.format(z_m_list[i]))[0:,1]*aaa
	nu1 = np.loadtxt('result/antony/nanoom/21T_nu1_{}.txt'.format(z_m_list[i]))[0:,1]*aaa
	#nu1 = np.loadtxt('result/antony/nanoom/21T_nonu_Neff_{}.txt'.format(z_m_list[i]))[0:,1]*aaa
	nu2 = np.loadtxt('result/antony/nanoom/21T_nu2_{}.txt'.format(z_m_list[i]))[0:,1]*aaa
	nonu21 = np.loadtxt('result/antony/nanoom/21_nonu_{}.txt'.format(z_m_list[i]))[0:,1]*aaa
	nonu_Neff21 = np.loadtxt('result/antony/nanoom/21_nonu_Neff_{}.txt'.format(z_m_list[i]))[0:,1]*aaa
	nu21_1 = np.loadtxt('result/antony/nanoom/21_nu1_{}.txt'.format(z_m_list[i]))[0:,1]*aaa
	#nu21_1 = np.loadtxt('result/antony/nanoom/21_nonu_Neff_{}.txt'.format(z_m_list[i]))[0:,1]*aaa
	nu21_2 = np.loadtxt('result/antony/nanoom/21_nu2_{}.txt'.format(z_m_list[i]))[0:,1]*aaa

	nonu = np.interp (ll, l_list, nonu)
	nonu_Neff = np.interp (ll, l_list, nonu_Neff)
	nu1 = np.interp (ll, l_list, nu1)
	nu2 = np.interp (ll, l_list, nu2)
	nonu21 = np.interp (ll, l_list, nonu21)
	nonu_Neff21 = np.interp (ll, l_list, nonu_Neff21)
	nu21_1 = np.interp (ll, l_list, nu21_1)
	nu21_2 = np.interp (ll, l_list, nu21_2)

	cl21T_nonu.append (nonu)
	cl21T_nonu_Neff.append (nonu_Neff)
	cl21T_nu1.append (nu1)
	cl21T_nu2.append (nu2)
	cl21_nonu.append (nonu21)
	cl21_nonu_Neff.append (nonu_Neff21)
	cl21_nu1.append (nu21_1)
	cl21_nu2.append (nu21_2)
print (len(cl_TT), len(cl21T_nonu[0]))

#cl_TT = np.interp(l_list, ll, cl_TT)
#cl_TE = np.interp(l_list, ll, cl_TE)
#cl_EE = np.interp(l_list, ll, cl_EE)
#cl_TT_nonu = np.interp(l_list, ll, cl_TT_nonu)
#cl_TE_nonu = np.interp(l_list, ll, cl_TE_nonu)
#cl_EE_nonu = np.interp(l_list, ll, cl_EE_nonu)
l_list = ll
cl_TT_N = (1*0.000290888*10**-3)**2 *np.e**(ll*(ll+1)*(10**-3*0.000290888)**2/(8*np.log(2)))
cl_TT_N = cl_TT_N * ll*(ll+1)/(2*np.pi)

#plt.axis([0,5000,0,6000])
#plt.xlabel (r'$l$')
#plt.ylabel (r'$l(l+1)C_l^{TT}$ [$\mu K ^2$]')
infile = 'clTT_sevem.txt'
ell = np.loadtxt (infile)[0:,0]
cl = np.loadtxt (infile)[0:,1]/10**6*10**12
plt.figure(100)
plt.plot(l_list, np.sqrt(cl_TT), label = 'clTT')
#plt.plot(ell, np.sqrt(cl), label = 'Planck + interp')
#plt.plot(l_list, np.sqrt(cl21T_nonu[0]))
plt.plot(l_list, np.sqrt(cl_TT_N), label = 'noise')
plt.yscale('log')
#plt.plot(l_list, cl_TT*10**6,'r', label = r'$m_{\nu}=0$eV')
#plt.legend ()
#plt.savefig ('clTT_interp.pdf',format='pdf')

plt.show()
cl_TT = cl
cl_TT_nonu = cl

#z_m = [1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]
z_m = [1,2,4,8,10,15,20,30,40,60,80,100]
l_bin_index = [2,4,7,13,27,50,91,190,320,620,1000]
l_bin_index = [2,4,7,13,26,50,91,190,320,620,1000]
l_bin_index = [2,9,43,204,956,4472,20912,97793,457305,2138469,10000000]
l_bin_index = [2, 4, 9, 20, 45, 100, 218, 478, 1045, 2286, 5000]
#l_bin_index = [2,5,12,31,79,200,502,1261,3169,7962,20000]
#l_bin_index = [2,4,8,16,34,70,144,294,600,1225,2500]
#fig, axs = plt.subplots(nrows=4, ncols=4, sharex=True)
plt.figure(1)
for i in range(9):
	fig = plt.subplot(3,3,i+1)
	if i == 8:
		plt.plot (l_list, 10**10*np.sqrt(cl21T_nonu[0]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0$eV')
		plt.plot (l_list, 10**10*np.sqrt(cl21T_nu1[0]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0.1$eV')
		plt.plot (l_list, 10**10*np.sqrt(cl21T_nu2[0]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0.2$eV')
		plt.plot (l_list, 10**10*np.sqrt(cl21T_nonu_Neff[0]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0$eV, $\mathrm{N}_{\mathrm{eff}}=3.2$')
		plt.legend(prop={'size':10})	
	
	else:
		#plt.plot (l_list, np.sqrt(cl21T_nonu[i]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0$eV')
		#plt.plot (l_list, np.sqrt(cl21T_nu1[i]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0.1$eV')
		#plt.plot (l_list, np.sqrt(cl21T_nu2[i]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0.2$eV')
		plt.plot (l_list, np.sqrt(cl21T_nonu[i]), label = r'$m_{\nu,i}=0$eV')
		plt.plot (l_list, np.sqrt(cl21T_nu1[i]), label = r'$m_{\nu,i}=0.1$eV')
		plt.plot (l_list, np.sqrt(cl21T_nu2[i]), label = r'$m_{\nu,i}=0.2$eV')
		plt.plot (l_list, np.sqrt(cl21T_nonu_Neff[i]), label = r'$m_{\nu,i}=0$eV, $\mathrm{N}_{\mathrm{eff}}=3.2$')
		plt.text(30,  3*10**-9, r'z$_m=${}'.format(z_m_list[i])) 
	plt.xscale ('log')
	plt.yscale ('log')
	plt.axis([2,10**7,10**-9,10**-2])
	if i not in [0,3,6]:	
		fig.axes.get_yaxis().set_visible(False)
	if i not in [6,7,8]:
		fig.axes.get_xaxis().set_visible(False)
	
	if i in [0,3,6]:
		plt.ylabel (r'$[l(l+1)C_l/2\pi]^{1/2}$ mK', size = 10)
	if i in [6,7,8]:
		plt.xlabel (r'$l$')
	if i == 0:
		fig.set_yticks([10**-8,10**-5,10**-2])
	if i in [3,6]:
		fig.set_yticks([10**-8,10**-5, 10**-2])
	
	if i == 8:
		fig.set_xticks([10, 10**3,10**5,10**7])
	if i in [6,7]:
		fig.set_xticks([10, 10**3,10**5])
	
	plt.subplots_adjust(wspace=0, hspace=0)
#plt.savefig ('cl21T.pdf',format='pdf')
plt.suptitle (r'$C_l^{21,\mathrm{ISW}}$',size = 15)

for i in [0,7]:#range(2):
	plt.figure(i+10)
	#plt.plot (l_list, np.sqrt(cl21T_nonu[i]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0$eV')
	#plt.plot (l_list, np.sqrt(cl21T_nu1[i]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0.1$eV')
	#plt.plot (l_list, np.sqrt(cl21T_nu2[i]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0.2$eV')
	plt.plot (l_list, np.sqrt(cl21T_nonu[i]), label = r'$m_{\nu,i}=0$eV')
	plt.plot (l_list, np.sqrt(cl21T_nu1[i]), label = r'$m_{\nu,i}=0.1$eV')
	plt.plot (l_list, np.sqrt(cl21T_nu2[i]), label = r'$m_{\nu,i}=0.2$eV')
	plt.plot (l_list, np.sqrt(cl21T_nonu_Neff[i]), label = r'$m_{\nu,i}=0$eV, $\mathrm{N}_{\mathrm{eff}}=3.2$')
	plt.text(30,  3*10**-9, r'z$_m=${}'.format(z_m_list[i])) 
	plt.xscale ('log')
	plt.yscale ('log')
	plt.axis([2,10**7,10**-9,10**-2])
	plt.legend (prop={'size':12})	
	plt.ylabel (r'$[l(l+1)C_l/2\pi]^{1/2}$ mK')
	plt.xlabel (r'$l$')
	fig.set_yticks([10**-8,10**-5,10**-2])
	fig.set_xticks([10, 10**3,10**5,10**7])
	#plt.savefig ('cl21T_{}.pdf'.format(z_m_list[i]), format='pdf')

plt.figure(2)
for i in range(9):
	fig = plt.subplot(3,3,i+1)
	if i == 8:
		plt.plot (l_list, 10**10*np.sqrt(cl21_nonu[0]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0$eV')
		plt.plot (l_list, 10**10*np.sqrt(cl21_nu1[0]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0.1$eV')
		plt.plot (l_list, 10**10*np.sqrt(cl21_nu2[0]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0.2$eV')
		plt.plot (l_list, 10**10*np.sqrt(cl21T_nonu_Neff[0]*l_list*(l_list+1)/(2*np.pi)), label = r'$m_{\nu,i}=0$eV, $\mathrm{N}_{\mathrm{eff}}=3.2$')
		plt.legend(prop={'size':10})	
	
	else:
		plt.plot (l_list, np.sqrt(cl21_nonu[i]), label = r'$m_{\nu,i}=0$eV')
		plt.plot (l_list, np.sqrt(cl21_nu1[i]), label = r'$m_{\nu,i}=0.1$eV')
		plt.plot (l_list, np.sqrt(cl21_nu2[i]), label = r'$m_{\nu,i}=0.2$eV')
		plt.plot (l_list, np.sqrt(cl21_nonu_Neff[i]), label = r'$m_{\nu,i}=0$eV, $\mathrm{N}_{\mathrm{eff}}=3.2$')
		plt.text(30,  3*10**-4, r'z$_m=${}'.format(z_m_list[i])) 
	plt.xscale ('log')
	plt.yscale ('log')
	plt.axis([2,10**7,10**-4,1])
	
	if i not in [0,3,6]:	
		fig.axes.get_yaxis().set_visible(False)
	if i not in [6,7,8]:
		fig.axes.get_xaxis().set_visible(False)


	if i in [0,3,6]:
		plt.ylabel (r'$[l(l+1)C_l/2\pi]^{1/2}$ mK', size = 10)
	if i in [6,7,8]:
		plt.xlabel (r'$l$')
	if i == 0:
		fig.set_yticks([10**-4,10**-2,1])
	if i in [3,6]:
		fig.set_yticks([10**-4,10**-2])
	
	if i == 8:
		fig.set_xticks([10, 10**3,10**5,10**7])
	if i in [6,7]:
		fig.set_xticks([10, 10**3,10**5])
	
	plt.subplots_adjust(wspace=0, hspace=0)
plt.suptitle (r'$C_l^{21,21}$',size = 15)
#plt.savefig ('cl21.pdf',format='pdf')

s_to_n = []
s_to_n_TE = []
s_to_n_nonu = []
s_to_n_nonu_TE = []
plt.figure(3)
for i in range(8):
	fig = plt.subplot(3,3,i+1)
	l_bin_yerr1 = []
	yy1 = []
	xx1 = []
	l_bin_yerr2 = []
	yy2 = []
	xx2 = []
	
	f_sky = 0.65
	f_sky = 1
	print (len(cl_TT), len(cl21_nu1[i]), len(l_list))
	yerr1 = np.sqrt ((cl_TT*(cl21_nu1[i]) + cl21T_nu1[i]**2)/(2*l_list+1)/f_sky)
	
	#s_to_n.append(sum(cl21T_nu1[i]**2/yerr1**2))
	s_to_n.append (sum(np.interp (np.arange(2,5001), l_list, cl21T_nu1[i]**2/yerr1**2)))
	yerr1_nonu = np.sqrt ((cl_TT_nonu*(cl21_nonu[i]) + cl21T_nonu[i]**2)/(2*l_list+1)/f_sky)
	#s_to_n_nonu.append(sum(cl21T_nonu[i]**2/yerr1_nonu**2))
	s_to_n_nonu.append (sum(np.interp (np.arange(2,5001), l_list, cl21T_nonu[i]**2/yerr1_nonu**2)))
	
	f_sky = 0.25
	f_sky = 1
	yerr2 = np.sqrt (((cl_TT-cl_TE**2/cl_EE)*(cl21_nu1[i]) + cl21T_nu1[i]**2)/(2*l_list+1)/f_sky)
	s_to_n_TE.append(sum(cl21T_nu1[i]**2/yerr2**2))

	yerr2_nonu = np.sqrt (((cl_TT-cl_TE**2/cl_EE)*(cl21_nonu[i]) + cl21T_nonu[i]**2)/(2*l_list+1)/f_sky)
	s_to_n_nonu_TE.append(sum(cl21T_nonu[i]**2/yerr2_nonu**2))
	w1 = 1/yerr1**2 
	w2 = 1/yerr2**2
	if i == 0:
		print (yerr1)
		print (cl_TT[0]-cl_TE[0]**2/cl_EE[0])
	real_var1 = []
	real_var2 = []
	interp_gT = interp1d (l_list, cl21T_nu1[i])
	for j in range(len(l_bin_index)-1):
		l_bin = []
			
		for mm in range(len(l_list)):
			if j == 0:
				if l_list[mm] < (l_bin_index[j+1]+1):
					l_bin.append (l_list[mm])
			else:
				if l_list[mm] < (l_bin_index[j+1]+1) and l_list[mm] > (l_bin_index[j]-1):
					l_bin.append (l_list[mm])
		print (j, len(l_bin), (l_bin[-1]-l_bin[0]))
		l_bin = np.array (l_bin)
		
		#l_bin = np.arange(l_bin_index[j], l_bin_index[j+1]+1)
		
		ww1 = 0
		ww2 = 0
		print (l_bin[-1])
		#for k in l_bin[:-1]:
		for k in np.arange(l_bin[:-1][0],l_bin[:-1][-1]+1):
			ww1 += np.interp (k, l_list, w1)
			ww2 += np.interp (k, l_list, w2)
			#ww1 += w1[k-2]
			#ww2 += w2[k-2]

		#ww1 = ww1 *((l_bin[-1]-l_bin[0])/len(l_bin-1))
		#ww2 = ww2 *((l_bin[-1]-l_bin[0])/len(l_bin-1))
		
		real_var1.append (np.sqrt(1/ww1))
		real_var2.append (np.sqrt(1/ww2))
	real_var1 = np.array (real_var1)
	real_var2 = np.array (real_var2)
	
	interp_gT = interp1d (l_list, cl21T_nu1[i])
	for j in range(len(l_bin_index)-1):
		yerr_bin = 0
		l_bin = []
			
		for mm in range(len(l_list)):
			if j == 0:
				if l_list[mm] < (l_bin_index[j+1]+1):
					l_bin.append (l_list[mm])
			else:
				if l_list[mm] < (l_bin_index[j+1]+1) and l_list[mm] > (l_bin_index[j]-1):
					l_bin.append (l_list[mm])
		l_bin = np.array (l_bin)
		
		#l_bin = np.arange(l_bin_index[j], l_bin_index[j+1]+1)
		
		l_mid = 10**((np.log10(l_bin[0]) + np.log10(l_bin[-1]))/2)
		#l_bin_yerr1.append (real_var1[j])
		#l_bin_yerr2.append (real_var2[j])
		l_bin_yerr1 = real_var1
		l_bin_yerr2 = real_var2
		xx1.append (l_mid)
		yy1.append (interp_gT(l_mid))
		xx2.append (l_mid)
		yy2.append (interp_gT(l_mid))
	plt.xscale('log')
	plt.yscale('log')
	xx1 = np.array(xx1)
	yy1 = np.array(yy1)
	l_bin_yerr1 = np.array (l_bin_yerr1)
	l_bin_yerr2 = np.array (l_bin_yerr2)

	yy1 = np.sqrt(yy1)
	l_bin_yerr1 = np.sqrt(l_bin_yerr1)
	l_bin_yerr2 = np.sqrt(l_bin_yerr2)
	cl21T_nu1[i] = np.sqrt(cl21T_nu1[i])
	cl21T_nonu[i] = np.sqrt(cl21T_nonu[i])

	plt.errorbar(xx1, yy1, yerr=l_bin_yerr1, fmt = None, ecolor='r')
	plt.errorbar(xx1, yy1, yerr=l_bin_yerr2, fmt = None, ecolor='b')
	plt.plot(l_list, cl21T_nu1[i], color = 'k')
	plt.plot(l_list, cl21T_nonu[i], '--k')
	#plt.axis([2,10**7,10**-9, 10**-2])
	plt.axis([2,5000,10**-5, 10**-2])
	plt.text(30,  3*10**-5, r'z$_m=${}'.format(z_m_list[i]))	
	#plt.axis([2,2*10**4,3*10**-5, 10**-2])
	#plt.text(30,  5*10**-5, r'z$_m=${}'.format(z_m_list[i]))	
	"""
	if i in [6,7,8]:
		plt.xlabel (r'$l$')
		if not i == 7:
			fig.set_xticks([10**1,10**3,10**5])
		else:
			fig.set_xticks([10**1,10**3,10**5,10**7])
	if i == 5:
		fig.set_xticks([10**1,10**3,10**5,10**7])
	if i == 0:
		fig.set_yticks([10**-8,10**-5,10**-2])
	if i != 0:
		fig.set_yticks([10**-8,10**-5])
	"""
	
	
	if i in [6,7,8]:
		plt.xlabel (r'$l$')
		if not i == 7:
			fig.set_xticks([10**1,10**2,10**3])
		else:
			fig.set_xticks([10**1,10**2,10**3])
	if i == 5:
		fig.set_xticks([10**1,10**2,10**3])
	if i == 0:
		fig.set_yticks([10**-4,10**-3,10**-2])
	if i != 0:
		fig.set_yticks([10**-4,10**-3])
		#fig.set_yticks([10**-5,10**-4,10**-3])
	if i in [0,3,6]:
		plt.ylabel (r'$[l(l+1)C_l/2\pi]^{1/2}$ mK', size = 10)
	if i not in [5,6,7]:
		fig.axes.get_xaxis().set_visible(False)
	if i not in [0,3,6]:
		fig.axes.get_yaxis().set_visible(False)
	#plt.locator_params(axis='x', nticks=4)
	plt.tick_params(labelsize=11)
	
	plt.subplots_adjust(wspace=0, hspace=0)
#plt.suptitle (r'$C_l^{21,\mathrm{ISW}}$, $m_{\nu, i} = 0\mathrm{eV}$, $\mathrm{N}_{\mathrm{eff}}=3.2$',size = 15)
plt.suptitle (r'$C_l^{21,\mathrm{ISW}}$, $m_{\nu, i} = 0.1\mathrm{eV}$',size = 15)
plt.savefig('err_cl21T_nu1.pdf',format='pdf')
#plt.savefig('err_cl21T_nonu_Neff.pdf',format='pdf')

plt.figure(4)
plt.plot(z_m_list, s_to_n_nonu, label = r'$m_{\nu}=0$eV')
plt.plot(z_m_list, s_to_n, label = r'$m_{\nu}=0.1$eV')
#plt.plot(z_m_list,[1.1217755862667766, 3.0468472146985195, 7.049576545797926, 13.304004517145021, 22.135147942892807, 33.676234190087641, 47.971560407080993, 65.231561221914916], label = r'$m_{\nu}=0.2$eV')
#plt.plot(z_m_list,[0.307455407322478, 1.4005556951122922, 4.3625313435401525, 9.580029883700405, 17.399184175674836, 27.962835373719969, 41.284964281147985, 57.490515401043972], label = r'$m_{\nu}=0$eV, $\mathrm{N}_{\mathrm{eff}}=3.2$')
#plt.plot(z_m_list, [0.79560796753887286, 2.9110498874214477, 8.4238642210555028, 18.120741689720099, 32.527323546931981, 51.537798983531644, 74.668364138010062, 101.61733722763323], label = r'$m_{\nu}=0.2$eV')
#plt.plot(z_m_list, [0.61780600904851457, 2.6946376667916478, 8.0560942832535147, 16.745382909165823, 28.353969954731703, 41.976397215955231, 56.596616267157074, 71.56972159537888],label = r'$m_{\nu,i}=0$eV, $\mathrm{N}_{\mathrm{eff}}=3.2$') 

#plt.yscale('log')
#plt.xscale('log')

#plt.plot(z_m_list, s_to_n_TE, label = 'nu TE')
#plt.plot(z_m_list, s_to_n_nonu_TE, label = 'nonu TE')

plt.grid ()
plt.ylabel (r'$(S/N)^2$')
plt.xlabel (r'$z_m$')
#plt.axis([0,200,0,450])
#plt.axis([0,30,0,18])
plt.legend(bbox_to_anchor=(0.4,1),prop={'size':12})
plt.title (r'Signal to noise of $C_l^{21,\mathrm{ISW}}$')
#plt.axis([0,100,0,180])
plt.savefig ('SN_cl21T.pdf', format = 'pdf')

print (s_to_n)

"""
plt.figure(123)
for i in range(12):
	fig = plt.subplot(3,4,i+1)
	plt.plot(l_list, limber_no_nu_tau[i]/limber_no_nu[i])
	plt.axis([2,1000,0.9,1])
	if i in [11,10,9,8]:
		plt.xlabel (r'$l$')
	if i not in [11,10,9,8]:
		fig.axes.get_xaxis().set_visible(False)
	if i not in [0,4,8]:
		fig.axes.get_yaxis().set_visible(False)
	plt.tick_params(labelsize=11)
	if i == 0:
		fig.set_yticks([0.9,0.95,1])
	if i != 0:
		fig.set_yticks([0.9,0.95])
	if i == 11:
		fig.set_xticks([200,400,600,800,1000])
	if i in [8,9,10]:
		fig.set_xticks([200,400,600,800])
	
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.text(400,  0.91, r'z$_m=${}'.format(z_m[i]))	
plt.suptitle(r'$C_l^{gT}(\tau)/C_l^{gT}(\tau=0)$', size = 15)
plt.savefig('clgT_tau_notau.pdf',format='pdf')
"""

"""
plt.figure(1000)
for i in range(18):
	fig = plt.subplot(5,4,i+1)
	plt.plot(l_list, limber_nu_3[i]**2*(2*l_list+1)/(limber_gg_3[i]*cl_TT + limber_nu_3[i]**2), label = 'w/o TE') 
	#plt.plot(l_list, limber_nu_3[i]**2*(2*l_list+1)/(limber_gg_3[i]*(cl_TT-cl_TE**2/cl_EE) + limber_nu_3[i]**2), label = 'w/ TE') 
	#print (sum( limber_nu_3[i]**2*(2*l_list+1)/(limber_gg_3[i]*cl_TT + limber_nu_3[i]**2) ), sum( (limber_nu_3[i]**2*(2*l_list+1)/(limber_gg_3[i]*cl_TT + limber_nu_3[i]**2))[60:]))
	#plt.plot(l_list, limber_no_nu[i]**2*(2*l_list+1)/(limber_gg[i]*cl_TT_nonu + limber_no_nu[i]**2), label = 'w/o TE') 
	#plt.plot(l_list, limber_no_nu[i]**2*(2*l_list+1)/(limber_gg[i]*(cl_TT_nonu-cl_TE_nonu**2/cl_EE_nonu) + limber_no_nu[i]**2), label = 'w/ TE') 
	#print (sum( limber_no_nu[i]**2*(2*l_list+1)/(limber_gg[i]*cl_TT_nonu + limber_no_nu[i]**2) ), sum( (limber_no_nu[i]**2*(2*l_list+1)/(limber_gg[i]*cl_TT_nonu + limber_no_nu[i]**2))[20:]))
	#plt.legend()
	plt.axis([2,1000,0,0.03])
	if i in [14,15,16,17]:
		plt.xlabel (r'$l$')
	
	if i not in [14,15,16,17]:
		fig.axes.get_xaxis().set_visible(False)
	if i not in [0,4,8,12,16]:
		fig.axes.get_yaxis().set_visible(False)
	#plt.locator_params(axis='x', nticks=4)
	plt.tick_params(labelsize=11)
	
	if i == 0:
		fig.set_yticks([0,0.01,0.02,0.03])
	if i != 0:
		fig.set_yticks([0,0.01,0.02])
	plt.text(400,  0.025, r'z$_m=${}'.format(z_m[i]))	
	plt.subplots_adjust(wspace=0, hspace=0)
plt.suptitle(r'$(S_l/N_l)^2$, $f_{\nu}=0$')
#plt.savefig('sn_nonu.pdf',format='pdf')

plt.figure(1001)
for i in range(n):
	fig = plt.subplot(5,4,i+1)
	plt.plot(l_list, limber_no_nu[i], label = '{}'.format(i)) 
var_gT = []
var_gg = []
cov = []
for i in range(n):
	var_gT.append(1/(2*l_list+1)*(cl_TT*limber_gg_3[i] + limber_nu_3[i]**2))
	var_gg.append(2/(2*l_list+1)*limber_gg_3[i]**2)
	cov.append(2/(2*l_list+1)*limber_nu_3[i]*limber_gg_3[i])

expectation = []
var = []
for i in range(n):
	expectation.append (limber_nu_3[i]/limber_gg_3[i] - cov[i]/limber_gg_3[i]**2 + var_gg[i]*limber_nu_3[i]/limber_gg_3[i]**3)
	var.append (limber_nu_3[i]**2/limber_gg_3[i]**2 * (var_gT[i]/limber_nu_3[i]**2 - 2*cov[i]/(limber_nu_3[i]*limber_gg_3[i]) + var_gg[i]/limber_gg_3[i]**2))

sn = []
for i in range(n):
	#w = 1/var[i]
	sn.append(sum ( (limber_nu_3[i]/limber_gg_3[i] - limber_no_nu[i]/limber_gg[i])**2 / var[i]))


z_m = [1,2,4]

plt.figure(3)
for i in range(3):
	fig = plt.subplot(1,3,i+1)
	plt.xscale('log')
	plt.yscale('log')
	
	plt.plot(abs(limber_nu_3[i]*limber_gg_3[i]-limber_no_nu[i]*limber_gg[i]), label = 'nu')
	#plt.axis([2,1000,0,9*10**-12])
	plt.axis([2,1000,10**-19,10**-11])
	
	#if i == 0:
	#	plt.ylabel(r'$|C_l^{gT}(f_{\nu}=0.1)C_l^{gg}(f_{\nu}=0.1)-C_l^{gT}(f_{\nu}=0)C_l^{gg}(f_{\nu}=0)|$  $\times (l(l+1)/2\pi)^2$')
	if i == 1:
		plt.title(r'$|C_l^{gT}(f_{\nu}=0.1)C_l^{gg}(f_{\nu}=0.1)-C_l^{gT}(f_{\nu}=0)C_l^{gg}(f_{\nu}=0)|$  $\times (l(l+1)/2\pi)^2$')
	plt.xlabel (r'$l$')
	plt.text(30,  10**-18, r'z$_m=${}'.format(z_m[i]))	
	if i in [1,2]:
		fig.axes.get_yaxis().set_visible(False)
plt.savefig ('gTgg_diff.pdf', format = 'pdf')




plt.figure(4)
for i in range(3):
	fig = plt.subplot(1,3,i+1)
	plt.xscale('log')
	plt.yscale('log')

	plt.plot(limber_nu_3[i]*limber_gg_3[i], label = r'$f_{\nu}=0.1$')
	plt.plot(limber_no_nu[i]*limber_gg[i], label = r'$f_{\nu}=0$')
	#plt.axis([2,1000,0,9*10**-12])
	plt.axis([2,1000,10**-17,10**-11])
	#if i == 0:
	#	plt.ylabel(r'$C_l^{gT}(f_{\nu})C_l^{gg}(f_{\nu})$ $\times (l(l+1)/2\pi)^2$')
	if i == 1:
		plt.title(r'$C_l^{gT}(f_{\nu})C_l^{gg}(f_{\nu})$ $\times (l(l+1)/2\pi)^2$')
	plt.xlabel (r'$l$')
	plt.text(30,  2*10**-17, r'z$_m=${}'.format(z_m[i]))	
	if i in [1,2]:
		fig.axes.get_yaxis().set_visible(False)
plt.legend ()
plt.savefig ('gTgg.pdf', format = 'pdf')

plt.figure(5)
for i in range(3):
	fig = plt.subplot(1,3,i+1)
	plt.xscale('log')
	plt.yscale('log')
	plt.plot(limber_nu_3[i], label = r'$f_{\nu}=0.1$')
	plt.plot(limber_no_nu[i], label = r'$f_{\nu}=0$')
	#plt.axis([2,1000,0,1.4*10**-8])
	plt.axis([2,1000,10**-10,2*10**-8])
	#if i == 0:
	#	plt.ylabel(r'$C_l^{gT}(f_{\nu})$ $\times l(l+1)/2\pi$')
	if i == 1:
		plt.title(r'$C_l^{gT}(f_{\nu})$ $\times l(l+1)/2\pi$')

	plt.xlabel (r'$l$')
	plt.text(30,  1.5*10**-10, r'z$_m=${}'.format(z_m[i]))	
	if i in [1,2]:
		fig.axes.get_yaxis().set_visible(False)
	
plt.legend ()
plt.savefig ('gT.pdf', format = 'pdf')

plt.figure(6)
for i in range(3):
	fig = plt.subplot(1,3,i+1)
	plt.xscale('log')
	plt.yscale('log')
	plt.plot(limber_gg_3[i], label = r'$f_{\nu}=0.1$')
	plt.plot(limber_gg[i], label = r'$f_{\nu}=0$')
	#plt.axis([2,1000,0,0.0035])
	plt.axis([2,1000,10**-6,10**-2])
	#if i == 0:
	#	plt.ylabel(r'$C_l^{gg}(f_{\nu})$ $\times l(l+1)/2\pi$')
	if i == 1:
		plt.title(r'$C_l^{gg}(f_{\nu})$ $\times l(l+1)/2\pi$')

	plt.xlabel (r'$l$')
	plt.text(30,  2*10**-6, r'z$_m=${}'.format(z_m[i]))	
	if i in [1,2]:
		fig.axes.get_yaxis().set_visible(False)
plt.legend ()
plt.savefig ('gg.pdf', format = 'pdf')

plt.figure(7)
for i in range(9):
	fig = plt.subplot(3,3,i+1)
	#plt.xscale('log')
	#plt.yscale('log')
	plt.plot(l_list, limber_nu_3[i]/limber_gg_3[i])
	plt.plot(l_list, limber_no_nu[i]/limber_gg[i])
	#plt.axis([2,1000,0,0.0035])
	#plt.axis([2,10,0,0.002])
	#plt.axis([2,1000,10**-6,10**-2])
	#if i == 0:
	#	plt.ylabel(r'$C_l^{gg}(f_{\nu})$ $\times l(l+1)/2\pi$')
	if i == 1:
		plt.title(r'$C_l^{gg}(f_{\nu})$ $\times l(l+1)/2\pi$')

	plt.xlabel (r'$l$')
	#plt.text(30,  2*10**-6, r'z$_m=${}'.format(z_m[i]))	
	if i in [1,2]:
		fig.axes.get_yaxis().set_visible(False)
plt.savefig ('gg.pdf', format = 'pdf')

plt.figure(70)
plt.loglog(l_list, cl_TE**2/cl_EE, label = 'nu')
plt.loglog(l_list, cl_TE_nonu**2/cl_EE_nonu, label = 'nonu')
plt.grid()
plt.legend()

z_m = [1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,50,100]
plt.figure(71)
for i in range(18):
	plt.loglog(l_list, abs(limber_no_nu[i]*aa), label= '{}'.format(z_m[i]))
plt.figure(72)
for i in range(18):
	plt.loglog(l_list, abs(limber_nu_3[i]*aa), label = '{}'.format(z_m[i]))

plt.figure(73)
plt.loglog(l_list, cl_TT)
"""
plt.show()

