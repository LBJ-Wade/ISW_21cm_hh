#from D import *
#from isw_transfer import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator
aa = 1/(10**6*2.7255**2)
#cl_TT_nonu = np.loadtxt(infile)[0:number_of_l,1]/aa
#cl_TE_nonu = np.loadtxt(infile)[0:number_of_l,3]/aa
#cl_EE_nonu = np.loadtxt(infile)[0:number_of_l,2]/aa

infile = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/data/cl_0_BBN.dat"
cl_TT = np.loadtxt(infile)[0:2999,1]/10**6
cl_TE = np.loadtxt(infile)[0:2999,3]/10**6
cl_EE = np.loadtxt(infile)[0:2999,2]/10**6


cl21T_nonu = []
cl21T_nonu_Neff = []
cl21T_nu1 = []
cl21T_nu2 = []
cl21_nonu = []
cl21_nonu_Neff = []
cl21_nu1 = []
cl21_nu2 = []
n = 8


result_dir = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/result_Yp_BBN/"
infile1 = result_dir + "cl21T_0.txt" 
z_m_list = [30,50,75,100,125,150,175,200]
cl21T = []
l_list = np.loadtxt (infile1) [0:2999,0]
aaa = 10**6* l_list*(l_list+1)/(2*np.pi) * 2.7255
for i in range(n):
	#cl21T_nonu.append ()
	cl21T_nu1.append (np.loadtxt (infile1)[0:2999,i+1]*aaa)
	#cl21_nonu.append (nonu21)
	cl21_nu1.append (np.loadtxt (result_dir + "cl21zmzl_{0}{1}.txt".format(i,i))[0:2999,1]*aaa)
ll = l_list
cl_TT_N = (1*0.000290888*10**-3)**2 *np.e**(ll*(ll+1)*(10**-3*1*0.000290888)**2/(8*np.log(2)))

cl_TT_N = cl_TT_N * ll*(ll+1)/(2*np.pi)

cl_TT += cl_TT_N
#cl_TT_nonu += cl_TT_N
#z_m = [1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]
z_m = [1,2,4,8,10,15,20,30,40,60,80,100]
l_bin_index = [2,4,7,13,27,50,91,190,320,620,1000]
l_bin_index = [2,4,7,13,26,50,91,190,320,620,1000]
l_bin_index = [2,9,43,204,956,4472,20912,97793,457305,2138469,10000000]
l_bin_index = [2, 4, 9, 20, 45, 100, 218, 478, 1045, 2286, 5000]
l_bin_index = [2, 4, 9, 18, 37, 77, 161, 334, 695, 1444, 3000]
#l_bin_index = [2,5,12,31,79,200,502,1261,3169,7962,20000]
#l_bin_index = [2,4,8,16,34,70,144,294,600,1225,2500]
#fig, axs = plt.subplots(nrows=4, ncols=4, sharex=True)

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
	f_sky = 0.7
	print (len(cl_TT), len(cl21_nu1[i]), len(l_list))
	yerr1 = np.sqrt ((cl_TT*(cl21_nu1[i]) + cl21T_nu1[i]**2)/(2*l_list+1)/f_sky)
	w1 = 1/yerr1**2 
	
	s_to_n.append (sum(np.interp (np.arange(2,3001), l_list, cl21T_nu1[i]**2/yerr1**2)))
	#yerr1_nonu = np.sqrt ((cl_TT_nonu*(cl21_nonu[i]) + cl21T_nonu[i]**2)/(2*l_list+1)/f_sky)
	#s_to_n_nonu.append (sum(np.interp (np.arange(2,3001), l_list, cl21T_nonu[i]**2/yerr1_nonu**2)))

	
	f_sky = 0.25
	f_sky = 0.4
	yerr2 = np.sqrt ((cl_TT*(cl21_nu1[i]) + cl21T_nu1[i]**2)/(2*l_list+1)/f_sky)
	s_to_n_TE.append(sum(np.interp (np.arange(2,3001), l_list, cl21T_nu1[i]**2/yerr2**2)))
	
	#yerr2_nonu = np.sqrt (((cl_TT-cl_TE**2/cl_EE)*(cl21_nonu[i]) + cl21T_nonu[i]**2)/(2*l_list+1)/f_sky)
	#yerr2_nonu = np.sqrt ((cl_TT_nonu*(cl21_nonu[i]) + cl21T_nonu[i]**2)/(2*l_list+1)/f_sky)
	#s_to_n_nonu_TE.append(sum(np.interp (np.arange(2,3001), l_list, cl21T_nonu[i]**2/yerr2_nonu**2)))
	w2 = 1/yerr2**2
	
	if i == 0:
		print (yerr1)
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
		#xx2.append (l_mid)
		#yy2.append (interp_gT(l_mid))
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
	#cl21T_nonu[i] = np.sqrt(cl21T_nonu[i])

	plt.errorbar(xx1, yy1, yerr=l_bin_yerr2, fmt = None, ecolor='r')
	plt.errorbar(xx1, yy1, yerr=l_bin_yerr1, fmt = None, ecolor='b')
	plt.plot(l_list, cl21T_nu1[i], color = 'k')
	#plt.plot(l_list, cl21T_nonu[i], '--k')
	#plt.axis([2,10**7,10**-9, 10**-2])
	plt.axis([2,3000,10**-5, 10**-2])
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
	
	
	if i in [5,6,7,8]:
		plt.xlabel (r'$l$',size=15)
		if not i == 7:
			fig.set_xticks([10**1,10**2,10**3])
		else:
			fig.set_xticks([10**1,10**2,10**3])
	if i == 0:
		fig.set_yticks([10**-4,10**-3,10**-2])
	if i != 0:
		fig.set_yticks([10**-4,10**-3])
		#fig.set_yticks([10**-5,10**-4,10**-3])
	if i in [3]:
		plt.ylabel (r'$[l(l+1)C_l^{21,T}/2\pi]^{1/2}$ mK', size = 15)
	if i not in [5,6,7]:
		fig.axes.get_xaxis().set_visible(False)
	if i not in [0,3,6]:
		fig.axes.get_yaxis().set_visible(False)
	#plt.locator_params(axis='x', nticks=4)
	plt.tick_params(labelsize=11)
	
	plt.subplots_adjust(wspace=0, hspace=0)
#plt.suptitle (r'$C_l^{21,\mathrm{ISW}}$, $m_{\nu, i} = 0\mathrm{eV}$, $\mathrm{N}_{\mathrm{eff}}=3.2$',size = 15)
#plt.suptitle (r'$C_l^{21,T}$, $\sum m_{\nu} = 0.06\mathrm{eV}$',size = 12)
plt.savefig('err_cl21T.pdf',format='pdf')
#plt.savefig('err_cl21T_nonu_Neff.pdf',format='pdf')

plt.figure(4)
#plt.plot(z_m_list, s_to_n_nonu, label = r'$m_{\nu}=0$eV, $f_{\mathrm{sky}}=1$')
#plt.plot(z_m_list, [0.56720228967111386, 2.731606882172207, 9.0185069290780131, 20.745855126580814, 39.29962729884511, 65.655389715112321, 100.58699927848964, 145.1996402201006], label = r'$m_{\nu}=0$eV, $N_{\mathrm{eff}}=3.2$, $f_{\mathrm{sky}}=1$')
#plt.plot(z_m_list, [7.4273895654087108, 17.367165363662984, 34.710340384545077, 58.777780356495121, 90.535434783192756, 130.62483459979222, 179.59688299697544, 238.77418108021374], label = r'$m_{\nu}=0.1$eV, $f_{\mathrm{sky}}=1$')
plt.plot(z_m_list, s_to_n, 'b')#, label = r'$m_{\nu}=0.1$eV, $f_{\mathrm{sky}}=1$')
plt.plot(z_m_list, s_to_n_TE, 'r')#, label = r'$m_{\nu}=0.1$eV, $f_{\mathrm{sky}}=1$')
#plt.plot(z_m_list, s_to_n_nonu_TE, label = r'$m_{\nu}=0$eV, $f_{\mathrm{sky}}=0.4$')
#plt.plot(z_m_list, [0.22676776904733198, 1.0920752117510693, 3.6054661588421695, 8.2937662481677776, 15.710964784667064, 26.246971219732995, 40.211042881272249, 58.044883188166757], label = r'$m_{\nu}=0$eV, $N_{\mathrm{eff}}=3.2$, $f_{\mathrm{sky}}=0.4$')
#plt.plot(z_m_list, [2.9686241989554452, 6.9415468695067526, 13.873825137703216, 23.494123813879106, 36.188579529719085, 52.213612006184682, 71.789393024152361, 95.444546620564168], label = r'$m_{\nu}=0.1$eV, $f_{\mathrm{sky}}=0.4$')
#plt.plot(z_m_list, s_to_n_TE, label = r'$m_{\nu}=0.1$eV, $f_{\mathrm{sky}}=0.4$')
print (s_to_n)
print (s_to_n_TE)

#plt.yscale('log')
#plt.xscale('log')

#plt.plot(z_m_list, s_to_n_TE, label = 'nu TE')
#plt.plot(z_m_list, s_to_n_nonu_TE, label = 'nonu TE')

plt.grid ()
plt.ylabel (r'$(S/N)^2$ for each window function',size=15)
plt.xlabel (r'$z_m$',size=15)
#plt.axis([0,200,0,450])
#plt.axis([0,30,0,18])
#plt.legend(bbox_to_anchor=(0.5,1),prop={'size':12})
#plt.title (r'Signal to Noise of $C_l^{21,T}$')
#plt.axis([0,100,0,180])
plt.savefig ('SN_cl21T.pdf', format = 'pdf')

plt.show ()


