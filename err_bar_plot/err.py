#from D import *
#from isw_transfer import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator
aa = 1/(10**6*2.7255**2)
lmax = 2999
#lmax = 4999

infile = "m0_cl.dat"
cl_TT = np.loadtxt(infile)[0:lmax,1]/10**6
cl_TE = np.loadtxt(infile)[0:lmax,3]/10**6
cl_EE = np.loadtxt(infile)[0:lmax,2]/10**6


cl21T_nonu = []
cl21T_nonu_Neff = []
cl21T_nu1 = []
cl21T_nu2 = []
cl21T_N = []
cl21_nonu = []
cl21_nonu_Neff = []
cl21_nu1 = []
cl21_nu2 = []
cl21_N = []
n = 8


infile_nonu1 = "cl21T_m0.txt" 
infile_nonu2 = "cl21_m0.txt"
infile_nu1 = "cl21T_Neff1.txt" 
infile_nu2 = "cl21_Neff1.txt"
infile_nu1_2 = "cl21T_Neff2.txt"
infile_nu2_2 = "cl21T_Neff2.txt"
infile_N1 = "cl21T_Neff2.txt"
infile_N2 = "cl21_Neff2.txt"

z_m_list = [30,50,75,100,125,150,175,200]
cl21T = []
l_list = np.loadtxt (infile_nonu1) [0:lmax,0]
aaa = 10**6* l_list*(l_list+1)/(2*np.pi) * 2.7255
for i in range(n):
	#cl21T_nonu.append ()
	cl21T_nonu.append (np.loadtxt (infile_nonu1)[0:lmax,i+1]*aaa)
	cl21_nonu.append (np.loadtxt (infile_nonu2)[0:lmax,i+1]*aaa)
	cl21T_nu1.append (np.loadtxt (infile_nu1)[0:lmax,i+1]*aaa)
	cl21_nu1.append (np.loadtxt (infile_nu2)[0:lmax,i+1]*aaa)
	cl21T_nu2.append (np.loadtxt (infile_nu1_2)[0:lmax,i+1]*aaa)
	cl21_nu2.append (np.loadtxt (infile_nu2_2)[0:lmax,i+1]*aaa)
	cl21T_N.append (np.loadtxt (infile_N1)[0:lmax,i+1]*aaa)
	cl21_N.append (np.loadtxt (infile_N2)[0:lmax,i+1]*aaa)

	#cl21_nonu.append (nonu21)
ll = l_list
cl_TT_N = (1*0.000290888*10**-3)**2 *np.e**(ll*(ll+1)*(10**-3*1*0.000290888)**2/(8*np.log(2)))

cl_TT_N = cl_TT_N * ll*(ll+1)/(2*np.pi)

cl_TT += cl_TT_N
#cl_TT_nonu += cl_TT_N
l_bin_index = [2, 4, 9, 20, 45, 100, 218, 478, 1045, 2286, 5000]
l_bin_index = [2, 4, 9, 18, 37, 77, 161, 334, 695, 1444, 3000]

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
	f_sky = 0.4
	print (len(cl_TT), len(cl21_nu1[i]), len(l_list))
	
	#yerr1 = np.sqrt ((cl_TT*(cl21_nu1[i]) + cl21T_nu1[i]**2)/(2*l_list+1)/f_sky)
	#w1 = 1/yerr1**2 	
	#s_to_n.append (sum(np.interp (np.arange(2,3001), l_list, cl21T_nu1[i]**2/yerr1**2)))
	
	yerr1 = np.sqrt ((cl_TT*(cl21_nonu[i]) + cl21T_nonu[i]**2)/(2*l_list+1)/f_sky)
	w1 = 1/yerr1**2 	
	s_to_n.append (sum(np.interp (np.arange(2,3001), l_list, cl21T_nonu[i]**2/yerr1**2)))
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
	
	#interp_gT = interp1d (l_list, cl21T_nu1[i])
	interp_gT = interp1d (l_list, cl21T_nonu[i])
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

	#yy1 = np.sqrt(yy1)
	#l_bin_yerr1 = np.sqrt(l_bin_yerr1)
	#l_bin_yerr2 = np.sqrt(l_bin_yerr2)
	#cl21T_nu1[i] = np.sqrt(cl21T_nu1[i])
	#cl21T_nu2[i] = np.sqrt(cl21T_nu2[i])
	#cl21T_nonu[i] = np.sqrt(cl21T_nonu[i])
	#cl21T_N[i] = np.sqrt(cl21T_N[i])

	#plt.errorbar(xx1, yy1, yerr=l_bin_yerr2, fmt = None, ecolor='k')
	plt.errorbar(xx1, yy1, yerr=l_bin_yerr1, fmt = None, ecolor='steelblue', elinewidth = 1.5)
	plt.plot(l_list, cl21T_nu2[i], color = 'blue')
	plt.plot(l_list, cl21T_nu1[i], color = 'y')
	plt.plot(l_list, cl21T_nonu[i], 'red')
	#plt.axis([2,10**7,10**-9, 10**-2])
	
	plt.axis([2,lmax+1,8*10**-8, 6*10**-5])
	
	plt.text(3,  2.5*10**-5, r'z$_m=${}'.format(z_m_list[i]))	
	
			
	if i in [5,6,7,8]:
		plt.xlabel (r'$l$',size=15)
		if not i == 7:
			fig.set_xticks([10**1,10**2,10**3])
		else:
			fig.set_xticks([10**1,10**2,10**3])
	if i == 0:
		fig.set_yticks([10**-7,10**-6,10**-5])
	if i != 0:
		fig.set_yticks([10**-7,10**-6,10**-5])
		#fig.set_yticks([10**-5,10**-4,10**-3])
	if i in [3]:
		plt.ylabel (r'$l(l+1)C_l^{21,T}/2\pi$ $\mathrm{mK}^2$', size = 15)
	if i not in [5,6,7]:
		fig.axes.get_xaxis().set_visible(False)
	if i not in [0,3,6]:
		fig.axes.get_yaxis().set_visible(False)
	
	#plt.locator_params(axis='x', nticks=4)
	plt.tick_params(labelsize=11)
	
fig = plt.subplot(3,3,9)
fig.axes.get_yaxis().set_visible(False)
plt.xlabel (r'$l$',size=15)
label_list = [r'$\sum m_\nu = 0\mathrm{eV}$',r'$\sum m_\nu = 0.1\mathrm{eV}$',r'$\sum m_\nu = 0.2\mathrm{eV}$']
label_listN = [r'$\Delta N_\mathrm{eff} = 0$',r'$\Delta N_\mathrm{eff} = 1$',r'$\Delta N_\mathrm{eff} = 2$']
#plt.plot (1,1,label = label_list[2], color = 'blue')
#plt.plot (1,1,label = label_list[1], color = 'y')
#plt.plot (1,1,label = label_list[0], color = 'red')
#plt.legend(prop={'size':12})
plt.plot (1,1,label = label_listN[2], color = 'blue')
plt.plot (1,1,label = label_listN[1], color = 'y')
plt.plot (1,1,label = label_listN[0], color = 'red')
plt.legend(prop={'size':15})
plt.axis([2,3000,10**-4, 10**-2])
fig.set_xticks([10**1,10**2,10**3])
plt.xscale ('log')



plt.subplots_adjust(wspace=0, hspace=0)
#plt.suptitle (r'$C_l^{21,\mathrm{ISW}}$, $m_{\nu, i} = 0\mathrm{eV}$, $\mathrm{N}_{\mathrm{eff}}=3.2$',size = 15)
#plt.suptitle (r'$C_l^{21,T}$, $\sum m_{\nu} = 0.06\mathrm{eV}$',size = 12)
plt.savefig('err_cl21T.pdf',format='pdf')
#plt.savefig('err_cl21T_nonu_Neff.pdf',format='pdf')

plt.figure(4)
#plt.plot(z_m_list, np.array(s_to_n)/f_sky, 'b')#, label = r'$m_{\nu}=0.1$eV, $f_{\mathrm{sky}}=1$')
plt.plot(z_m_list, np.array([0.50562674098962734, 1.9380106813842937, 5.694145937689818, 12.295061582412195, 22.225559761127947, 35.67384244926442, 52.624493040563969, 72.918805398155257])/f_sky, 'blue', label = label_list[2])
plt.plot(z_m_list, np.array([0.39923288211650765, 1.6776739520147197, 5.1853564442414291, 11.475073370957279, 21.037884706292953, 34.073560315243434, 50.573380502489968, 70.384710090732867])/f_sky, 'y', label = label_list[1])
#plt.plot(z_m_list, np.array([0.30309031064112607, 1.4481490849045111, 4.744970803932488, 10.781740905331233, 20.061935480527218, 32.797851771685551, 48.990334447755309, 68.495072079537394])/f_sky, 'red', label = label_list[0])

plt.plot(z_m_list, np.array([0.37344715293089364, 1.8003485962948005, 5.9171941742411089, 13.415790506175139, 24.85835737888673, 40.437296363000549, 60.084781621625623, 83.56578022744435])/f_sky, '--',color = 'blue', label = label_listN[2])
plt.plot(z_m_list, np.array([0.33955944601160387, 1.6299524218641799, 5.3495162936873131, 12.142710254628378, 22.546856475037284, 36.767458814432992, 54.772648226609938, 76.372633018494753])/f_sky, '--',color = 'y', label = label_listN[1])
plt.plot(z_m_list, np.array([0.30309031064112607, 1.4481490849045111, 4.744970803932488, 10.781740905331233, 20.061935480527218, 32.797851771685551, 48.990334447755309, 68.495072079537394])/f_sky, 'red', label = label_list[0] + ', '+label_listN[0])
plt.legend (bbox_to_anchor=(0.4,1),prop={'size':12})
print (s_to_n)
print (s_to_n_TE)

#plt.yscale('log')
#plt.xscale('log')

#plt.plot(z_m_list, s_to_n_TE, label = 'nu TE')
#plt.plot(z_m_list, s_to_n_nonu_TE, label = 'nonu TE')

plt.grid ()
plt.ylabel (r'$(S/N)^2 f_\mathrm{sky}^{-1}$ for each window function',size=15)
plt.xlabel (r'$z_m$',size=15)
#plt.axis([0,200,0,450])
#plt.axis([0,30,0,18])
#plt.legend(bbox_to_anchor=(0.5,1),prop={'size':12})
#plt.title (r'Signal to Noise of $C_l^{21,T}$')
plt.axis([20,200,0,200])
plt.savefig ('SN_cl21T.pdf', format = 'pdf')
print (f_sky)
plt.show ()


