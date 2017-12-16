import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
z_m_list = [30,50,75,100,125,150,175,200]
infile_list = ['cl21T_m0.txt', 'cl21T_m1.txt', 'cl21T_m2.txt']
color = ['red','y','blue']
for j in range(len(infile_list)):
	infile = infile_list[j]
	l_list = np.loadtxt (infile)[0:,0]
	aaa = 10**6* l_list*(l_list+1)/(2*np.pi) * 2.7255
	for i in range(8):
		fig = plt.subplot(3,3,i+1)
		if j == 0:
			plt.text(30,  3*10**-4, r'z$_m=${}'.format(z_m_list[i]))
		plt.plot(l_list, np.sqrt(np.loadtxt (infile)[0:,i+1]*aaa), color = color[j])
		plt.axis([2,3000,10**-4, 10**-2])
		plt.xscale ('log')
		plt.yscale ('log')
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
		if i in [3]:
			plt.ylabel (r'$[l(l+1)C_l^{21,T}/2\pi]^{1/2}$ mK', size = 15)
		if i not in [5,6,7]:
			fig.axes.get_xaxis().set_visible(False)
		if i not in [0,3,6]:
			fig.axes.get_yaxis().set_visible(False)
		plt.tick_params(labelsize=11)
fig = plt.subplot(3,3,9)
fig.axes.get_yaxis().set_visible(False)
plt.xlabel (r'$l$',size=15)
label_list = [r'$\sum m_\nu = 0\mathrm{eV}$',r'$\sum m_\nu = 0.1\mathrm{eV}$',r'$\sum m_\nu = 0.2\mathrm{eV}$']
for i in range(len(label_list)):
	plt.plot (1,1,label = label_list[-1-i], color = color[-1-i])
plt.legend(prop={'size':12})
plt.axis([2,3000,10**-4, 10**-2])
fig.set_xticks([10**1,10**2,10**3])
plt.xscale ('log')
#bbox_to_anchor=(1,1),prop={'size':15}
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig ('cl21T.pdf')
plt.show()

