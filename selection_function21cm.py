import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
z_m = [30, 50,75,100,125,150,175,200]
w = [1.27689,1, 0.8215, 0.7081, 0.6299, 0.5719, 0.5265, 0.4897]
# width = [14.957, 27.580, 42.409, 59.295,  78.075, 98.588, 120.732]
w = [3.31685, 2.47355 ,1.93014, 1.60434, 1.38246, 1.21989, 1.09467, 1]
def run_sel(w,z_m):
	c = 299792458
	wavelength = 0.21106
	#z_m = 75
	f0 = 1420.4*10**6
	freq = f0/(1+z_m)
	width = w*10**6#/(1+z_m)
	freq1 = freq-4*width
	#print (freq1 >0)
	freq2 = freq+4*width
	f = np.linspace(freq1,freq2,5000)
	l = c/f
	z = l/wavelength-1
	z = f0/f-1
	wl1 = c/freq1
	wl2 = c/freq2
	z1 = wl1/wavelength-1
	z2 = wl2/wavelength-1
	sigma_wl = wl2-wl1
	
	zz = np.linspace(400,10,5000)

	sel_func = 1/np.sqrt(2*np.pi*width**2) * np.exp (-(f-freq)**2/(2*width**2))
	sel_func = f0/(1+zz)**2 * 1/np.sqrt(2*np.pi*width**2) * np.exp (-(f0/(1+zz)-f0/(1+z_m))**2/(2*width**2))
	
	#plt.plot(z[::-1], sel_func[::-1], label = 'f')
	#plt.axis([0,200,0,max(sel_func)])
	#plt.xlabel (r'$z$')
	#plt.title ('Window Function')
	#print (simps(sel_func[::-1], z[::-1]))
	#print (max(z), min(z), max(z)-min(z))	
	#plt.show()
	return zz[::-1], sel_func[::-1], max(z)-min(z)

#for i in range(len(w)):
#	run_sel(w[i],z_m[i])
#plt.grid()
#plt.axis([0,300,0,0.18])
#plt.savefig ('window_21.pdf',format='pdf')
#plt.show()

#run_sel(1.27689,30)
#run_sel(1.27689,50)
#run_sel(0.5,50)
#plt.show()
"""
plt.figure(1)
w = [0.01, 0.1, 1, 6]
for i in range(4):
	z, sel = run_sel(w[i])
	fig = plt.subplot (2,2,i+1)
	plt.plot(z, sel, label = '{}MHz'.format(w[i]))
	plt.xlabel (r'$z$')
	if i == 0:
		fig.set_xticks([49.9,50,50.1])
	if i == 1:
		fig.set_xticks([49.5,50,50.5])
	if i == 3:
		plt.axis([0,150,0,0.04])
		fig.set_xticks([0,50,100,150])
	plt.legend(bbox_to_anchor=(1,0.25),prop={'size':15})
plt.suptitle('Gaussian Window Function', size = 15)
plt.savefig ('sel21.pdf',format='pdf')
plt.show()
"""
