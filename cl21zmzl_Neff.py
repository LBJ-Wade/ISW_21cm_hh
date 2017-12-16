import numpy as np
import matplotlib.pyplot as plt

#infile_dev = "result/cl21T_deriv_Neff.txt"
infile_dev = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_Yp_data_21/cl21T_deriv_Neff.txt"
infile_dev = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/result_Yp_BBN/cl21T_deriv_Neff.txt"

#infile_0 = "result/cl21T_0.txt"
#infile_TT = "result/cl_0.dat"
infile_0 = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_Yp_data_21/cl21T_0.txt"
infile_TT = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_Yp_data_21/cl_0.dat"
infile_0 = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/result_Yp_BBN/cl21T_0.txt"
infile_TT = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/data/cl_0_BBN.dat"
l = np.loadtxt (infile_TT)[0:,0]
aa = 10**-12 * 2*np.pi / (l*(l+1))
clTT = np.loadtxt (infile_TT)[0:,1]*aa
clTT_N = (2*10**-6*0.000290888)**2 *np.e**(l*(l+1)*(0.000290888**2)/(8*np.log(2)))
clTT += clTT_N
sigma_list = []
f_sky = 0.4
for i in range(8):
	dev_cl = np.loadtxt (infile_dev) [0:,i+1] *2.7255
	cl21T = np.loadtxt (infile_0) [0:,i+1] *2.7255
	#infile2 = "result/cl21zmzl_{0}{1}.txt".format(i,i)
	#infile2 = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_Yp_data_21/cl21zmzl_{0}{1}.txt".format(i,i)
	infile2 = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/ISW_21cm/result_Yp_BBN/cl21zmzl_{0}{1}.txt".format(i,i)	
	cl21 = np.loadtxt (infile2) [0:,1]
	cov = (cl21*clTT + cl21T**2)/(2*l+1)/f_sky
	result = dev_cl**2 / cov
	result = sum(result[28:2999])
	sigma_list.append (np.sqrt(1/result))
print (sigma_list)

z_m_list = [30, 50,75,100,125,150,175,200,225]
z_m_list = np.array(z_m_list)
sigma_list_ori = [3.8957e-01, 2.9930e-01, 2.5857e-01, 2.4081e-01, 2.3287e-01, 2.2995e-01, 2.2951e-01, 2.3008e-01, 1.1900e-01]
sigma_list_modi = [3.8969e-01, 2.9953e-01, 2.5922e-01, 2.4199e-01, 2.3446e-01, 2.3156e-01, 2.3143e-01, 2.3834e-01,1.1966e-01]
sigma_list_2 = [3.8959e-01, 2.9930e-01, 2.5875e-01, 2.4140e-01, 2.3414e-01, 2.3224e-01, 2.3322e-01, 2.3562e-01, 1.1992e-01]
sigma_list_1 = [3.8961e-01, 3.0021e-01, 2.6240e-01, 2.5127e-01, 2.5767e-01,2.8253e-01, 2.6599e-01, 1.1517e-01, 6.0989e-02]
plt.figure(1)
fiducial = np.ones(9)*0.06
my_xticks = ['30','50','75','100','125','150','175','200','All']
plt.xticks(z_m_list, my_xticks)
plt.plot (z_m_list, fiducial)

plt.errorbar (z_m_list-3, fiducial, yerr=sigma_list_ori, fmt = None, ecolor='red', label = r'$k$=7e-5')
plt.errorbar (z_m_list-1, fiducial, yerr=sigma_list_2, fmt = None, ecolor='y',label = r'$k$=1e-2' )
plt.errorbar (z_m_list+1, fiducial, yerr=sigma_list_1, fmt = None, ecolor='blue', label = r'$k$=1e-1')
plt.errorbar (z_m_list+3, fiducial, yerr=sigma_list_modi, fmt = None, ecolor='purple', label = r'$k$-dep')
plt.legend (prop={'size':12})

#plt.yscale ('log')
plt.xlabel (r'$z_m$',size=20)
plt.ylabel (r'$\sum m_\nu$',size=20)
plt.axis ([20,235,-0.5,0.6])
plt.savefig ('mnu_z.pdf')



sigma_list_ori = [1.2756e+01, 6.5728e+00, 3.3955e+00, 2.1072e+00, 1.4927e+00, 1.1491e+00, 9.3280e-01, 7.8479e-01, 6.711e-01]	# ori
sigma_list_modi = [1.2746e+01, 6.5781e+00, 3.3956e+00, 2.1027e+00,1.4868e+00,1.1421e+00,9.1892e-01, 7.6967e-01, 6.5392e-01]	#modi
sigma_list_2 = [1.2748e+01, 6.5763e+00, 3.3757e+00, 2.0714e+00, 1.4494e+00, 1.1032e+00, 8.8733e-01,7.4179e-01,6.3600e-01]
sigma_list_1 = [1.2769e+01, 6.5761e+00, 3.4553e+00, 2.2323e+00, 1.6667e+00, 1.2220e+00, 6.1641e-01,2.7773e-01, 2.0598e-01]
plt.figure(2)
fiducial = np.ones(9)*3.046
my_xticks = ['30','50','75','100','125','150','175','200','All']
plt.xticks(z_m_list, my_xticks)

plt.plot (z_m_list, fiducial)
plt.errorbar (z_m_list-3, fiducial, yerr=sigma_list_ori, fmt = None, ecolor='red', label = r'$k$=7e-5')
plt.errorbar (z_m_list-1, fiducial, yerr=sigma_list_2, fmt = None, ecolor='y',label = r'$k$=1e-2' )
plt.errorbar (z_m_list+1, fiducial, yerr=sigma_list_1, fmt = None, ecolor='blue', label = r'$k$=1e-1')
plt.errorbar (z_m_list+3, fiducial, yerr=sigma_list_modi, fmt = None, ecolor='purple', label = r'$k$-dep')
plt.legend (prop={'size':12})

#plt.yscale ('log')
plt.xlabel (r'$z_m$',size=20)
plt.ylabel (r'$N_{\mathrm{eff}}$',size=20)
#plt.axis ([20,210,1,5])
plt.axis ([20,235,1,5])
plt.savefig ('Neff_z.pdf')

"""
fiducial = np.ones(8)*0.06
sigma_mnu = [0.329, 0.269, 0.239, 0.225, 0.220, 0.219, 0.218, 0.214]
plt.figure(3)
plt.plot (z_m_list, fiducial)
plt.errorbar (z_m_list, fiducial, yerr=sigma_mnu, fmt = None, ecolor='r')
plt.xlabel (r'$z_m$')
plt.ylabel (r'$\sum m_\nu$')
plt.axis ([20,210,-0.3,0.4])
plt.savefig ('mnu_z.pdf')
"""


plt.show ()
