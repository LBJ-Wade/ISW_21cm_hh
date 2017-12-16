from cl_21 import *
from path import *
import time


Cl = set_cl_21('0',Yp_BBN = True)
l_list = np.arange(2,5001)
cl21T_list = []
#z_m_list = [30, 50,75,100,125,150,175,200]
#w_list = [3.1891, 2.4779, 2.0211, 1.7456, 1.5560, 1.4151, 1.3050, 1.2156]

z_m_list = [30]
w_list = [3.1891]

for i in range(len(z_m_list)):

	# Cross-correlation
	cl_zm = Cl.cl21T (z_m_list[i], w_list[i])
	
	plt.loglog (Cl.l_list, cl_zm)
	plt.show()
	# Auto-correlation
	#cl_zm = Cl.cl21 ([z_m_list[i],z_m_list[i]], [w_list[i],w_list[i]])
	cl21T_list.append (cl_zm)
