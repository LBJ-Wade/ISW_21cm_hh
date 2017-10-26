from path import *
from cl_21 import *
from run import *

params_list = np.loadtxt (params_input)[0:,]
tag = "test0"
run_21cm (params_list, tag)
Cl1 = set_cl_21 (tag)
l_list = Cl1.l_list
cl1_zm = Cl1.cl21T (30, 3.31685)
data = np.column_stack(( l_list, cl1_zm))
np.savetxt('ns_{0}.txt'.format(tag), data, fmt = '%1.6e')
plt.loglog (l_list, np.sqrt(l_list*(l_list+1)/(2*np.pi)*cl1_zm))
plt.show()
