import numpy as np
import matplotlib.pyplot as plt
syn_pk = 'params_no_nu_nanoom_00_pk.dat'
new_pk = 'params_no_nu_nanoom_00_new_pk.dat'

syn_k = np.loadtxt(syn_pk)[0:,0]
syn_Pm = np.loadtxt(syn_pk)[0:,1]
new_k = np.loadtxt(new_pk)[0:,0]
new_Pm = np.loadtxt(new_pk)[0:,1]

plt.loglog(syn_k, syn_Pm, label = 'syn')
plt.loglog(new_k, new_Pm, label = 'new')
plt.legend()
plt.show()

