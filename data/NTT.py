import numpy as np
import matplotlib.pyplot as plt

infile = 'cl_0.dat'
l = np.loadtxt (infile)[0:,0]
clTT = np.loadtxt (infile)[0:,1] 
clTT_N = (0.1*0.000290888)**2 *np.e**(l*(l+1)*(0.000290888**2)/(8*np.log(2))) * l*(l+1)/(2*np.pi)

plt.loglog (l, clTT)
plt.loglog (l, clTT_N)
plt.show()
