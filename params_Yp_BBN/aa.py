import numpy as np
infile = 'params_nonu.dat'

a = np.loadtxt (infile)[0:,]
print type(a[0])
