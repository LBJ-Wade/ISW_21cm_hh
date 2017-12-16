import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
from scipy.interpolate import interp1d

def derivative (x, f):
	if len(x) != len(f):
		error = 'Dimension error'
		print (error)
		return None
	else:
		dif_f = []
		x_mid = []
		for i in range(len(x)-1):
			dx = x[i+1] - x[i]
			if dx == 0:
				pass
			else:
				dif = (f[i+1] - f[i]) / (x[i+1] - x[i])
				dif_f.append(dif)

				x_mid.append ( (x[i+1] + x[i])/2 )

		dif_f = np.interp(x, x_mid, dif_f)
		
		return dif_f

def integrate (x, f):
	integ_f = 0
	for i in range(len(x)-1):
		integ_f += (f[i+1] + f[i])/2 * (x[i+1] - x[i])
	return integ_f

def integrate2 (x, f):
	if len(x) != len(f):
		error = 'Dimension error'
		print (error)
		return None
	else:
		new_x = []
		for i in range(len(x)-1):
			new_x.append((x[i+1]+x[i])/2)
		new_f = np.interp (new_x, x, f)

		integ_f = 0
		for i in range(len(x)-1):
			integ_f += (f[i+1] + 4*new_f[i] + f[i]) * (x[i+1] - x[i])/6
		return integ_f
