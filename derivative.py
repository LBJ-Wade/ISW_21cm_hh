import numpy as np

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

