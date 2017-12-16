import numpy as np

NTR = 100
NTM = 40
TR_MIN = 0.004
TR_MAX = 0.4
TM_TR_MIN = 0.1
TM_TR_MAX = 1.0

def maketab (xmin, xmax, Nx):
	h = (xmax - xmin)/(Nx-1.0)
	xtab = np.ones(Nx)*xmin + np.arange(Nx)*h
	return xtab	

def interpolate_rates (TR, TM_TR, logAlpha_tab, logR2p2s_tab, DTM_TR, DlogTR):
	global TM_TR_MIN, TR_MIN, NTM, NTR
	coeff1 = np.zeros(4)
	coeff2 = np.zeros(4)
	temp = np.zeros(4)
	logTR = np.log(TR)
	Alpha = np.zeros(2)
	Beta = np.zeros(2)
	R2p2s = None
	iTM = int(np.floor ( (TM_TR - TM_TR_MIN)/DTM_TR ))
	if iTM < 1:
		iTM = 1
	if iTM > NTM-3:
		iTM = NTM-3
	frac1 = (TM_TR - TM_TR_MIN)/DTM_TR - iTM
	coeff1[0] = frac1*(frac1-1.)*(2.-frac1)/6.
	coeff1[1] = (1.+frac1)*(1.-frac1)*(2.-frac1)/2.
	coeff1[2] = (1.+frac1)*frac1*(2.-frac1)/2.
	coeff1[3] = (1.+frac1)*frac1*(frac1-1.)/6.

	iTR = int(np.floor ( (logTR - np.log(TR_MIN))/DlogTR ))
	if iTR < 1:
		iTR = 1
	if iTR > NTR-3:
		iTR = NTR-3
	frac2 = (logTR - np.log(TR_MIN))/DlogTR - iTR
	coeff2[0] = frac2*(frac2-1.)*(2.-frac2)/6.
	coeff2[1] = (1.+frac2)*(1.-frac2)*(2.-frac2)/2.
	coeff2[2] = (1.+frac2)*frac2*(2.-frac2)/2.
	coeff2[3] = (1.+frac2)*frac2*(frac2-1.)/6.

	for l in range(2):
		for k in range(4):
			temp[k] = logAlpha_tab[l][iTM-1+k][iTR-1]*coeff2[0]\
					+ logAlpha_tab[l][iTM-1+k][iTR]*coeff2[1]\
					+ logAlpha_tab[l][iTM-1+k][iTR+1]*coeff2[2]\
					+ logAlpha_tab[l][iTM-1+k][iTR+2]*coeff2[3]
		Alpha[l] = np.exp(temp[0]*coeff1[0]+temp[1]*coeff1[1]\
					+temp[2]*coeff1[2]+temp[3]*coeff1[3])
		Beta[l] = np.exp(logAlpha_tab[l][NTM-1][iTR-1]*coeff2[0]\
				+ logAlpha_tab[l][NTM-1][iTR]*coeff2[1]\
				+ logAlpha_tab[l][NTM-1][iTR+1]*coeff2[2]\
				+ logAlpha_tab[l][NTM-1][iTR+2]*coeff2[3])

	factor = 3.016103031869581e21 *  TR*np.sqrt(TR) * np.exp(-3.399571517984581/TR)
	Beta[0] *= factor
	Beta[1] *= factor/3.

	R2p2s = np.exp(logR2p2s_tab[iTR-1]*coeff2[0]\
			+ logR2p2s_tab[iTR]*coeff2[1]\
			+ logR2p2s_tab[iTR+1]*coeff2[2]\
			+ logR2p2s_tab[iTR+2]*coeff2[3])
	return Alpha, Beta, R2p2s	
