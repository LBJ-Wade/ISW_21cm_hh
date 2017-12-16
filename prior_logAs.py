import numpy as np
from numpy.linalg import inv

infile_planck = "/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/prior_planck_logAs.txt"
#infile_cmb = '/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_result/fisher_matrix_woNeff_07.txt'
#infile_cl21T = '/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_result/fisher_matrix_cl21T_l3000_fsky0.7.txt'
#infile_cmb = '/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_result/fisher_matrix_woNeff_07.txt'
#infile_cmb = 'fisher_matrix_new_nl_woNeff.txt'
#infile_cl21T = '/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_result/fisher_matrix_cl21T_l3000_fsky0.7.txt'

# Without Yp
#infile_cmb = 'fisher_matrix_new_nl_noise_woNeff.txt'
#infile_cl21T = '/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_result/fisher_matrix_cl21T_l3000_fsky0.7.txt'
#infile_zm = '/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_result/fisher_matrix_cl21T_200.txt'

# With Yp
infile_cmb = 'fisher_matrix_new_nl_Yp_logAs_noise_woNeff.txt'
infile_cl21T = '/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_Yp_result/fisher_matrix_cl21T_Yp_l3000_fsky0.7.txt'
infile_zm = '/Users/NanoomLee/Documents/Master@StonyBrook/ISW_21cm_Code/logAs_Yp_result/fisher_matrix_cl21T_Yp_150.txt'

prior_planck = np.zeros([8,8])
prior_cmb = np.zeros([8,8])
fisher_cl21T = []
prior_zm = []
sub_prior_planck = []
sub_prior_cmb = []
for i in range(6):
	a = np.loadtxt (infile_planck)[0:,i]
	prior_planck[i][:6] = a.copy() 
#	sub_prior_planck.append (a)

for i in [0,1,2,3,4,5]:
	a = []
	for j in [0,1,2,3,4,5]:
		a.append (np.loadtxt(infile_planck)[0:,i][j])
	sub_prior_planck.append (np.array(a))


for i in range(7):
	a = np.loadtxt (infile_cmb)[0:,i]
	prior_cmb[i][:7] = a.copy()
	sub_prior_cmb.append (a)

for i in range(8):
	a = np.loadtxt (infile_cl21T)[0:,i]
	fisher_cl21T.append (a)
	b = np.loadtxt (infile_zm)[0:,i]
	prior_zm.append (b)

sub_prior_planck = np.array (sub_prior_planck)
sub_prior_cmb = np.array (sub_prior_cmb)
fisher_cl21T = np.array (fisher_cl21T)
prior_zm = np.array (prior_zm)

fisher1 = []
for i in range(7):
	fisher1.append (prior_planck[i][:-1] + prior_cmb[i][:-1])





inv_prior_planck = inv (sub_prior_planck)
inv_prior_cmb = inv (sub_prior_cmb)
inv_fisher_cl21T = inv (fisher_cl21T)
inv_prior_zm = inv (prior_zm)

"""
# remove tau
print (fisher_cl21T)
fisher_cl21T[3] *= 0
fisher_cl21T = fisher_cl21T.transpose()
fisher_cl21T[3] *= 0
fisher_cl21T = fisher_cl21T.transpose()
print (fisher_cl21T)
"""

fisher2 = prior_planck + prior_cmb + fisher_cl21T
inv_fisher1 = inv (fisher1)
inv_fisher2 = inv (fisher2)

sigma1 = []
sigma2 = []
for i in range(7):
	sigma1.append (np.sqrt(inv_fisher1[i,i]))
	sigma2.append (np.sqrt(inv_fisher2[i,i]))
sigma1 = np.array (sigma1)
sigma2 = np.array (sigma2)
improve = (sigma1-sigma2)/sigma1 * 100

for i in [0,1,5,3,4,6,2,7]:#range(len(inv_fisher1)):
	if i < 4:
		print ('%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e' % (np.sqrt(inv_prior_planck[i,i]), np.sqrt(inv_prior_cmb[i,i]), np.sqrt(inv_fisher_cl21T[i,i]), np.sqrt(inv_fisher1[i,i]), np.sqrt(inv_fisher2[i,i]), improve[i]))
	elif i == 4:
		print ('%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e' % (np.sqrt(inv_prior_planck[i,i]), np.sqrt(inv_prior_cmb[i,i]), np.sqrt(inv_fisher_cl21T[i,i]), np.sqrt(inv_fisher1[i,i]), np.sqrt(inv_fisher2[i,i]), improve[i]))
	elif i == 5:
		print ('%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e' % (np.sqrt(inv_prior_planck[i,i]), np.sqrt(inv_prior_cmb[i,i]), np.sqrt(inv_fisher_cl21T[i,i]), np.sqrt(inv_fisher1[i,i]), np.sqrt(inv_fisher2[i,i]), improve[i]))
	elif i == 6:
		print ('%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e' % (0, np.sqrt(inv_prior_cmb[i,i]), np.sqrt(inv_fisher_cl21T[i,i]), np.sqrt(inv_fisher1[i,i]), np.sqrt(inv_fisher2[i,i]), improve[i]))
	else:
		print ('%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e' % (0, 0, np.sqrt(inv_fisher_cl21T[i,i]), 0, np.sqrt(inv_fisher2[i,i]),0))

fisher3 = prior_planck + prior_cmb + prior_zm
inv_fisher3 = inv (fisher3)
for i in [0,1,5,3,4,6,2,7]:#range(len(inv_fisher1)):
	print ('%1.4e %1.4e' % (np.sqrt(inv_prior_zm[i,i]), np.sqrt(inv_fisher3[i,i])))
