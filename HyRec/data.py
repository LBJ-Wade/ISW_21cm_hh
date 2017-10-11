params = 'nanoom_cb'
#params = 'nanoom2_cb'
#params = 'nanoom_Neff'

if params == 'nanoom_cb':
	transfer_new = "params/delta_ncdm3_nanoom1_new.dat"
	transfer_syn = "params/delta_ncdm3_nanoom1.dat"
	transfer_21 = "transfer_21_nu1.txt"
	power_spectrum = "params/params_nu3_nanoom1_cb_00_pk.dat"
	cl = "params/params_nu3_nanoom1_cb_00_cl.dat"
elif params == 'nanoom2_cb':
	transfer_new = "params/delta_ncdm3_nanoom2_new.dat"
	transfer_syn = "params/delta_ncdm3_nanoom2.dat"
	transfer_21 = "transfer_21_nu2.txt"
	power_spectrum = "params/params_nu3_nanoom2_cb_00_pk.dat"
	cl = "params/params_nu3_nanoom2_cb_00_cl.dat"
elif params == 'nanoom_Neff':
	transfer_new = "params/delta_no_nu_nanoom_Neff_new.dat"
	transfer_syn = "params/delta_no_nu_nanoom_Neff.dat"
	transfer_21 = "transfer_21_nonu_Neff.txt"
	power_spectrum = "params/params_no_nu_nanoom_Neff_00_pk.dat"
	cl = "params/params_no_nu_nanoom_Neff_00_cl.dat"

