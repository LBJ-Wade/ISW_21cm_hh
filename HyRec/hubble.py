import numpy as np
import data_no_nu as data
import selection_function21cm as sel

infile = data.transfer_syn

z = np.loadtxt(infile)[0:,0]
print ('z data', len(z))
k = np.loadtxt(infile)[0:,1]
print ('k data', len(k))
zlist2 = np.array(sorted(set(z)))
klist2 = np.array(sorted(set(k)))
number_of_z2 = len(zlist2)
number_of_k2 = len(klist2)
hubble_class = np.loadtxt(infile)[0:number_of_z2,10][::-1]
print ('hubble data', len(hubble_class))

z_m = [75, 100, 125, 150, 175, 200]
z_m = [30,75]
#z_m = [50,100,125,150,175,200]
_, _, width_30 = sel.run_sel (1, 200)
print (width_30)
new_width_list = []
for i in range(len(z_m)):
	#new_width = 14.9579339279 / np.interp (50, zlist2, hubble_class) * np.interp(z_m[i], zlist2, hubble_class)
	#new_width = 0.0541256396731 / np.interp (30, zlist2, hubble_class) * np.interp(z_m[i], zlist2, hubble_class)
	new_width = width_30 / np.interp (200, zlist2, hubble_class) * np.interp(z_m[i], zlist2, hubble_class)

	error = 100
	w = 1
	while error > 10**-2*new_width:
		_, _, width = sel.run_sel(w, z_m[i])
		error = abs(width-new_width)
		w += 0.00001
		print (z_m[i],w, new_width, error)
	new_width_list.append ([w, width, new_width, error])
print (new_width_list)
