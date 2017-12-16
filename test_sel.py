import numpy as np
from selection_function21cm import *

z_m_list = [30, 50,75,100,125,150,175,200]
w_list = [3.31685, 2.47355 ,1.93014, 1.60434, 1.38246, 1.21989, 1.09467, 1]

width = []
max_z = []
min_z = []
for i in range(8):
	_,_,w, z1,z2 = run_sel(w_list[i],z_m_list[i])
	width.append (w)
	max_z.append (z1)
	min_z.append (z2)
	print (w,z1,z2)

print (max_z)
print (min_z)
