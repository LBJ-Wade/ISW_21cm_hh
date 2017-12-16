#!/usr/bin/env python
# make a horizontal bar chart

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as mpatches

# log10^10A_s
val = 3.089*np.ones(3)				  # the bar lengths
pos = np.arange(3)*0.9+.8    # the bar centers on the y axis
fig,ax = plt.subplots(1)
plt.scatter (3.089,pos[2], color='green')
plt.scatter (3.089,pos[1], color='b')
plt.scatter (3.089,pos[0], color='r')
err = np.array([9.126e-2,1.98e-1,6.600e-2])
err = np.array([2.564e-2,7.72e-1,1.930e-2])
err = np.array([3.25e-2,2.564e-2,1.930e-2])
plt.barh(pos[0],val[0], xerr=err[0], ecolor='r', align='center',color='none',edgecolor='none',label ='Planck')
plt.barh(pos[1],val[1], xerr=err[1], ecolor='b', align='center',color='none',edgecolor='none', label = 'Planck + S-4')
plt.barh(pos[2],val[2], xerr=err[2], ecolor='green', align='center',color='none',edgecolor='none', label = 'Planck + S-4 + ISW-21')
plt.xlabel(r'$\mathrm{log}10^{10}A_s$')
ax.axes.get_yaxis().set_visible(False)
plt.xticks([3.057,3.089,3.121])
plt.axis([3.05,3.13,0,3])
ax.invert_yaxis()
plt.grid()
classes = [r'Planck($l<30$)',r'Planck($l<30$) + S-4',r'Planck($l<30$) + S-4 + ISW-21cm']
class_colours = ['r','b','green']
recs = []
for k in range(0,len(class_colours)):
	recs.append(mpatches.Rectangle((0,0),1,1,fc=class_colours[k]))
plt.legend(recs,classes,bbox_to_anchor=(1,1),prop={'size':12})
plt.savefig ('bar_logAs.pdf')


val = 0.06*np.ones(3)				  # the bar lengths
val = 0.06*np.ones(4)				  # the bar lengths
pos = 0.7*np.arange(4)+0.3    # the bar centers on the y axis
#pos = [0.7,1.5]
fig,ax = plt.subplots(1)
plt.scatter (.06,pos[0], color='red')
plt.scatter (.06,pos[1], color='y')
plt.scatter (.06,pos[2], color='blue')
plt.scatter (.06,pos[3], color='k')
err = np.array([9.126e-2,1.98e-1,6.600e-2])
err = np.array([2.564e-2,7.72e-1,1.930e-2])
err = np.array([2.564e-2,1.930e-2])
err = np.array([1.94e-1, 1.64e-1, 9.126e-2, 6.554e-2])
plt.barh(pos[0],val[0], xerr=err[0], ecolor='red', align='center',color='none',edgecolor='none',label ='ISW-21')
plt.barh(pos[1],val[1], xerr=err[1], ecolor='y', align='center',color='none',edgecolor='none', label = 'S-4')
plt.barh(pos[2],val[2], xerr=err[2], ecolor='blue', align='center',color='none',edgecolor='none',label ='Planck + S-4')
plt.barh(pos[3],val[3], xerr=err[3], ecolor='k', align='center',color='none',edgecolor='none', label = 'Planck + S-4 + ISW-21')
plt.xlabel(r'$\sum m_\nu$')
ax.axes.get_yaxis().set_visible(False)
#plt.xticks([-0.14,-0.04,0.06,0.16,0.26])
plt.xticks([-0.1,0,0.06,0.12,0.22])
plt.axis([-0.145,0.26,-1,3])
plt.grid()
ax.invert_yaxis()
classes = ['Planck + S-4','Planck + S-4 + ISW-21cm']
class_colours = ['b','r']
classes = [r'ISW-21cm', r'S-4',r'Planck($l<30$) + S-4',r'Planck($l<30$) + S-4 + ISW-21cm']
class_colours = ['red','y','blue','k']
recs = []
for k in range(0,len(class_colours)):
	recs.append(mpatches.Rectangle((0,0),1,1,fc=class_colours[k]))
plt.legend(recs,classes,bbox_to_anchor=(1,1),prop={'size':13})
plt.savefig ('bar_m_nu.pdf')
plt.show()
