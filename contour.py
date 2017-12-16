import matplotlib.pyplot as plt
import numpy as np
from matplotlib import mlab, cm

# Default delta is large because that makes it fast, and it illustrates
# the correct registration between image and contours.
delta = 0.00001

extent = (-3, 4, -4, 3)

x = np.arange(-0.002, 0.002, delta)
y = np.arange(-0.001, 0.001, delta)
X, Y = np.meshgrid(x, y)
Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0,0.9)
Z1 = mlab.bivariate_normal(X, Y, 0.000692610679333, 3.7124235704e-05, 0.0, 0.0,0.0166808041712)
Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
Z = (Z1 - Z2) * 10
Z = Z1
print (np.amax(Z))
# Boost the upper limit to avoid truncation errors.
levels = np.arange(-2.0, 1.601, 0.4)
levels = np.arange(0,0.1,0.05)
levels1 = [np.amax(Z1)*np.e**(-3.09),np.amax(Z1)*np.e**(-1.15),np.amax(Z1)]
levels2 = [0.02,0.05,np.amax(Z2)]
norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
cmap = cm.PRGn

fig = plt.figure()
fig.subplots_adjust(hspace=0.3)


plt.subplot(2, 2, 1)

cset1 = plt.contourf(X, Y, Z1, levels1,
                     cmap=cm.get_cmap(cmap, len(levels) - 1), norm=norm)
#cset2 = plt.contourf(X, Y, Z2, levels2,
#                     cmap=cm.get_cmap(cmap, len(levels) - 1), norm=norm)
# It is not necessary, but for the colormap, we need only the
# number of levels minus 1.  To avoid discretization error, use
# either this number or a large number such as the default (256).

# If we want lines as well as filled regions, we need to call
# contour separately; don't try to change the edgecolor or edgewidth
# of the polygons in the collections returned by contourf.
# Use levels output from previous call to guarantee they are the same.

cset2 = plt.contour(X, Y, Z, cset1.levels, colors='k')

# We don't really need dashed contour lines to indicate negative
# regions, so let's turn them off.

#for c in cset2.collections:
#	    c.set_linestyle('solid')

# It is easier here to make a separate call to contour than
# to set up an array of colors and linewidths.
# We are making a thick green line as a zero contour.
# Specify the zero level as a tuple with only 0 in it.

#cset3 = plt.contour(X, Y, Z, (0,), colors='g', linewidths=2)
plt.title('Filled contours')
plt.colorbar(cset1)
plt.show()
