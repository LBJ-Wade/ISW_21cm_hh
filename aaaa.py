import selection_function21cm as sf
import matplotlib.pyplot as plt

z, sel,_ =sf.run_sel(3.31685,30)
print (z[0],z[-1])
plt.plot(z, sel)
plt.show()
