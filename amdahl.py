import numpy as np
import matplotlib.pyplot as plt

p1 = 0.1985
p2 = 0.5942
N = np.arange(2,10**3)
speed = 1/( (1-p1-p2) + p1/2. + p2/N)
N = np.array([1] + list(N))
speed = np.array([1] + list(speed))
plt.plot (N, speed)
plt.xscale ('log')
plt.scatter (2, 1.5032479249368458)
plt.axis ([1, 10**3, 0.5, 3.5])
plt.xlabel ('Number of Processors',size = 15)
plt.ylabel ('Speedup', size = 15)
plt.title ('Amdahl\'s Law', size = 15)
plt.savefig ('amdahl.pdf')
plt.show()
