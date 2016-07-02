import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15


h, m = np.loadtxt("b.dat", unpack=True)
x = np.linspace(-5, 5, 1e4)
plt.plot(h, m, color='blue', label="Numerisch")
plt.plot(x, np.tanh(x), color='red', ls='--', label="Analytisch")
plt.axhline(ls='--', color='k')
plt.axvline(ls='--', color='k')
plt.xlabel("Magnetfeld")
plt.ylabel("Magnetisierung")
plt.legend(loc="best")
plt.savefig("b.pdf")
plt.clf()
