import matplotlib.pyplot as plt
import numpy as np

Ns = [ 2, 5, 10 ]
for N in Ns:
    theta, E_f, E_af = np.loadtxt("a1_N="+str(N)+".dat", unpack=True,
    skiprows=1);
    plt.plot(theta, E_f)
    plt.plot(theta, E_af)
    plt.show()

