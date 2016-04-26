import matplotlib.pyplot as plt
import numpy as np
#from tabulate import tabulate
from matplotlib.pyplot import rc
plt.rcParams['figure.figsize'] = (20, 10)
Ns = [ 2, 5, 10 ]
fig = plt.figure()
for i, N in enumerate(Ns):
    ax = fig.add_subplot(1, len(Ns), i + 1)
    theta, E_f, T_f_a, T_f_k, E_af, T_af_a, T_af_k = np.loadtxt("N="+str(N)+".dat", unpack=True,
    skiprows=1)
    ax.plot(theta, T_af_a, 'r-', label="über die Ableitung")
    ax.plot(theta, T_af_k, 'gx', label="über das Kreuzprodukt")
    ax.legend(loc=1)
    ax.set_title(r"$N = {}$ für den antiferromgn Fall".format(N))
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 0.15)
    ax.set_xlabel(r"$\theta")
    ax.set_ylabel(r"$E({}, \theta)")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$E({}, \theta)$".format(N))
    ax.grid()
    #f = open("N="+str(N)+".tex", "w")
    #f.write(tabulate(np.loadtxt("N="+str(N)+".dat", skiprows=1),
    #                 tablefmt="latex"))
    #f.close()

plt.show();
fig.savefig('Aufgabe1b2.pdf')
