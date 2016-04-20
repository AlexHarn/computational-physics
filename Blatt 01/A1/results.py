import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

Ns = [ 2, 5, 10 ]
fig = plt.figure()
for i, N in enumerate(Ns):
    ax = fig.add_subplot(1, len(Ns), i + 1)
    theta, E_f, T_f_a, T_f_k, E_af, T_af_a, T_af_k = np.loadtxt("N="+str(N)+".dat", unpack=True,
    skiprows=1)
    ax.plot(theta, E_f, label="Ferro")
    ax.plot(theta, E_af, label="Antiferro")
    ax.legend(loc=1)
    ax.set_title(r"$N = {}$".format(N))
    ax.set_xlim(0, 360)
    ax.set_ylim(-0.45, 0.45)
    ax.set_xlabel(r"$\theta")
    ax.set_ylabel(r"$E({}, \theta)")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$E({}, \theta)$".format(N))
    ax.grid()
    f = open("N="+str(N)+".tex", "w")
    # f.write(tabulate({"$\theta$" : theta,
                      # "$"},
                     # tablefmt="latex"))

plt.show();
# fig.savefig('Aufgabe 1a.pdf') # sieht in dem Format nicht gut aus habe keine
# Zeit mehr das zu fixen also einfach manuell speichern

