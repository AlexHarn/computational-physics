import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15

names = [ "b1", "b2", "b3", "b4" ]
Ns = dict(zip(names, [ 6074, 254, 2147483646, 2147483646 ]))
# Die Frage ist, ob es viel Sinn macht immer so viele zu plotten

for name in names:
    x = np.loadtxt(name+".dat", unpack=True)
    plt.hist(x) # 10 Bins sind default
    plt.xlabel("Zufallszahl")
    plt.ylabel("Häufigkeit")
    plt.savefig(name+"_hist.pdf")
    plt.clf()

    plt.plot(x[:Ns[name]][1::2], x[:Ns[name]][::2], ".")
    plt.xlabel(r"$r_n$")
    plt.ylabel(r"$r_{n-1}$")
    plt.savefig(name+"_korr.pdf")
    plt.clf()
