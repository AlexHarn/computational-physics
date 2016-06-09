import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15

# A1
names = [ "1b1", "1b2", "1b3", "1b4" ]
Ns = [ 6074, 254, 2147483646, 2147483646 ]
# Die Frage ist, ob es viel Sinn macht immer so viele zu plotten

for name, N in zip(names, Ns):
    x = np.loadtxt(name+".dat", unpack=True)
    plt.hist(x) # 10 Bins sind default
    plt.xlabel("Zufallszahl")
    plt.ylabel("Häufigkeit")
    plt.savefig(name+"_hist.pdf")
    plt.clf()

    plt.plot(x[:N][1::2], x[:N][::2], ".")
    plt.xlabel(r"$r_n$")
    plt.ylabel(r"$r_{n-1}$")
    plt.savefig(name+"_korr.pdf")
    plt.clf()

# A2
names = [ "2a", "2b" ]
ps = [ lambda x: np.exp(-x**2/2)/( np.sqrt(2*np.pi) ),
       lambda x: np.exp(-x**2/2)/( np.sqrt(2*np.pi) )]

for name, p in zip(names, ps):
    n = np.loadtxt(name+".dat", unpack=True)
    plt.hist(n, normed=True, alpha=0.5, label="Generierte Zahlen") # 10 Bins sind default
    x = np.linspace(min(n), max(n), 1e3)
    plt.plot(x, p(x), linewidth=2, label="Analytisch")
    plt.legend()
    plt.xlabel(r"$x$")
    plt.ylabel(r"$p(x)$")
    plt.savefig(name+".pdf")
    plt.clf()
