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
names = [ "2a", "2b", "2c", "2d" ]
ps = [ lambda x: np.exp(-x**2/2)/( np.sqrt(2*np.pi) ),
       lambda x: np.exp(-x**2/2)/( np.sqrt(2*np.pi) ),
       lambda x: np.sin(x)/2,
       lambda x: 3*x**2]

for name, p in zip(names, ps):
    n = np.loadtxt(name+".dat", unpack=True)
    plt.hist(n, 100, normed=True, alpha=0.5, label="Generierte Zahlen") # 10 Bins sind default
    x = np.linspace(min(n), max(n), 1e3)
    plt.plot(x, p(x), linewidth=2, label="Analytisch")
    plt.legend()
    plt.xlabel(r"$x$")
    plt.ylabel(r"$p(x)$")
    plt.savefig(name+".pdf")
    plt.clf()

def f(x):
    return np.exp(-x**2/2)/( np.sqrt(2*np.pi) )
    
def g(x):
    return np.sin(x)/2
    
x = np.linspace(0, np.pi)
plt.plot(x, 1.5*f(x-np.pi/2))
plt.plot(x, g(x))
plt.savefig("2c_test.pdf")