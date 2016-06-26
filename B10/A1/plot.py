import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15

# a
r = np.loadtxt("a.dat", unpack=True)
x = np.linspace(0, np.pi, 1e4)
plt.hist(r, bins=100, range=(0, np.pi), normed=True, label="Zufallszahlen")
plt.plot(x, 0.5*np.sin(x), linewidth=2, label="Analytisch")
plt.xlabel("Zufallszahl")
plt.ylabel("normierte Häufigkeit")
plt.savefig("a_hist.pdf")
plt.clf()

# b
r = np.loadtxt("b.dat", unpack=True)
plt.hist(r, bins=100, range=(-0.5, 1), normed=True, label="Zufallszahlen")
plt.xlabel("Zufallszahl")
plt.ylabel("normierte Häufigkeit")
plt.savefig("b_hist.pdf")
plt.clf()
