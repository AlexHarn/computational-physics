import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15

# b errors
n, pi = np.loadtxt("b_errors.dat", unpack=True)
plt.loglog(n, pi, ".")
plt.xlabel(r"$N$")
plt.ylabel(r"$\left|\pi - pi_N\right|$")
plt.savefig("b_errors.pdf")
plt.clf()

# b histogram
r = np.loadtxt("b_hist.dat", unpack=True)
plt.hist(r, bins=20, label="Daten")
plt.xlabel(r"$\pi$")
plt.ylabel("Häufigkeit")
plt.savefig("b_hist.pdf")
plt.clf()

# c plot
a, b, area = np.loadtxt("c.dat", unpack=True)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(a, b, area, cmap=cm.jet)
ax.set_xlabel(r"$a$")
ax.set_ylabel(r"$b$")
ax.set_zlabel("Fläche")

# plt.show()
plt.savefig("c.pdf")
