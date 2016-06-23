from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['figure.figsize'] = 15, 10


dataX = np.loadtxt("2d_x.dat")
dataE = np.loadtxt("2d_e.dat")

nX = dataX.T[0][1:]
nE = dataE.T[0][1:]
absX = [np.linalg.norm(dataX[i][1:] - dataX[i-1][1:]) for i in range(1, len(nX)+1)]
absE = [np.linalg.norm(dataE[i][1:] - dataE[i-1][1:]) for i in range(1, len(nE)+1)]

fig, ax = plt.subplots()

# ax.axhline(color="k")
ax.plot(nX, absX, ".", label="Ortsdarstellung")
ax.plot(nE, absE, ".", label="Besetzungszahldarstellung")
plt.xlim(9, 51)
ax.set_ylim(-50, 1400)
ax.set_xlabel(r"$N$")
ax.set_ylabel(r"$\left|V_N - V_{N-1}\right|$")
ax.legend()

axins = inset_axes(ax, 2, 2, loc=5)
mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")
axins.plot(nX, absX, ".")
axins.plot(nE, absE, ".")
axins.set_xlim(40, 50)
axins.set_ylim(-1, 6)
# plt.yticks(visible=False)
# plt.xticks(visible=False)

axins1 = inset_axes(ax, 2, 3, loc=9)
mark_inset(ax, axins1, loc1=3, loc2=4, fc="none", ec="0.5")
axins1.plot(nE, absE, "g.")
axins1.set_xlim(10, 20)
axins1.set_ylim(-1, 17)

plt.savefig("d.pdf")
