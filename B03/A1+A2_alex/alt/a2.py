import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15

for file in os.listdir():
    if file.endswith(".dat"):
        data = np.loadtxt(file, unpack=True)

        d = round(( len(data) - 3 )/2)

        fig, ax = plt.subplots(d+1, sharex=True)
        for i in range(0, d):
            ax[i].plot(data[0], data[i+1], label=r"$r^{("+str(i)+")}$")
            ax[i].plot(data[0], data[i+d+1], label=r"$v^{("+str(i)+")}$")
            ax[i].legend()

        ax[d].plot(data[0], data[-2], label=r"$U$")
        ax[d].plot(data[0], data[-1], label=r"$T$")
        E = data[-1]+data[-2]
        E_0 = data[-1][0]+data[-2][0]
        ax[d].plot(data[0], E, label=r"$E_{\rm{Gesamt}}$")
        # ax[d].plot(data[0], E - E_0, label=r"$\Delta E_{\rm{Gesamt}}$")
        # ax[d].set_ylim(-0.1*E[0], E[0]+0.1*E[0])
        ax[d].legend()

        print("Delta E: "+str(E[-1] - E[0]))

        # plt.show()
        plt.savefig(file[:-4]+".pdf")
