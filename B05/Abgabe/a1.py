import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15
rcParams['legend.numpoints'] = 1

labels = [ r"$v_{\rm{S}}$", r"$T$", r"$E_{\rm{Pot}}$", r"$E_{\rm{Kin}}$",
          r"$E_{\rm{Gesamt}}$" ]

for file in os.listdir():
    if file.endswith(".dat"):
        if file[file.index("_")+1:file.index(".")]=="savedata":
            data = list(np.loadtxt(file, unpack=True))
            fig, ax = plt.subplots(5, sharex=True)
            data.append(data[-1] + data[-2])
            data = np.asarray(data)
            for i in range(len(data)-1):
                ax[i].plot(data[0], data[i+1], label=labels[i])
                ax[i].legend()
            plt.xlabel(r"$t$")
            plt.savefig(file[:-4]+".pdf")
            plt.clf()
            # plt.show()
        else:
            x, y = np.loadtxt(file, unpack=True)
            plt.plot(x, y)
            # plt.show()
            plt.savefig(file[:-4]+".pdf")
            plt.clf()
