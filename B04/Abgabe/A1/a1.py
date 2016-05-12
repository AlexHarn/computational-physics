import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15
rcParams['legend.numpoints'] = 1

color = ["g", "r"];

f = open("out.txt", "w")
for file in os.listdir():
    if file.endswith(".dat"):
        data = np.loadtxt(file, unpack=True)
        d = 2; # ist hier nicht dynamsich
        fig, ax = plt.subplots(d+1, sharex=True)

        # Winkel und ihre Ableitungen
        for i in range(0, d):
            ax[i].plot(data[0], data[i+1], color[i], label=r"$\theta_{("+str(i+1)+")}$")
            ax[i].plot(data[0], data[i+d+1], color[i]+"--", label=r"$\dot\theta_{("+str(i+1)+")}$")
            ax[i].legend()

        # Energie
        for i in range(0, d):
             ax[-1].plot(data[0], data[-2*d+i-1], color[i], label=r"$U_{("+str(i+1)+")}$")
             ax[-1].plot(data[0], data[-d+i-1], color[i]+"--", label=r"$T_{("+str(i+1)+")}$")
        ax[-1].plot(data[0], data[-1], label=r"$E_{\rm{Gesamt}}$")
        ax[-1].legend()
        plt.xlabel(r"$t$")

        # plt.show()
        fig.savefig(file[:-4]+"_zeitachse.pdf")
        plt.clf()
        f.write(file+"\n")
        f.write("Delta E: "+str(np.max(abs(data[-1] - data[-1][0])))+"\n")

        # Trajektorien Plot (die Animation ist viel schönr)
        plt.plot(data[5], -data[6]-1, "g--", label=r"$(x,y-1)_1$")
        plt.plot(data[5][0], -data[6][0]-1, "go", label=r"$(x,y-1)_1(t = 0)$")
        plt.plot(data[7], -data[8], "r--", label=r"$(x,y)_2$")
        plt.plot(data[7][0], -data[8][0], "ro", label=r"$(x,y)_2(t = 0)$")
        plt.legend(loc=3)
        plt.savefig(file[:-4]+"_trajektorie.pdf")

f.close()
