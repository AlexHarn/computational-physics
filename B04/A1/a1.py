import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15

f = open("out.txt", "w")
for file in os.listdir():
    if file.endswith(".dat"):
        data = np.loadtxt(file, unpack=True)
        d = round(( len(data) - 6 )/2) # 6: Energien und Zeit
                                       # 2: theta 1, 2
        fig, ax = plt.subplots(d+1, sharex=True)

        # Winkel und ihre Ableitungen
        for i in range(0, d):
            ax[i].plot(data[0], data[i+1], label=r"$\theta_{("+str(i)+")}$")
            ax[i].plot(data[0], data[i+d+1], label=r"$\dot\theta_{("+str(i)+")}$")
            ax[i].legend()

        # Energie
        # for i in range(0, d):
            # ax[-1].plot(data[0], data[-2*d+i-1], label=r"$U_{("+str(i)+")}$")
            # ax[-1].plot(data[0], data[-d+i-1], label=r"$T_{("+str(i)+")}$")
        ax[-1].plot(data[0], data[-1], label=r"$E_{\rm{Gesamt}}$")
        # ax[-1].plot(data[0], E - E_0, label=r"$\Delta E_{\rm{Gesamt}}$")
        ax[-1].legend()
        plt.xlabel(r"$t$")

        plt.show()
        # fig.savefig(file[:-4]+"_zeitachse.pdf")
        plt.clf()

        f.write(file+"\n")
        f.write("Delta E: "+str(data[-1][-1] - data[-1][0])+"\n")

        # # Ein 2D Plott der Bahn darf hier natürlich nicht fehlen
        # sun = plt.Circle((0,0), 0.1, color='y')
        # fig = plt.gcf()
        # ax = plt.gca()
        # ax.add_artist(sun)

        # ax.set_xlim((-1.5,1.5))
        # ax.set_ylim((-1.5,1.5))
        # plt.grid(True)
        # plt.axes().set_aspect('equal', 'datalim')

        # ax.plot(data[1][0], data[2][0], "g.", label="Start")
        # ax.plot(data[1], data[2], ".", label=r"$xy$-Bahn", markersize=0.5)
        # ax.legend()

        # # plt.show()
        # plt.savefig(file[:-4]+"_bahn.png") # pdf versagt hier
        # plt.clf()

    # elif file.endswith(".dat1"):
        # n, dr = np.loadtxt(file, unpack=True)

        # plt.plot(n, dr, ".", label="Ortsfehler")
        # plt.legend()
        # plt.xlabel(r"$n$")
        # # plt.show()
        # plt.savefig(file[:-4]+"pdf")
        # plt.clf()
f.close()
