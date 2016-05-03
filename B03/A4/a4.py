import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15

f = open("out.txt", "w")
for file in os.listdir():
    if file.endswith(".dat"):
        data = np.loadtxt(file, unpack=True)
        d = round(( len(data) - 4 )/3) # 4: Energien, Fläche, Drehimpuls
                                       # 3: r,v und Lenz Runge
        fig, ax = plt.subplots(d-1+4)

        # Orte und Geschwindigkeiten
        for i in range(0, d-1):
            ax[i].plot(data[0], data[i+1], label=r"$r^{("+str(i)+")}$")
            ax[i].plot(data[0], data[i+d+1], label=r"$v^{("+str(i)+")}$")
            ax[i].legend()

        # Lenz-Runge-Vektor
        # man erkennt nix wenn man alle 3 Komponenten auf einmal darstellt
        # darum nur mal die x Komponente repräsentativ z ist eh 0 und y verhält
        # sich wie x
        ax[-4].plot(data[0], data[-7], label=r"$\Lambda_x$")
        # ax[-4].plot(data[0], data[-6], label=r"$\Lambda_y$")
        # ax[-4].plot(data[0], data[-5], label=r"$\Lambda_z$")
        ax[-4].legend()

        # Drehimpuls
        ax[-3].plot(data[0], data[-4], label=r"$L(t)$")
        ax[-3].legend()

        # Überstrichene Fläche
        ax[-2].plot(data[0], data[-3], label=r"$A(t)$") # sollte Gerade sein
        ax[-2].legend()

        # Energie
        ax[-1].plot(data[0], data[-2], label=r"$U$")
        ax[-1].plot(data[0], data[-1], label=r"$T$")
        E = data[-1]+data[-2]
        E_0 = data[-1][0]+data[-2][0]
        ax[-1].plot(data[0], E, label=r"$E_{\rm{Gesamt}}$")
        # ax[-1].plot(data[0], E - E_0, label=r"$\Delta E_{\rm{Gesamt}}$")
        ax[-1].legend()
        plt.xlabel(r"$t$")

        # plt.show()
        fig.savefig(file[:-4]+"_zeitachse.pdf")
        plt.clf()

        f.write(file+"\n")
        f.write("Delta E: "+str(E[-1] - E[0])+"\n")
        f.write("Delta L: "+str(data[-4][-1] - data[-4][0])+"\n")
        f.write("Delta LR_x: "+str(data[-7][-1] - data[-7][0])+"\n")
        f.write("Delta LR_y: "+str(data[-6][-1] - data[-6][0])+"\n")
        f.write("Delta LR_z: "+str(data[-5][-1] - data[-5][0])+"\n")

        # Ein 2D Plott der Bahn darf hier natürlich nicht fehlen
        sun = plt.Circle((0,0), 0.1, color='y')
        fig = plt.gcf()
        ax = plt.gca()
        ax.add_artist(sun)

        ax.set_xlim((-1.5,1.5))
        ax.set_ylim((-1.5,1.5))
        plt.grid(True)
        plt.axes().set_aspect('equal', 'datalim')

        ax.plot(data[1][0], data[2][0], "g.", label="Start")
        ax.plot(data[1], data[2], ".", label=r"$xy$-Bahn", markersize=0.5)
        ax.legend()

        # plt.show()
        plt.savefig(file[:-4]+"_bahn.png") # pdf versagt hier
        plt.clf()

    elif file.endswith(".dat1"):
        n, dr = np.loadtxt(file, unpack=True)

        plt.plot(n, dr, ".", label="Ortsfehler")
        plt.legend()
        plt.xlabel(r"$n$")
        # plt.show()
        plt.savefig(file[:-4]+"pdf")
        plt.clf()
f.close()
