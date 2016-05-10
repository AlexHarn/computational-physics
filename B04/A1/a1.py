import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15


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

         #Energie
        for i in range(0, d):
             ax[-1].plot(data[0], data[-2*d+i-1], label=r"$U_{("+str(i)+")}$")
             ax[-1].plot(data[0], data[-d+i-1], label=r"$T_{("+str(i)+")}$")
        ax[-1].plot(data[0], data[-1], label=r"$E_{\rm{Gesamt}}$")
        ax[-1].legend()
        plt.xlabel(r"$t$")

        #plt.show()
        fig.savefig(file[:-4]+"_zeitachse.pdf")
        plt.clf()
        f = open("out.txt", "w")
        f.write(file+"\n")
        f.write("Delta E: "+str(data[-1][-1] - data[-1][0])+"\n")
        f.close()
