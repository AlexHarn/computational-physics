import numpy as np
import matplotlib.pyplot as plt
import os

for file in os.listdir():
    if file.endswith(".dat"):
        r, g = np.genfromtxt(file, unpack=True)
        plt.xlim(min(r), max(r))
        plt.plot(r, g)
        plt.axhline(color="k")
        plt.xlabel(r"$r$")
        plt.ylabel(r"$g(r)$")
        # plt.show()
        plt.savefig(file[:-4]+".pdf")
        plt.clf()
