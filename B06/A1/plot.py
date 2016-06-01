import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
np.seterr(divide='ignore', invalid='ignore')
rcParams['figure.figsize'] = 15, 15
delta = 0.05
L = 1

for file in os.listdir():
    if file.endswith(".dat"):
        pre = file[:file.index("_")]
        if file[file.index("_")+1:-4] == "phi":
            phi = np.loadtxt(file);
            ex = np.loadtxt(pre+"_ex.dat")
            ey = np.loadtxt(pre+"_ey.dat")

            fig = plt.figure();
            ax = fig.add_subplot(1, 1, 1);
            image = ax.imshow(phi, origin="lower",
                              extent=[0, L, 0, L]);
            fig.colorbar(image, ax=ax)
            scale = np.arange(delta/2, L+delta/2, delta)
            # ax.quiver(scale, scale, ex, ey)
            ax.streamplot(scale, scale, ex, ey, color="Black")
            fig.savefig(pre+".pdf");
            # plt.show()

            # Weil es explizit gefordert ist auch die Beträge:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            im = ax.imshow(np.sqrt(ex**2 + ey**2), origin="lower",
                           extent=[0, L, 0, L])
            fig.colorbar(im, ax=ax)
            fig.savefig(pre+"_absE.pdf");
