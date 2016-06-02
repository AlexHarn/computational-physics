import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
np.seterr(divide='ignore', invalid='ignore')
rcParams['figure.figsize'] = 15, 15
delta = 0.05
L = 1
n = 1
scale = np.arange(delta/2, L+delta/2, delta)

def f(x,y):
	return ( np.sin(n*np.pi*x) * np.sinh(n*np.pi*y) )/( np.sin(n*np.pi*0.5) * np.sinh(n*np.pi) )
theorie = np.zeros((20, 20))
for i in range(0, 20, 1):
	for j in range(0, 20, 1):
		theorie[i][j] = f(0.025+j*0.05, 0.025+i*0.05)

for file in os.listdir():
    if file.endswith(".dat"):
        pre = file[:file.index("_")]
        if file[file.index("_")+1:-4] == "phi":
            phi = np.loadtxt(file);
            ex = np.loadtxt(pre+"_ex.dat")
            ey = np.loadtxt(pre+"_ey.dat")

            fig = plt.figure();
            if pre == "b":
                ax = fig.add_subplot(1, 2, 1);
                image1 = ax.imshow(phi, origin="lower",
                                  extent=[0, L, 0, L]);
                # fig.colorbar(image, ax=ax)
                ax.streamplot(scale, scale, ex, ey, color="Black")
                ax = fig.add_subplot(1, 2, 2)
                image2 = ax.imshow(theorie, origin="lower",
                                  extent=[0, 1, 0, 1])
                # fig.colorbar(image, ax=ax)
                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                fig.colorbar(image1, cax=cbar_ax)
                fig.savefig(pre+".pdf");
            else:
                ax = fig.add_subplot(1, 1, 1);
                image = ax.imshow(phi, origin="lower",
                                  extent=[0, L, 0, L]);
                fig.colorbar(image, ax=ax)
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
