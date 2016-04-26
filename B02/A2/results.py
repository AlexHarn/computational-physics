import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
import numpy as np
plt.rcParams['figure.figsize'] = (20, 10)
#plt.rc('text', usetex=True)

x,y = np.genfromtxt("data2a+b.dat", unpack = True)

plt.subplot(121)
plt.title('Blatt 02, Aufgabe 2a, 2b')
plt.plot(x,y, 'bx')
plt.xlabel(r'$x/a$')
plt.ylabel(r'$\phi\cdot 4\pi\epsilon_0 / \rho$')

plt.axvspan(0, 1, color = 'grey')

plt.grid()


x,y = np.genfromtxt("data2c.dat", unpack = True)

plt.subplot(122)
plt.title('Blatt 02, Aufgabe 2c')
plt.plot(x,y, 'bx')
plt.xlabel(r'$x/a$')
plt.ylabel(r'$\phi\cdot 4\pi\epsilon_0 / \rho$')

plt.axvspan(0, 1, color = 'grey')

plt.grid()
plt.savefig('Aufgabe2.pdf')
plt.show()
