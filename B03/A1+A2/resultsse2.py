import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
import numpy as np

plt.rcParams['figure.figsize'] = (20, 10)

t,x,y,z,vx,vy,vz, Ek, Ep = np.genfromtxt("data2a2.dat", unpack = True)

plt.title('Blatt 03, Aufgabe 2a2')
plt.subplot(411)
plt.plot(t,x, 'bx', label = 'Ort')
plt.plot(t,vx, 'rx', label = 'Geschwindigkeit')
plt.xlabel(r't')
plt.ylabel(r'x')
plt.grid() #
plt.legend()
plt.subplot(412)
plt.plot(t,y, 'bx', label = 'Ort')
plt.plot(t,vy, 'rx', label = 'Geschwindigkeit')
plt.xlabel(r't')
plt.ylabel(r'y')
plt.grid() #
plt.legend()
plt.subplot(413)
plt.plot(t,z, 'bx', label = 'Ort')
plt.plot(t,vz, 'rx', label = 'Geschwindigkeit')
plt.xlabel(r't')
plt.ylabel(r'z')
plt.grid() #
plt.legend()
plt.subplot(414)
plt.plot(t,Ek, 'bx', label = 'Ek')
plt.plot(t,Ep, 'gx', label = 'Ep')
plt.plot(t,Ek+Ep, 'rx', label = 'Eg')
plt.xlabel(r't')
plt.ylabel(r'E')
plt.grid() #
plt.legend()
plt.show()
plt.savefig('Aufgabe2a2.pdf')
