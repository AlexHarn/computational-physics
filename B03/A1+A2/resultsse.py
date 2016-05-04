import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
import numpy as np

plt.rcParams['figure.figsize'] = (20, 10)

t,x,v = np.genfromtxt("data2.dat", unpack = True)

plt.title('Blatt 03, Aufgabe 2')
plt.plot(t,x, 'bx', label = 'Ort')
plt.plot(t,v, 'rx', label = 'Geschwindigkeit')
plt.xlabel(r'Zeit')
# plt.ylabel(r'Ort')

plt.legend()
plt.grid()
# plt.savefig('Aufgabe2.pdf')
plt.show()