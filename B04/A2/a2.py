import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
rcParams['figure.figsize'] = 15, 15

#########################
# Teilaufgabe a
#########################

a1t, a1th1, a1th2, a1tp1, a1tp2 = np.genfromtxt("a1.dat", unpack = True)
a2t, a2th1, a2th2, a2tp1, a2tp2 = np.genfromtxt("a2.dat", unpack = True)
a3t, a3th1, a3th2, a3tp1, a3tp2 = np.genfromtxt("a3.dat", unpack = True)


plt.subplot(311)
plt.title('periodisch')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$\dot{\theta_1}$')
plt.plot(a1th1, a1tp1, 'r-')
plt.subplot(312)
plt.title('quasiperiodisch')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$\dot{\theta_1}$')
plt.plot(a2th1, a2tp1, 'r-')
plt.subplot(313)
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$\dot{\theta_1}$')
plt.title('chaotisch')
plt.plot(a3th1, a3tp1, 'r-')
plt.tight_layout()
plt.savefig("teila.pdf")
plt.clf()

#########################
# Teilaufgabe b
#########################

b2t, b2th1, b2th2, b2tp1, b2tp2 = np.genfromtxt("b2.dat", unpack = True)
b3t, b3th1, b3th2, b3tp1, b3tp2 = np.genfromtxt("b3.dat", unpack = True)

plt.subplot(221)
plt.title('quasiperiodisch')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$|\dot{\theta_1}^u - \dot{\theta_1}^g|$') # u = ungestört, g = gestört
plt.plot(b2th1,abs(a2tp1-b2tp1)) # quasiperiodisch
plt.subplot(222)
plt.title('quasiperiodisch')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$|\dot{\theta_2}^u - \dot{\theta_2}^g|$') # u = ungestört, g = gestört
plt.plot(b2th2,abs(a2tp2-b2tp2))
plt.subplot(223)
plt.title('chaotisch')
plt.plot(b3th1,abs(a3tp1-b3tp1)) # chaotisch
plt.xlabel(r'$\theta_2$')
plt.ylabel(r'$|\dot{\theta_1}^u - \dot{\theta_1}^g|$') # u = ungestört, g = gestört
plt.subplot(224)
plt.title('chaotisch')
plt.plot(b3th2,abs(a3tp2-b3tp2))
plt.xlabel(r'$\theta_2$')
plt.ylabel(r'$|\dot{\theta_2}^u - \dot{\theta_2}^g|$') # u = ungestört, g = gestört
plt.tight_layout()
plt.savefig("teilb.pdf")
plt.clf()

#########################
# Teilaufgabe c
#########################

c1th1, c1tp1 = np.genfromtxt("c1.dat", unpack = True)
c2th1, c2tp1 = np.genfromtxt("c2.dat", unpack = True)
c3th1, c3tp1 = np.genfromtxt("c3.dat", unpack = True)

plt.subplot(311)
plt.title('3 Joule')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$\dot{\theta_1}$')
plt.plot(c1th1, c1tp1, 'b.')
plt.subplot(312)
plt.title('10 Joule')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$\dot{\theta_1}$')
plt.plot(c2th1, c2tp1, 'b.')
plt.subplot(313)
plt.title('20 Joule')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$\dot{\theta_1}$')
plt.plot(c3th1, c3tp1, 'b.')
plt.tight_layout()
plt.savefig("teilc.pdf")
plt.clf()