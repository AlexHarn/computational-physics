import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import rcParams
import os
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

if(os.stat("c1_1.dat").st_size != 0):
    c11th1, c11tp1 = np.genfromtxt("c1_1.dat", unpack = True)
if(os.stat("c1_2.dat").st_size != 0):
    c12th1, c12tp1 = np.genfromtxt("c1_2.dat", unpack = True)
if(os.stat("c1_3.dat").st_size != 0):
    c13th1, c13tp1 = np.genfromtxt("c1_3.dat", unpack = True)
if(os.stat("c1_4.dat").st_size != 0):
    c14th1, c14tp1 = np.genfromtxt("c1_4.dat", unpack = True)
if(os.stat("c2_1.dat").st_size != 0):
    c21th1, c21tp1 = np.genfromtxt("c2_1.dat", unpack = True)
if(os.stat("c2_2.dat").st_size != 0):
    c22th1, c22tp1 = np.genfromtxt("c2_2.dat", unpack = True)
if(os.stat("c2_3.dat").st_size != 0):
    c23th1, c23tp1 = np.genfromtxt("c2_3.dat", unpack = True)
if(os.stat("c2_4.dat").st_size != 0):
    c24th1, c24tp1 = np.genfromtxt("c2_4.dat", unpack = True)
if(os.stat("c3_1.dat").st_size != 0):
    c31th1, c31tp1 = np.genfromtxt("c3_1.dat", unpack = True)
if(os.stat("c3_2.dat").st_size != 0):
    c32th1, c32tp1 = np.genfromtxt("c3_2.dat", unpack = True)
if(os.stat("c3_3.dat").st_size != 0):
    c33th1, c33tp1 = np.genfromtxt("c3_3.dat", unpack = True)
if(os.stat("c3_4.dat").st_size != 0):
    c34th1, c34tp1 = np.genfromtxt("c3_4.dat", unpack = True)

plt.subplot(311)
plt.title('3 Joule')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$\dot{\theta_1}$')
if(os.stat("c1_1.dat").st_size != 0):
    plt.plot(c11th1, c11tp1, 'r.')
if(os.stat("c1_2.dat").st_size != 0):
    plt.plot(c12th1, c12tp1, 'g.')
if(os.stat("c1_3.dat").st_size != 0):
    plt.plot(c13th1, c13tp1, 'b.')
if(os.stat("c1_4.dat").st_size != 0):
    plt.plot(c14th1, c14tp1, 'y.')
plt.subplot(312)
plt.title('10 Joule')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$\dot{\theta_1}$')
if(os.stat("c2_1.dat").st_size != 0):
    plt.plot(c21th1, c21tp1, 'r.')
if(os.stat("c2_2.dat").st_size != 0):
    plt.plot(c22th1, c22tp1, 'g.')
if(os.stat("c2_3.dat").st_size != 0):
    plt.plot(c23th1, c23tp1, 'b.')
if(os.stat("c2_4.dat").st_size != 0):
    plt.plot(c24th1, c24tp1, 'y.')
plt.subplot(313)
plt.title('20 Joule')
plt.xlabel(r'$\theta_1$')
plt.ylabel(r'$\dot{\theta_1}$')
if(os.stat("c3_1.dat").st_size != 0):
    plt.plot(c31th1, c31tp1, 'r.')
if(os.stat("c3_2.dat").st_size != 0):
    plt.plot(c32th1, c32tp1, 'g.')
if(os.stat("c3_3.dat").st_size != 0):
    plt.plot(c33th1, c33tp1, 'b.')
if(os.stat("c3_4.dat").st_size != 0):
    plt.plot(c34th1, c34tp1, 'y.')
plt.tight_layout()
plt.savefig("teilc.pdf")
plt.clf()