# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
j, rho, u, e, P =np.loadtxt("UpwindGodunov_step_5.dat", unpack=True)
#j, rho, u, e, P = data[:,0], data[:,1] , data[:,2] , data[:,3] , data[:,4]
plt.plot(j,rho) 
plt.xlabel(u"Posición") 
plt.ylabel("Densidad") 
plt.title(u"Densidad en función de la posición para tiempo final (1 s)") 
plt.savefig("x-rho.png", dpi=400) 
plt.clf()

plt.plot(j,u) 
plt.xlabel(u"Posición") 
plt.ylabel("Velocidad") 
plt.title(u"Velocidad en función de la posición para tiempo final (1 s)") 
plt.savefig("x-u.png", dpi=400) 
plt.clf()

plt.plot(j,e) 
plt.xlabel(u"Posición") 
plt.ylabel(u"Energía") 
plt.title(u"Energía en función de la posición para tiempo final (1 s)") 
plt.savefig("x-e.png", dpi=400) 
plt.clf()

plt.plot(j,P) 
plt.xlabel(u"Posición") 
plt.ylabel(u"Presión") 
plt.title(u"Presión en función de la posición para tiempo final (1 s)") 
plt.savefig("x-P.png", dpi=400) 
plt.clf()
