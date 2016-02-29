#This script is based on the one provided in An Introduction to Scientific Computing (Springer, 2005) by I. Danaila, P. Joly, S. M. Kaber & M. Postel

def mach(x, gamma, aL, pL, rhoL, aR, pR, rhoR, dum1, dum2, dum3):
	return (x-1./x)-aL/dum2*(1-(pR/pL*(dum1*gamma*x*x-dum2))**(dum3/gamma))

def exacta(x,x0,t, gamma, aL, pL, rhoL, aR, pR, rhoR, dum1, dum2, dum3):
	
	import numpy as np
	from scipy.optimize import fsolve
	from math import sqrt
#	f = open("exact.dat", "w")
#	f.write( str(yEst)  )      # str() converts to string

	M=len(x)
	uex=np.zeros((3,M), float)# this might be better transposed

	Ms=fsolve(mach,2.,args=( gamma, aL, pL, rhoL, aR, pR, rhoR, dum1, dum2, dum3))
#	fprintf('Shock Mach number Ms=%f \n',Ms)  #irrelevant?

	dumm=Ms*Ms
 
	p1  = pR*(dum1*gamma*dumm-dum2)
	rho1= rhoR/(dum1/dumm+dum2)
	U1  = dum1*(Ms-1./Ms)
	a1  = sqrt(gamma*p1/rho1)
	
	a2  = aL-dum3*U1
	rho2= rhoL*(p1/pL)**(1./gamma)
 
	x1  = x0-aL*t
	x2  = x0+(U1-a2)*t
	x3  = x0+U1*t
	x4  = x0+Ms*t
	
	idum = (x<=x1).ravel().nonzero()  #indexes are transposed
	uex[0,idum] = rhoL
	uex[1,idum] = 0.
	uex[2,idum] = pL

	idum = np.array((x1<x) & (x<=x2)).ravel().nonzero()
	uex[1,idum] = dum1*(aL+     (x[idum]-x0)/t)
       	adet =        dum1*(aL-dum3*(x[idum]-x0)/t)
	uex[2,idum] = pL*(adet/aL)**(2.*gamma/(gamma-1.))
	uex[0,idum] = gamma*(uex[2,idum])/(adet*adet)

	idum = np.array((x2<x) & (x<=x3)).ravel().nonzero()
	uex[0,idum] = rho2
	uex[1,idum] = U1
	uex[2,idum] = p1

	idum = np.array((x3<x) & (x<=x4)).ravel().nonzero()
	uex[0,idum] = rho1
	uex[1,idum] = U1
	uex[2,idum] = p1

	idum = (x4<x).ravel().nonzero()
	uex[0,idum] = rhoR	#1.	#rhoR
	uex[1,idum] = 0.	
	uex[2,idum] = pR	#1./gamma	#pR

	final=np.vstack((x,uex))
	np.savetxt('exact.dat', final.T, delimiter='\t')

#	f.close()	
	return uex
	
import numpy as np
from math import sqrt
x=np.linspace(0., 4., 5000)
x0=2.
t=1.
gamma=1.4
pL=1.
rhoL=1.
pR=0.1 
rhoR=0.125
aL=sqrt(gamma*pL/rhoL)
aR=sqrt(gamma*pR/rhoR)
dum1=2./(gamma+1.)
dum2=(gamma-1.)/(gamma+1.)
dum3=(gamma-1.)/2.
exacta(x,x0,t, gamma, aL, pL, rhoL, aR, pR, rhoR, dum1, dum2, dum3)
