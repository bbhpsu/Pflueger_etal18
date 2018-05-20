from functions import *
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import yaml
stream = open("input.dat","r").read()
doc = yaml.load(stream)

np.seterr(divide='ignore')


# Probability of a binary residing at separation a: \rho_a
def rhoa(a,q,m,mdot):
	if(ttot(a,q,m,mdot)/(365.25*24*60*60*1E+9) > 13.6):
		rhoa=0
	else:
		if(tevolv(q,m,mdot)/(365.25*24*60*60*1E+9) > 13.6):
			rhoa=-1.0/aregimes(a,q,m,mdot)/(13.6*365.25*24*60*60*1E+9)
		else:	
			rhoa=-1.0/aregimes(a,q,m,mdot)/tevolv(q,m,mdot)
	return rhoa


# Mass ratio probability distribution: \rho_q
def rhoq(q):
	return 1.0/(doc['xsia'])*q**(-0.3)*(1-q)*np.exp(-doc['beta']/(q**2.0))


# Probability of simultaneous measurement of V and \Delta V: P_{V,\Delta V}
def extrinsic(a,q,m):
	if ((np.fabs(c1(a,q))>1) or (np.fabs(c2(a,q,m))>1)):
		return 0.0
	else:
		lower=min(np.arcsin(c1(a,q)),np.arccos(c2(a,q,m)))
		upper=min(np.arccos(c2(a,q,m)),np.pi/2.0)
		return quad(lambda phi: 2.0/np.pi*(1.0-cA(phi,a,q)/2.0-cB(phi,a,q,m)/2.0-np.fabs(cB(phi,a,q,m)/2.0-cA(phi,a,q)/2.0)),lower,upper)[0]


#Likelihood of Detection: \rho_a * \rho_q * P_{V,\Delta V}
def likelihood(a,q,m,mdot):
	return rhoa(a,q,m,mdot)*rhoq(q)*extrinsic(a,q,m)



# 2D plot of the Likelihood function (q vs. a):

qi = np.linspace(doc['qmin'],doc['qmax'],doc['res'][0])
ai = np.logspace(doc['amin'],doc['amax'],doc['res'][1])

zi = np.empty((qi.size,ai.size))
for y in range(qi.size):
	for x in range(ai.size):
		zi[y,x]=np.log10(likelihood(ai[x],qi[y],doc['mlist'][0],doc['mdotlist'][0])*rg(doc['mlist'][0])) 

#Setting contour levels in both plot and color bar
v=[-9,-8,-7.5,-7.2,-6.95]
norml = colors.BoundaryNorm(boundaries=v, ncolors=256)

cmap=plt.cm.summer
cmap.set_under('white')
cmap.set_over('yellow')

plt.contour(np.log10(ai), qi, zi, v, linewidths=0.5, colors='k')
plt.contourf(np.log10(ai), qi, zi, levels=v, extend='max', cmap=cmap, norm=norml)
x = plt.colorbar(ticks=v)
plt.xlabel('log($a/r_g$)')
plt.ylabel('q')
plt.title('Likelihood of detection: log($\mathcal{L}r_g$)')
plt.savefig('Likelihood_q_a.pdf')












