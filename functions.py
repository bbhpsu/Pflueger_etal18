import numpy as np
import yaml
stream = open("input.dat","r").read()
doc = yaml.load(stream)

#Gravitational radius
def rg(m): 	
	return doc['G']*m*doc['Msolar']/(doc['c']**2.0)

#Transistion Distances
#Gas pressure dominated regime to radiation pressure regime
def Rt(q,m,mdot): 
	return 2849.74*doc['alpha']**(2.0/21)*(m/1E+7)**(2.0/21)*(mdot)**(16.0/21)

#Radiation pressure dominated regime to gravitational wave regime
def Rgrav(q,m,mdot): 
	return 0.01125*(q**2.0/((1.0+q)**4.0))*doc['alpha']**(-2.0)*mdot**(-4.0)

#Gas pressure dominated regime to gravitational wave regime
def Rgtogw(q,m,mdot): 
	return 260.15*doc['alpha']**(-8.0/26)*mdot**(-2.0/13)*(m/(1.0E+7))**(1.0/13)*((q**(10.0/26))/((1+q)**(10.0/13)))

#Inspiral rates for each regime
# Inspiral rate in gas pressure dominated portion of the disk
def Agas(a,q,m,mdot):
	return -3.0/2*doc['alpha']*(0.00266*doc['alpha']**(-1.0/10)*mdot**(1.0/5)*(m/1E+7)**(-1.0/10)*(a/3)**(1.0/20))**2*(doc['c']**2/(2*a))**(1.0/2)

# Inspiral rate in radiation pressure dominated portion of the disk
def Arad(a,q,m,mdot):
	return -3.0/2*doc['alpha']*(3.56*(mdot)*(a/3.0)**(-1.0))**2.0*(doc['c']**2/(2*a))**(1.0/2)

# Gravitational wave inspiral rate
def Agrav(a,q,m,mdot):
	return -64.0/5*doc['c']*(q/(1+q)**2)*a**(-3.0)

# Total time spent within each regime:
# Time in gas pressure dominated regime:
def Tg1(q,m,mdot):
	return 5.24852E+6*doc['alpha']**(-4.0/5)*(m/1.0E+7)**(6.0/5)*mdot**(-2.0/5)*(doc['a0']**(7.0/5)-Rt(q,m,mdot)**(7.0/5))

def Tg2(q,m,mdot):
	return 5.24852E+6*doc['alpha']**(-4.0/5)*(m/1.0E+7)**(6.0/5)*mdot**(-2.0/5)*(doc['a0']**(7.0/5)-Rgtogw(q,m,mdot)**(7.0/5))

# Time in radiation pressure dominated regime:
def Trd1(q,m,mdot):
	return 0.116682*doc['alpha']**(-1.0)*(m/1.0E+7)*mdot**(-2.0)*(Rt(q,m,mdot)**(7.0/2)-Rgrav(q,m,mdot)**(7.0/2))

def Trd2(q,m,mdot):
	return 0.116682*doc['alpha']**(-1.0)*(m/1.0E+7)*mdot**(-2.0)*(Rgtogw(q,m,mdot)**(7.0/2)-Rgrav(q,m,mdot)**(7.0/2))

#Time in gravitational wave regime
def Tgw1(q,m,mdot):
	return 0.964988*((1+q)**2.0/q)*(m/1.0E+7)*Rgrav(q,m,mdot)**(4.0)

def Tgw2(q,m,mdot):
	return 0.964988*((1+q)**2.0/q)*(m/1.0E+7)*Rt(q,m,mdot)**(4.0)

def Tgw3(q,m,mdot):
	return 0.964988*((1+q)**2.0/q)*(m/1.0E+7)*Rgtogw(q,m,mdot)**(4.0)

# Total time of evolution through all regimes:
def tevolv(q,m,mdot):	
	if ((Rgrav(q,m,mdot)>doc['a0']) and (doc['a0']>Rt(q,m,mdot)) and (Rt(q,m,mdot)>Rgtogw(q,m,mdot))):
		tevolv = Tg1(q,m,mdot)+0+Tgw2(q,m,mdot)
	if ((Rgrav(q,m,mdot)>doc['a0']) and (doc['a0']>Rgtogw(q,m,mdot)) and (Rgtogw(q,m,mdot)>Rt(q,m,mdot))):
		tevolv = Tg2(q,m,mdot)+0+Tgw3(q,m,mdot)
	if ((doc['a0']>Rt(q,m,mdot)) and (Rt(q,m,mdot)>Rgrav(q,m,mdot)) and (Rgrav(q,m,mdot)>Rgtogw(q,m,mdot))):
		tevolv = Tg1(q,m,mdot)+Trd1(q,m,mdot)+Tgw1(q,m,mdot)
	if ((doc['a0']>Rt(q,m,mdot)) and (Rt(q,m,mdot)>Rgtogw(q,m,mdot)) and (Rgtogw(q,m,mdot)>Rgrav(q,m,mdot))):
		tevolv = Tg1(q,m,mdot)+Trd1(q,m,mdot)+Tgw1(q,m,mdot)
	if ((doc['a0']>Rgrav(q,m,mdot)) and (Rgrav(q,m,mdot)>Rt(q,m,mdot)) and (Rt(q,m,mdot)>Rgtogw(q,m,mdot))):
		tevolv = Tg1(q,m,mdot)+0+Tgw2(q,m,mdot)
	if ((doc['a0']>Rgrav(q,m,mdot)) and (Rgrav(q,m,mdot)>Rgtogw(q,m,mdot)) and (Rgtogw(q,m,mdot)>Rt(q,m,mdot))):
		tevolv = Tg2(q,m,mdot)+0+Tgw3(q,m,mdot) 
	if ((doc['a0']>Rgtogw(q,m,mdot)) and (Rgtogw(q,m,mdot)>Rt(q,m,mdot)) and (Rt(q,m,mdot)>Rgrav(q,m,mdot))):
		tevolv = Tg2(q,m,mdot)+Trd2(q,m,mdot)+Tgw1(q,m,mdot)
	if ((doc['a0']>Rgtogw(q,m,mdot)) and (Rgtogw(q,m,mdot)>Rgrav(q,m,mdot)) and (Rgrav(q,m,mdot)>Rt(q,m,mdot))):
		tevolv = Tg2(q,m,mdot)+Trd2(q,m,mdot)+Tgw1(q,m,mdot)
	return tevolv

# Time it takes binary to evolve to a certain orbital separation:
def Tg1a(a,q,m,mdot):
	return 5.24852E+6*doc['alpha']**(-4.0/5)*(m/1.0E+7)**(6.0/5)*mdot**(-2.0/5)*(doc['a0']**(7.0/5)-a**(7.0/5))

def Trd1a(a,q,m,mdot):
	return 0.116682*doc['alpha']**(-1.0)*(m/1.0E+7)*mdot**(-2.0)*(Rt(q,m,mdot)**(7.0/2)-a**(7.0/2))

def Trd2a(a,q,m,mdot):
	return 0.116682*doc['alpha']**(-1.0)*(m/1.0E+7)*mdot**(-2.0)*(Rgtogw(q,m,mdot)**(7.0/2)-a**(7.0/2))

def Tgw1a(a,q,m,mdot):
	return 0.964988*((1+q)**2.0/q)*(m/1.0E+7)*(Rgrav(q,m,mdot)-a)**(4.0)

def Tgw2a(a,q,m,mdot):
	return 0.964988*((1+q)**2.0/q)*(m/1.0E+7)*(Rt(q,m,mdot)-a)**(4.0)

def Tgw3a(a,q,m,mdot):
	return 0.964988*((1+q)**2.0/q)*(m/1.0E+7)*(Rgtogw(q,m,mdot)-a)**(4.0)

def ttot(a,q,m,mdot):
	t=0
	if((Rgrav(q,m,mdot)>doc['a0'])and(doc['a0']>Rt(q,m,mdot))and(Rt(q,m,mdot)>Rgtogw(q,m,mdot))):
		if(a>Rt(q,m,mdot)):
			t=t+Tg1a(a,q,m,mdot)
		else:
			t=t+Tg1(q,m,mdot)+Tgw2a(a,q,m,mdot)
	if((Rgrav(q,m,mdot)>doc['a0'])and(doc['a0']>Rgtogw(q,m,mdot))and(Rgtogw(q,m,mdot)>Rt(q,m,mdot))):
		if(a>Rgtogw(q,m,mdot)):
			t=t+Tg1a(a,q,m,mdot)
		else:
			t=t+Tg2(q,m,mdot)+Tgw3a(a,q,m,mdot)
	if((doc['a0']>Rgrav(q,m,mdot))and(Rgrav(q,m,mdot)>Rt(q,m,mdot))and(Rt(q,m,mdot)>Rgtogw(q,m,mdot))):
		if(a>Rt(q,m,mdot)):
			t=t+Tg1a(a,q,m,mdot)
		else:
			t=t+Tg1(q,m,mdot)+Tgw2a(a,q,m,mdot)
	if((doc['a0']>Rgrav(q,m,mdot))and(Rgrav(q,m,mdot)>Rgtogw(q,m,mdot))and(Rgtogw(q,m,mdot)>Rt(q,m,mdot))):
		if(a>Rgtogw(q,m,mdot)):
			t=t+Tg1a(a,q,m,mdot)
		else:
			t=t+Tg2(q,m,mdot)+Tgw3a(a,q,m,mdot)
	if((doc['a0']>Rt(q,m,mdot))and(Rt(q,m,mdot)>Rgrav(q,m,mdot))and(Rgrav(q,m,mdot)>Rgtogw(q,m,mdot))):
		if(a>Rt(q,m,mdot)):
			t=t+Tg1a(a,q,m,mdot)
		if((Rt(q,m,mdot)>a)and(a>Rgrav(q,m,mdot))):
			t=t+Tg1(q,m,mdot)+Trd1a(a,q,m,mdot)
		if(Rgrav(q,m,mdot)>a):
			t=t+Tg1(q,m,mdot)+Trd1(q,m,mdot)+Tgw1a(a,q,m,mdot)
	if((doc['a0']>Rt(q,m,mdot))and(Rt(q,m,mdot)>Rgtogw(q,m,mdot))and(Rgtogw(q,m,mdot)>Rgrav(q,m,mdot))):
		if(a>Rt(q,m,mdot)):
			t=t+Tg1a(a,q,m,mdot)
		if((Rt(q,m,mdot)>a)and(a>Rgrav(q,m,mdot))):
			t=t+Tg1(q,m,mdot)+Trd1a(a,q,m,mdot)
		if(Rgrav(q,m,mdot)>a):
			t=t+Tg1(q,m,mdot)+Trd1(q,m,mdot)+Tgw1a(a,q,m,mdot)
	if((doc['a0']>Rgtogw(q,m,mdot))and(Rgtogw(q,m,mdot)>Rt(q,m,mdot))and(Rt(q,m,mdot)>Rgrav(q,m,mdot))):
		if(a>Rgtogw(q,m,mdot)):
			t=t+Tg1a(a,q,m,mdot)
		if((Rgtogw(q,m,mdot)>a)and(a>Rgrav(q,m,mdot))):
			t=t+Tg2(q,m,mdot)+Trd2a(a,q,m,mdot)
		if(Rgrav(q,m,mdot)>a):
			t=t+Tg2(q,m,mdot)+Trd2(q,m,mdot)+Tgw1a(a,q,m,mdot)
	if((doc['a0']>Rgtogw(q,m,mdot))and(Rgtogw(q,m,mdot)>Rgrav(q,m,mdot))and(Rgrav(q,m,mdot)>Rt(q,m,mdot))):
		if(a>Rgtogw(q,m,mdot)):
			t=t+Tg1a(a,q,m,mdot)
		if((Rgtogw(q,m,mdot)>a)and(a>Rgrav(q,m,mdot))):
			t=t+Tg2(q,m,mdot)+Trd2a(a,q,m,mdot)
		if(Rgrav(q,m,mdot)>a):
			t=t+Tg2(q,m,mdot)+Trd2(q,m,mdot)+Tgw1a(a,q,m,mdot)
	return t

# Inspiral rate for the 8 cases and where the regimes switch:
def aregimes(a,q,m,mdot):
	#r=0
	if((Rgrav(q,m,mdot)>doc['a0'])and(doc['a0']>Rt(q,m,mdot)and(Rt(q,m,mdot)>Rgtogw(q,m,mdot)))):
		if((doc['a0']>=a)and(a>Rt(q,m,mdot))):
			aregimes=Agas(a,q,m,mdot)
		if((Rt(q,m,mdot)>=a)and(a>0)):
			aregimes=Agrav(a,q,m,mdot)
	if((Rgrav(q,m,mdot)>doc['a0'])and(doc['a0']>Rgtogw(q,m,mdot)and(Rgtogw(q,m,mdot)>Rt(q,m,mdot)))):
		if((doc['a0']>=a)and(a>Rgtogw(q,m,mdot))):
			aregimes=Agas(a,q,m,mdot)
		if((Rgtogw(q,m,mdot)>=a)and(a>0)):
			aregimes=Agrav(a,q,m,mdot)
	if((doc['a0']>Rt(q,m,mdot))and(Rt(q,m,mdot)>Rgrav(q,m,mdot)and(Rgrav(q,m,mdot)>Rgtogw(q,m,mdot)))):
		if((doc['a0']>=a)and(a>Rt(q,m,mdot))):
			aregimes=Agas(a,q,m,mdot)
		if((Rt(q,m,mdot)>=a)and(a>Rgrav(q,m,mdot))):
			aregimes=Arad(a,q,m,mdot)
		if((Rgrav(q,m,mdot)>=a)and(a>0)):
			aregimes=Agrav(a,q,m,mdot)	
	if((doc['a0']>Rt(q,m,mdot))and(Rt(q,m,mdot)>Rgtogw(q,m,mdot)and(Rgtogw(q,m,mdot)>Rgrav(q,m,mdot)))):
		if((doc['a0']>=a)and(a>Rt(q,m,mdot))):
			aregimes=Agas(a,q,m,mdot)
		if((Rt(q,m,mdot)>=a)and(a>Rgrav(q,m,mdot))):
			aregimes=Arad(a,q,m,mdot)
		if((Rgrav(q,m,mdot)>=a)and(a>0)):
			aregimes=Agrav(a,q,m,mdot)
	if((doc['a0']>Rgrav(q,m,mdot))and(Rgrav(q,m,mdot)>Rt(q,m,mdot)and(Rt(q,m,mdot)>Rgtogw(q,m,mdot)))):
		if((doc['a0']>=a)and(a>Rt(q,m,mdot))):
			aregimes=Agas(a,q,m,mdot)
		if((Rt(q,m,mdot)>=a)and(a>0)):
			aregimes=Agrav(a,q,m,mdot)
	if((doc['a0']>Rgrav(q,m,mdot))and(Rgrav(q,m,mdot)>Rgtogw(q,m,mdot)and(Rgtogw(q,m,mdot)>Rt(q,m,mdot)))):
		if((doc['a0']>=a)and(a>Rgtogw(q,m,mdot))):
			aregimes=Agas(a,q,m,mdot)
		if((Rgtogw(q,m,mdot)>=a)and(a>0)):
			aregimes=Agrav(a,q,m,mdot)	
	if((doc['a0']>Rgtogw(q,m,mdot))and(Rgtogw(q,m,mdot)>Rt(q,m,mdot)and(Rt(q,m,mdot)>Rgrav(q,m,mdot)))):
		if((doc['a0']>=a)and(a>Rgtogw(q,m,mdot))):
			aregimes=Agas(a,q,m,mdot)
		if((Rgtogw(q,m,mdot)>=a)and(a>Rgrav(q,m,mdot))):
			aregimes=Arad(a,q,m,mdot)
		if((Rgrav(q,m,mdot)>=a)and(a>0)):
			aregimes=Agrav(a,q,m,mdot)
	if((doc['a0']>Rgtogw(q,m,mdot))and(Rgtogw(q,m,mdot)>Rgrav(q,m,mdot)and(Rgrav(q,m,mdot)>Rt(q,m,mdot)))):
		if((doc['a0']>=a)and(a>Rgtogw(q,m,mdot))):
			aregimes=Agas(a,q,m,mdot)
		if((Rgtogw(q,m,mdot)>=a)and(a>Rgrav(q,m,mdot))):
			aregimes=Arad(a,q,m,mdot)
		if((Rgrav(q,m,mdot)>=a)and(a>0)):
			aregimes=Agrav(a,q,m,mdot)
	return aregimes


# Functions for calculation of P_{V,\Delta V}
def c1(a,q):
	return (1.0+q)*np.sqrt(a)*doc['vlim']/doc['c']

def c2(a,q,m):
	return ((1.0+q)*np.sqrt(a))/(np.sin(doc['c']**3.0/doc['G']*doc['deltat']/(m*doc['Msolar']*(a**(3.0/2)))))*doc['deltavlim']/doc['c']

def cA(phi,a,q):
	return c1(a,q)/np.sin(phi)

def cB(phi,a,q,m):
	return c2(a,q,m)/np.cos(phi)

def f(a,q,m):
	return 2.0/np.pi*(1.0-cA(phi,a,q)/2.0-cB(phi,a,q,m)/2-np.fabs(cB(phi,a,q,m)/2.0-cA(phi,a,q)/2.0))

# Numerical integration with numpy:
def integrand(x,a,q,m):
    return f(x,a,q,m)













		
	








