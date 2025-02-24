#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib as plt
from matplotlib import pyplot
import numpy as np
from scipy.interpolate import UnivariateSpline


# In[2]:


plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('text', usetex=True)

pyplot.rc('text', usetex=True)
pyplot.rc('font', family='serif')

pyplot.rcParams['xtick.direction'] = 'in'
pyplot.rcParams['xtick.minor.visible'] = True
pyplot.rcParams['ytick.direction'] = 'in'
pyplot.rcParams['ytick.minor.visible'] = True
pyplot.rcParams['xtick.major.size'] = 5
pyplot.rcParams['ytick.major.size'] = 5
pyplot.rcParams['ytick.right'] = True
pyplot.rcParams['xtick.top'] = True


# In[3]:


from uncertainties import ufloat
from uncertainties.umath import *  # sin(), etc. 
from astropy import units as u
from astropy.coordinates import SkyCoord 


# In[4]:


vs = ufloat(430.4,5.4) #km


# In[5]:


from astropy.coordinates import Galactocentric
from astropy.coordinates import galactocentric_frame_defaults
import astropy.coordinates as coord
_ = galactocentric_frame_defaults.set('latest')
Galactocentric()


# In[6]:


dchi_c = ufloat(8.1,0.7) #kpc
rhochi_c = ufloat(5.22,.46)*1e-2 #solar mass per pc^3
rhochi = lambda d: rhochi_c*(dchi_c/d)* (1+d/dchi_c)**(-2)


rsc_disk = ufloat(3,0.0) #kpc
rhoc_disk = ufloat(15,0.0) #solar mass/ pc^3

rho_disk = lambda d: rhoc_disk * exp(-d/rsc_disk)


rho_total = lambda d: rho_disk(d) + rhochi(d)


# In[7]:


#defining the mass-fraction a ratio of the dark matter density over the sum of the densities
f_chi = lambda d: rhochi(d)/(rho_total(d))


# In[8]:


#J0740

#Sgr A* in the icrs frame,i.e, frame of the barycenter of the solar system, which is basically just outside of the sun
sgrA = coord.SkyCoord(ra =266.41681663*u.degree, dec =-29.00782497*u.degree, distance = 8.3*u.kpc,frame = 'icrs')



#PSR J0740+6620 see Fonseca et al. 2021 APJL 915 L12 for distance measurement
J0740 = coord.SkyCoord(ra =115.19082917*u.degree, dec = 66.34266667*u.degree, distance = 1.14*u.kpc,frame = 'icrs')


dist_0740toGC = sgrA.separation_3d(J0740).value


print(dist_0740toGC)

print('Fchi Estimate: '+str(f_chi(dist_0740toGC)*100))



# In[9]:


#J0437
psr_0437 =coord.SkyCoord(ra =69.31583108*u.degree, dec = -47.25237343*u.degree, distance = .157*u.kpc)


dist_0437toGC = sgrA.separation_3d(psr_0437).value

print(dist_0437toGC)

print('Fchi Estimate: '+str(f_chi(dist_0437toGC)*100))


# In[10]:


#J0030
psr_0030 =coord.SkyCoord(ra =7.61428208*u.degree, dec = 4.86103056*u.degree, distance = 0.325*u.kpc)


dist_0030toGC = sgrA.separation_3d(psr_0030).value

print(dist_0030toGC)

print('Fchi Estimate: '+str(f_chi(dist_0030toGC)*100))


# In[ ]:




