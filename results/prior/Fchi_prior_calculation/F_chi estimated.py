#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib as plt
from matplotlib import pyplot
import numpy as np


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


from astropy import units as u
from astropy.coordinates import SkyCoord 
from astropy.coordinates import Galactocentric
from astropy.coordinates import galactocentric_frame_defaults
import astropy.coordinates as coord


# In[4]:


dchi_c_upper = 8.1+0.7 # 8.1 +/- 0.7 #kpc
rhochi_c_upper =5.22e-2 + .46e-2 #5.22e-2,  +/-.46e-2 #solar mass per pc^3
rhochi_upper = lambda d: rhochi_c_upper*(dchi_c_upper/d)* (1+d/dchi_c_upper)**(-2)


rsc_disk = 3   #kpc
rhoc_disk = 15 #solar mass/ pc^3

rho_disk = lambda d: rhoc_disk * np.exp(-d/rsc_disk)


rho_total_upper = lambda d: rho_disk(d) + rhochi_upper(d)


# In[5]:


#defining the mass-fraction a ratio of the dark matter density over the sum of the densities
f_chi = lambda d: rhochi_upper(d)/(rho_total_upper(d))


# In[6]:


#Sgr A* in the icrs frame,i.e, frame of the barycenter of the solar system, which is basically just outside of the sun
sgrA = coord.SkyCoord(ra =266.41681663*u.degree, dec =29.00782497*u.degree, distance = 8.3*u.kpc,frame = 'icrs')


# In[7]:


#J0740
dist_0740toGC = 8.6 #kpc following 1910.09925

print('Distance 0740 to GC:', dist_0740toGC)

print('Fchi Estimate 0740: '+str(f_chi(dist_0740toGC)*100))



#J0437
psr_0437 =coord.SkyCoord(ra =69.31583108*u.degree, dec = -47.25237343*u.degree, distance = .157*u.kpc)


dist_0437toGC = sgrA.separation_3d(psr_0437).value

print('Distance 0437 to GC:', dist_0437toGC)

print('Fchi Estimate 0437: '+str(f_chi(dist_0437toGC)*100))




#J0030
psr_0030 =coord.SkyCoord(ra =7.61428208*u.degree, dec = 4.86103056*u.degree, distance = 0.325*u.kpc)


dist_0030toGC = sgrA.separation_3d(psr_0030).value

print('Distance 0030 to GC:', dist_0030toGC)

print('Fchi Estimate 0030: '+str(f_chi(dist_0030toGC)*100))


# In[8]:


dchi_c_lower = 8.1-0.7 # 8.1 +/- 0.7 #kpc
rhochi_c_lower =(5.22-0.46)*1e-2 #5.22e-2,  +/-.46e-2 #solar mass per pc^3
rhochi_lower = lambda d: rhochi_c_lower*(dchi_c_lower/d)* (1+d/dchi_c_lower)**(-2)


rsc_disk = 3   #kpc
rhoc_disk = 15 #solar mass/ pc^3

rho_disk = lambda d: rhoc_disk * np.exp(-d/rsc_disk)


rho_total_lower = lambda d: rho_disk(d) + rhochi_lower(d)

#defining the mass-fraction a ratio of the dark matter density over the sum of the densities
f_chilower = lambda d: rhochi_lower(d)/(rho_total_lower(d))


# In[9]:


#J0740
dist_0740toGC = 8.6 #kpc

print('Distance 0740 to GC:', dist_0740toGC)

print('Fchi Estimate 0740: '+str(f_chilower(dist_0740toGC)*100))



#J0437
psr_0437 =coord.SkyCoord(ra =69.31583108*u.degree, dec = -47.25237343*u.degree, distance = .157*u.kpc)


dist_0437toGC = sgrA.separation_3d(psr_0437).value

print('Distance 0437 to GC:', dist_0437toGC)

print('Fchi Estimate 0437: '+str(f_chilower(dist_0437toGC)*100))




#J0030
psr_0030 =coord.SkyCoord(ra =7.61428208*u.degree, dec = 4.86103056*u.degree, distance = 0.325*u.kpc)


dist_0030toGC = sgrA.separation_3d(psr_0030).value

print('Distance 0030 to GC:', dist_0030toGC)

print('Fchi Estimate 0030: '+str(f_chilower(dist_0030toGC)*100))


# In[ ]:




