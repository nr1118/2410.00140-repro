#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import neost
from neost.eos import polytropes
from neost.Star import Star
import numpy
import matplotlib as plt
from matplotlib import pyplot
import time
import pandas as pd
import pathlib


import neost.global_imports as global_imports

c = global_imports._c
G = global_imports._G
Msun = global_imports._M_s
pi = global_imports._pi
rho_ns = global_imports._rhons


# In[ ]:
run_name = 'Percent_diff_intermediate_newfchistep'

def rel_diff(x,y): #really this is the relative percent difference, where we are finding the relative change to x
    return (y-x)/x*100


# In[ ]:


def Radial_diff(eos1,eos2,mchi,fchi,num_stars):
    
    #Defining EoS 0 and 4 parameters for some fchi
    eos1.update({'gamma1':2.3, 'gamma2':4., 'gamma3':2.6, 'rho_t1':1.8, 'rho_t2':4,'mchi':mchi,'gchi_over_mphi':0, 'adm_fraction':fchi, 'ceft':2.6}, max_edsc=True)
    eos2.update({'gamma1':2.3, 'gamma2':4., 'gamma3':2.6, 'rho_t1':1.8, 'rho_t2':4,'mchi':mchi,'gchi_over_mphi':pow(10,-5), 'adm_fraction':fchi, 'ceft':2.6}, max_edsc=True)


    #We define a range of central densities for the baryonic NS and loop over each
    #We can only make one central_densities array for both EoS since they have the same baryonic EoS
    #14.87 was chosen because for the baryonic EoS, that yields ~1 M neutron star.

    central_densities = numpy.logspace(14.87, numpy.log10(eos1.max_edsc), num_stars)
    
    
    relcent_diff = []

    #Including a weight with each value will ensure that all halos have no influence on the average, thus they are effectivley dropped in our calculation!
    Weights = []

    for i,e in enumerate(central_densities):
        eps_centdm_0 = eos1.find_epsdm_cent(eos1.adm_fraction,e)
        Mass_0,Radius_0 = eos1.Mass_Radius(e,eps_centdm_0)

        eps_centdm_5 = eos2.find_epsdm_cent(eos2.adm_fraction,e)
        Mass_5,Radius_5 = eos2.Mass_Radius(e,eps_centdm_5)


        if Radius_0!=0 and Radius_5!=0: #(Testing for an ADM core because if there is an ADM halo then R_star = 0)
            if Mass_0 >=1. and Mass_5 >=1.:
                if numpy.round(Mass_0,2)==numpy.round(Mass_5,2) or numpy.round(Mass_0,3)==numpy.round(Mass_5,3) or numpy.round(Mass_0,4)==numpy.round(Mass_5,4):
                    #relative percent difference between 10^-5 and 0. Where we find the change relative to 0.
                    relcent_diff.append(abs(rel_diff(Radius_0,Radius_5)))
                    Weights.append(1)

                else: 
                    print('The masses do not match')
                    print(Mass_0,Mass_5)
                    print(abs(rel_diff(Mass_0,Mass_5)))

        else: 
            relcent_diff.append(999) #Flag for the halos
            Weights.append(0)
                  
    avg_relcent_diff = numpy.average(relcent_diff,weights = Weights)

    return avg_relcent_diff
    
    

#EOS_0 corresponds the EOS with gchi_over_mphi = 0
EOS_0 = polytropes.PolytropicEoS(crust='ceft-Hebeler', rho_t=2e14,adm_type = 'Fermionic')

#EOS_5 corresponds to the EOS with gchi_over_mphi = 10^-5
EOS_5 = polytropes.PolytropicEoS(crust='ceft-Hebeler', rho_t=2e14,adm_type = 'Fermionic')


# In[ ]:


start = time.time()

mchi_start = 500
mchi_end = 4500

mchi_step = 250
mchi_num_steps = (mchi_end-mchi_start)/mchi_step
mchi_num_steps = int(mchi_num_steps)+1
mchi_array = numpy.linspace(mchi_start,mchi_end,mchi_num_steps)


fchi_start = 0.
fchi_end = 3.

#Array of mass-fractions in which we increment by 0.1
step_size = 0.05
fchi_num_steps = (fchi_end-fchi_start)/step_size
fchi_num_steps = int(fchi_num_steps)+1
#Generating the ADM mass-fraction array of evenly spaces mass-fractions seperated by 0.1%
fchi_array = numpy.linspace(fchi_start,fchi_end,fchi_num_steps)


# In[ ]:


num_stars = 50

Array = numpy.zeros((len(mchi_array),len(fchi_array)))      
   
for i,mchi in enumerate(mchi_array):
    for j,fchi in enumerate(fchi_array):
        avg_reldiff = Radial_diff(EOS_0,EOS_5,mchi,fchi,num_stars)
        Array[i,j] = avg_reldiff
        
           


# In[ ]:



end = time.time()
print(numpy.shape(Array))
print(Array)
data = pd.DataFrame(Array.flatten())

results_directory = '../../../repro/{run_name}/'

pathlib.Path(results_directory).mkdir(parents=True, exist_ok=True) # Create the directory if it doesn't exist

numpy.save(results_directory + 'Relcent_diff_intermediately_stiff_baryonic_newfchistep.npy',Array)
print('The average of the relative percent differences is = '+str(numpy.average(Array.flatten())))
print('The max relative percent difference is = '+str(max(Array.flatten())), 'The min relative percent difference is = '+str(min(Array.flatten())))
print('The mode of the relative percent differences is = '+str(data.mode()))
print('The median of the relative percent differences is = ' +str(numpy.median(Array.flatten())))
print("Execution time of the MR is: " + str(end-start))


# In[ ]:
plots_directory = '../../../plots/'

plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.family'] ='serif'

pyplot.rcParams['xtick.direction'] = 'in'
pyplot.rcParams['xtick.minor.visible'] = True
pyplot.rcParams['ytick.direction'] = 'in'
pyplot.rcParams['ytick.minor.visible'] = True
pyplot.rcParams['xtick.major.size'] = 5
pyplot.rcParams['ytick.major.size'] = 5
pyplot.rcParams['ytick.right'] = True
pyplot.rcParams['xtick.top'] = True   

fig, ax = pyplot.subplots(1,1, figsize=(14,5))
X, Y = numpy.meshgrid(mchi_array, fchi_array)
pos1 = ax.pcolormesh(X, Y, Array.T,shading ='gouraud',cmap='viridis') #shading ='gouraud' used for pcolormesh
pos1.set_clim(0.0,.004)
cbar = fig.colorbar(pos1, ax=ax,ticks = [0.,.001,.002,.003,.004])
cbar.set_label('RRPD [$\%$]',rotation = 90,labelpad = 25,fontsize=24)
cbar.ax.tick_params(labelsize=22)
ax.set_ylabel('F$_\chi \, [\%]$',fontsize=28)
ax.set_xlabel('m$_\chi \, [\mathrm{MeV}]$',fontsize=28)
ax.tick_params(axis='both', which='major', labelsize=22)
ax.set_title('Intermediately Stiff EoS with $\Delta F_\chi = 0.05 \%$',fontsize = 24)
pyplot.savefig(plots_directory + 'Percent_diff_intermediate_stiff_plot_newfchistep.png')
pyplot.show()


