#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import neost
from neost.eos import polytropes
from neost.Star import Star
import numpy
import matplotlib as plt
from matplotlib import pyplot
import pathlib


import neost.global_imports as global_imports

c = global_imports._c
G = global_imports._G
Msun = global_imports._M_s
pi = global_imports._pi
rho_ns = global_imports._rhons




gamma1 = numpy.array([1.1,2.3,3.8])
gamma2 = numpy.array([5.7,4.,3.8])
gamma3 = numpy.array([3.,2.6,2.2])
rho_t1 = numpy.array([2,1.8,1.8])
rho_t2 = numpy.array([3.7,4.,2.7])




lns = []
fig, ax = pyplot.subplots(1,1, figsize=(12, 10))
for i in range(len(gamma1)):
    EOS = polytropes.PolytropicEoS(crust='ceft-Hebeler', rho_t=2e14,adm_type = 'None')
    EOS.update({'gamma1':gamma1[i], 'gamma2':gamma2[i], 'gamma3':gamma3[i], 'rho_t1':rho_t1[i], 'rho_t2':rho_t2[i], 'ceft':2.6}, max_edsc=True)
    central_densities = numpy.logspace(14.3,numpy.log10(EOS.max_edsc),50)
    MR = numpy.zeros((len(central_densities), 2))
    for k, eps in enumerate(central_densities):
        star = Star(eps)
        star.solve_structure(EOS.energydensities, EOS.pressures)
        MR[k] = star.Mrot, star.Req
    if (i==0):
        lns1 = ax.plot(MR[:,1], MR[:,0],label = 'Soft EoS', lw=2.5)
    elif(i==1):
        lns1 = ax.plot(MR[:,1], MR[:,0],label = 'Int. Stiff EoS', lw=2.5)
    else:
        lns1 = ax.plot(MR[:,1], MR[:,0],label ='Stiff EoS', lw=2.5)
        
    lns = lns + lns1

    


results_directory = '../../../repro/{run_name}/' 

pathlib.Path(results_directory).mkdir(parents=True, exist_ok=True) # Create the directory if it doesn't exist

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
    
ax.set_xlabel(r'Radius [km]', fontsize=28)
ax.set_ylabel(r'Mass [M$_\odot$]', fontsize=28)
ax.tick_params(width=2, labelsize=22, direction='in')
labs = [l.get_label() for l in lns]
ax.legend(lns,labs,loc='best', fontsize=26)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,fontsize=22)

for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)
pyplot.tight_layout()
pyplot.show()
pyplot.savefig(results_directory + 'MR_fixed.png')




