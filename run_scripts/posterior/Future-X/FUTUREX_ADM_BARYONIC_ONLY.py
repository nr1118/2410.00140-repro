#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import neost
from neost.eos import polytropes
from neost.Prior import Prior
from neost.Star import Star
from neost.Likelihood import Likelihood
from scipy.stats import multivariate_normal, gaussian_kde
from neost import PosteriorAnalysis
import numpy
import matplotlib
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot
from pymultinest.solve import solve
import time
import os
import corner as corner


# In[ ]:


import neost.global_imports as global_imports

c = global_imports._c
G = global_imports._G
Msun = global_imports._M_s
pi = global_imports._pi
rho_ns = global_imports._rhons


eos_name = 'polytropes'


EOS = polytropes.PolytropicEoS(crust = 'ceft-Hebeler', rho_t = 2e14, adm_type = 'None')
#EOS.update({'gamma1':2.3, 'gamma2':4., 'gamma3':2.6, 'rho_t1':1.8, 'rho_t2':4.,'mchi':15000, 'gchi_over_mphi':pow(10,-2), 'adm_fraction':1.5, 'ceft':2.6}, max_edsc=True)


# EOS.plot()
# EOS.plot_massradius()
# Here we implement synthetic data based on the future-x scenario with more sources


muM = 2.232   #STROBE-X SOURCE
muR = 10.999
sigM = muM*0.02   #2% uncertainty
sigR = muR*0.02 # 2 % uncertainty in radius
test = multivariate_normal(mean=[muM, muR], cov=[[sigM, 0.0], [0.0, sigR]])


muM2 = 2.001
muR2 = 11.418
sigM2 = muM2*0.02     #2% uncertainty
sigR2 = muR2*0.02   #2 % uncertainty in radius
test2 = multivariate_normal(mean=[muM2, muR2], cov=[[sigM2, 0.0], [0.0, sigR2]])


muM3 = 1.886   
muR3 = 11.474
sigM3 = muM3*0.02   #2% uncertainty in mass
sigR3 = muR3*0.02 # 2 % uncertainty in radius 
test3 = multivariate_normal(mean=[muM3, muR3], cov=[[sigM3, 0.0], [0.0, sigR3]])


muM4 = 1.636  
muR4 = 11.524
sigM4 = muM4*0.02 # 2 % uncerainty in mass 
sigR4 = muR4*0.02  # 2 % uncertainty in radius 
test4 = multivariate_normal(mean=[muM4, muR4], cov=[[sigM4, 0.0], [0.0, sigR4]])


muM5 = 1.457   
muR5 = 11.522
sigM5 = muM5*0.02  # 2% uncertainty
sigR5 = muR5*0.02  # 2 % uncertainty in radius 
test5 = multivariate_normal(mean=[muM5, muR5], cov=[[sigM5, 0.0], [0.0, sigR5]])


muM6 = 1.263   
muR6 = 11.499
sigM6 = muM6*0.02   #2% uncertianty
sigR6 = muR6*0.02  # 2 % uncertainty in radius 
test6 = multivariate_normal(mean=[muM6, muR6], cov=[[sigM6, 0.0], [0.0, sigR6]])


likelihood_functions = [test.pdf,test2.pdf,test3.pdf,test4.pdf,test5.pdf,test6.pdf]
likelihood_params = [['Mass', 'Radius'],['Mass','Radius'],['Mass','Radius'],['Mass','Radius'],['Mass','Radius'],['Mass','Radius']]

# This is not a GW event so we set chirp mass to None
chirp_mass = [None,None,None,None,None,None]
number_stars = len(chirp_mass)

run_name = "FUTUREX_ADM_BARYONIC_ONLY_"

#posterior on mchi is influences such that we don't same mchi < 100 (10^2) see prior corner plots
variable_params = {'gamma1':[1., 4.5], 'gamma2':[0., 8.], 'gamma3':[0.5, 8.], 'rho_t1':[1.5, 8.3], 'rho_t2':[1.5, 8.3],'ceft':[EOS.min_norm, EOS.max_norm]}
# variable_params.update({'ceft':[EOS.min_norm, EOS.max_norm]})
#variable_params = {'mchi':[-2, 9],'gchi_over_mphi': [-5,3],'adm_fraction':[0., 1.7]}
for i in range(number_stars):
	variable_params.update({'rhoc_' + str(i+1):[14.6, 16]})



#static_params = {'gamma1': 2.3, 'gamma2': 4., 'gamma3': 2.6, 'rho_t1': 1.8, 'rho_t2': 4., 'ceft': 2.6}
static_params = {}
# In[ ]:


prior = Prior(EOS, variable_params, static_params, chirp_mass)
likelihood = Likelihood(prior, likelihood_functions, likelihood_params, chirp_mass)

print("Bounds of prior are")
print(variable_params)
print("with model "+eos_name)
print("number of parameters is %d" %len(variable_params))

## TESTING ##
print("Testing prior and likelihood")
cube = numpy.random.rand(50, len(variable_params))
for i in range(len(cube)):
    par = prior.inverse_sample(cube[i])
    print(likelihood.call(par))
print("Testing done")


# In[ ]:


start = time.time()
result = solve(LogLikelihood=likelihood.call, Prior=prior.inverse_sample, n_live_points=2000, evidence_tolerance=0.1,
               n_dims=len(variable_params), sampling_efficiency=0.8, outputfiles_basename=run_name, verbose=True)
end = time.time()
print(end - start)

ewposterior = numpy.loadtxt(run_name +'post_equal_weights.dat')
scattered = []

for i in range(0, len(ewposterior), 1):
    pr = ewposterior[i][0:len(variable_params)]
    par = {e:pr[j] for j, e in enumerate(list(variable_params.keys()))}
    par.update(static_params)
    print(i)
    EOS.update(par, max_edsc=True)

    x = numpy.random.random()
    rhoc_EOS = x*(numpy.log10(EOS.max_edsc)-14.6)+14.6    
     
    rhocpar = 10**rhoc_EOS

    star = Star(rhocpar)
    # star.solve_structure(tabulated_example.energydensities, tabulated_example.pressures)
    star.solve_structure(EOS.energydensities, EOS.pressures) 
                             #ADM_fraction=par['adm_fraction'], mchi=par['mchi'], gchi_over_mphi=par['gchi_over_mphi'])
        
    scattered.append([rhocpar, EOS.eos(rhocpar), star.Mgrav, star.Req,star.tidal]) ## if Rdm_halo == 999 then we have a halo if Rdm_halo ==0 then we have a core



       
 

scattered = numpy.array(scattered)
    


#print(scattered)
numpy.savetxt(run_name + 'scattered.txt', scattered)

scattered = []

for i in range(0, len(ewposterior), 1):
    pr = ewposterior[i][0:len(variable_params)]
    par = {e:pr[j] for j, e in enumerate(list(variable_params.keys()))}
    par.update(static_params)
    print(i)
    EOS.update(par, max_edsc=True)

    tmp = []
    
    for j in range(number_stars):
        rhocpar = 10**par['rhoc_' + str(j+1)] 

        star = Star(rhocpar)
    # star.solve_structure(tabulated_example.energydensities, tabulated_example.pressures)
        star.solve_structure(EOS.energydensities, EOS.pressures) 
                             #ADM_fraction=par['adm_fraction'], mchi=par['mchi'], gchi_over_mphi=par['gchi_over_mphi'])
        
        tmp.append([rhocpar, EOS.eos(rhocpar), star.Mgrav, star.Req,star.tidal]) ## if Rdm_halo == 999 then we have a halo if Rdm_halo ==0 then we have a core
        
    scattered.append(tmp[0])



       
scattered = numpy.array(scattered)
numpy.savetxt(run_name + 'rhoc_scattered.txt', scattered)
PosteriorAnalysis.compute_auxiliary_data_baryonic(run_name, EOS,
                                         variable_params, static_params, chirp_mass)

print("Testing done")




