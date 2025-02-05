#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import neost
from neost.eos import polytropes
from neost.Prior import Prior
from neost.Star import Star
from neost.Likelihood import Likelihood
from scipy.stats import multivariate_normal, gaussian_kde
import numpy
import matplotlib
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot
from pymultinest.solve import solve
import time
import os
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import seaborn as sns
import corner as corner
import getdist
from getdist import plots


# In[ ]:


import neost.global_imports as global_imports

c = global_imports._c
G = global_imports._G
Msun = global_imports._M_s
pi = global_imports._pi
rho_ns = global_imports._rhons


eos_name = 'polytropes'


EOS = polytropes.PolytropicEoS(crust = 'ceft-Hebeler', rho_t = 2e14, adm_type = 'None')


# EOS.plot()
# EOS.plot_massradius()
# Here we implement old NICER data on J0740 and J0030 from Riley et al.


# PSR J0740+6620 (arXiv: 2105.06980)
muM = 2.07  
muR = 12.39
sigM = 0.07  #published uncertainty
sigR_plus = 1.30
sigR_minus = 0.98
J0740_mr_samples = numpy.load('Riley00740_mr_posterior_samples.npy').T
J0740_kde = gaussian_kde(J0740_mr_samples)

# PSR J0030+0451(arXiv: 1912.05702)
muM2 = 1.34 
muR2 = 12.71
sigM2 = 0.15 #published uncertainty
sigR2_plus = 1.13
sigR2_minus = 1.18
J0030_mr_samples = numpy.load('Riley0030_mr_posterior_samples.npy').T
J0030_kde = gaussian_kde(J0030_mr_samples)


likelihood_functions = [J0740_kde.pdf,J0030_kde.pdf]
likelihood_params = [['Mass', 'Radius'],['Mass','Radius']]

# This is not a GW event so we set chirp mass to None
chirp_mass = [None,None]
number_stars = len(chirp_mass)

run_name = "POSTERIOR_FERMIONIC_REAL_"

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

start = time.time()
result = solve(LogLikelihood=likelihood.call, Prior=prior.inverse_sample, n_live_points=3000, evidence_tolerance=0.1,
               n_dims=len(variable_params), sampling_efficiency=0.8, outputfiles_basename=run_name, verbose=True)
end = time.time()
print(end - start)
# In[ ]:


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
print("Testing done")

def calc_bands(x, y):
    miny = numpy.zeros((len(y),3))
    maxy = numpy.zeros((len(y),3))
    
    for i in range(len(y)):
        z = y[i][y[i]>0.0]
        if len(z)<200:
            print('sample too small for %.2f' %x[i])
            continue
        kde = gaussian_kde(z)
        testz = numpy.linspace(min(z),max(z), 1000)
        pdf = kde.pdf(testz)
        array = pdf
        index_68 = numpy.where(numpy.cumsum(numpy.sort(array)[::-1]) < sum(array)*0.6827)[0]
        index_68 = numpy.argsort(array)[::-1][index_68]
        index_95 = numpy.where(numpy.cumsum(numpy.sort(array)[::-1]) < sum(array)*0.95)[0]
        index_95 = numpy.argsort(array)[::-1][index_95]
        miny[i] =  x[i], min(testz[index_68]), min(testz[index_95])
        maxy[i] =  x[i], max(testz[index_68]), max(testz[index_95])
        
    miny = miny[~numpy.all(miny == 0, axis=1)]
    maxy = maxy[~numpy.all(maxy == 0, axis=1)]
    return miny, maxy


def compute_auxiliary_data(root_name, EOS, variable_params, static_params, chirp_masses): 
    ewposterior = numpy.loadtxt(root_name + 'post_equal_weights.dat')
    print("Total number of samples is %d" %(len(ewposterior)))

    masses = numpy.linspace(.2, 2.9, 50)
    energydensities = numpy.logspace(14.2, 16, 50)
    
    radii = numpy.zeros((len(masses), len(ewposterior)))
    pressures = numpy.zeros((len(masses), len(ewposterior)))
    minradii = numpy.zeros((3, len(masses)))
    maxradii = numpy.zeros((3, len(masses)))
    minpres = numpy.zeros((3, len(energydensities)))
    maxpres = numpy.zeros((3, len(energydensities)))
    MR_prpr_pp = numpy.zeros((len(ewposterior), 2))



    with alive_bar(len(ewposterior)) as bar:
        for i in range(0, len(ewposterior), 1):

            pr = ewposterior[i][0:len(variable_params)]
            par = {e:pr[j] for j, e in enumerate(list(variable_params.keys()))}
            par.update(static_params)
            EOS.update(par, max_edsc=True)

            rhopres = UnivariateSpline(EOS.massdensities, EOS.pressures, k=1, s=0)
            edsrho = UnivariateSpline(EOS.energydensities, EOS.massdensities, k=1, s=0)
            max_rhoc = edsrho(EOS.max_edsc)

            #pressures_rho[:,i][energydensities<max_rhoc] = rhopres(energydensities[energydensities<max_rhoc])
            pressures[:,i][energydensities<EOS.max_edsc] = EOS.eos(energydensities[energydensities<EOS.max_edsc])
            
            rhocs = numpy.logspace(14.5, numpy.log10(EOS.max_edsc), 40)
            M = numpy.zeros(len(rhocs))
            R = numpy.zeros(len(rhocs))
            for j, eps in enumerate(rhocs):
                star = Star(eps)
                star.solve_structure(EOS.energydensities, EOS.pressures)
                M[j] = star.Mgrav
                R[j] = star.Req

            M, indices = numpy.unique(M, return_index=True)
            MR = UnivariateSpline(M, R[indices], k=1, s=0, ext=1)
            rhocM = UnivariateSpline(M, rhocs[indices], k=1, s=0)
            
            # ADM is inlcuded here as this is the scattered portion of the run code
            #rhocpar = numpy.array([10**v for k,v in par.items() if 'rhoc' in k])
            #tmp = []
            #for j, e in enumerate(rhocpar):
             #   star = Star(e)
             #   star.solve_structure(EOS.energydensities, EOS.pressures)
             #   tmp.append([e, EOS.eos(e), star.Mrot, star.Req, star.tidal])


            #scattered.append(tmp)
            
            rhoc = numpy.random.rand() *(numpy.log10(EOS.max_edsc) - 14.6) + 14.6
            star = Star(10**rhoc)
            star.solve_structure(EOS.energydensities, EOS.pressures)
            MR_prpr_pp[i] = star.Mgrav, star.Req

            radii[:,i] = MR(masses)

            bar()


    numpy.save(root_name  + 'radii',radii)
    numpy.save(root_name + 'pressures', pressures)
    numpy.savetxt(root_name + 'MR_prpr.txt', MR_prpr_pp)

    minpres, maxpres = calc_bands(energydensities, pressures)
    #minpres_rho, maxpres_rho = calc_bands(energydensities, pressures_rho)
    minradii, maxradii = calc_bands(masses, radii)



    # save everything
    numpy.save(root_name + 'minradii', minradii)
    numpy.save(root_name + 'maxradii', maxradii)
    numpy.save(root_name + 'minpres', minpres)
    numpy.save(root_name + 'maxpres', maxpres)



root_name = run_name
root_name_B = run_name

compute_auxiliary_data(root_name, EOS, variable_params, static_params, chirp_mass)
