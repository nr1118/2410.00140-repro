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
import corner as corner
from alive_progress import alive_bar


# In[ ]:


import neost.global_imports as global_imports

c = global_imports._c
G = global_imports._G
Msun = global_imports._M_s
pi = global_imports._pi
rho_ns = global_imports._rhons


eos_name = 'polytropes'


EOS = polytropes.PolytropicEoS(crust = 'ceft-Hebeler', rho_t = 2e14, adm_type = 'Fermionic')


# EOS.plot()
# EOS.plot_massradius()
# Here we implement old NICER data on J0740 and J0030 from Riley et al.


# PSR J0740+6620 (arXiv: 2105.06980)
muM = 2.07  
muR = 12.39
sigM = 0.07  #published uncertainty
sigR_plus = 1.30
sigR_minus = 0.98


# PSR J0030+0451(arXiv: 1912.05702)
muM2 = 1.34 
muR2 = 12.71
sigM2 = 0.15 #published uncertainty
sigR2_plus = 1.13
sigR2_minus = 1.18



likelihood_functions = []
likelihood_params = [['Mass', 'Radius']]

# This is not a GW event so we set chirp mass to None
chirp_mass = [None]
number_stars = len(chirp_mass)

run_name = "FERMIONIC_REAL_DATA_PRIOR_"


variable_params = {'gamma1':[1., 4.5], 'gamma2':[0., 8.], 'gamma3':[0.5, 8.], 'rho_t1':[1.5, 8.3], 'rho_t2':[1.5, 8.3],
                  'mchi':[0, 9],'gchi_over_mphi': [-5,3],'adm_fraction':[0., 1.7],'ceft':[EOS.min_norm, EOS.max_norm]}
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


# Then we start the sampling, note the greatly increased number of livepoints, this is required because each livepoint terminates after 1 iteration
start = time.time()
result = solve(LogLikelihood=likelihood.loglike_prior, Prior=prior.inverse_sample, n_live_points=15000, evidence_tolerance=0.1,
               n_dims=len(variable_params), sampling_efficiency=0.8, outputfiles_basename=run_name, verbose=True)
end = time.time()
print(end - start)

# In[ ]:

print('Solving done, moving to Prior Analysis')
ewposterior = numpy.loadtxt(run_name +'post_equal_weights.dat')

energydensities = numpy.logspace(14.2, 16, 50)
pressures = numpy.zeros((len(energydensities), len(ewposterior)))
pressures_rho = numpy.zeros((len(energydensities), len(ewposterior)))
scattered = numpy.zeros((len(ewposterior),8))

with alive_bar(len(ewposterior)) as bar:
    for i in range(0, len(ewposterior), 1):
        pr = ewposterior[i][0:len(variable_params)]
        par = {e:pr[j] for j, e in enumerate(list(variable_params.keys()))}
        par.update(static_params)
    
        EOS.update(par, max_edsc=True)
        rhopres = UnivariateSpline(EOS.massdensities, EOS.pressures, k=1, s=0)
        edsrho = UnivariateSpline(EOS.energydensities,EOS.massdensities, k=1, s=0)
        max_rhoc = edsrho(EOS.max_edsc)
        pressures_rho[:,i][energydensities<max_rhoc] = rhopres(energydensities[energydensities<max_rhoc])
        pressures[:,i][energydensities<EOS.max_edsc] = EOS.eos(energydensities[energydensities<EOS.max_edsc])  
     
        rhocpar = 10**par['rhoc_1']

        epscent_dm = EOS.find_epsdm_cent(par['adm_fraction'],rhocpar)
        star = Star(rhocpar,epscent_dm)
    # star.solve_structure(tabulated_example.energydensities, tabulated_example.pressures)
        star.solve_structure(EOS.energydensities, EOS.pressures,
                                EOS.energydensities_dm, EOS.pressures_dm) 
                             #ADM_fraction=par['adm_fraction'], mchi=par['mchi'], gchi_over_mphi=par['gchi_over_mphi'])
        
        scattered[i] = par['mchi'], par['gchi_over_mphi'], par['adm_fraction'],rhocpar, EOS.eos(rhocpar), star.Mgrav, star.Req,star.Rdm_halo






numpy.save(run_name + 'Pressure_array.npy',pressures)
numpy.save(run_name + 'Pressures_Rho_array.npy',pressures_rho)

#print(scattered)
numpy.savetxt(run_name + 'scattered.txt', scattered)

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


minpres, maxpres = calc_bands(energydensities, pressures)
minpres_rho, maxpres_rho = calc_bands(energydensities, pressures_rho)
numpy.save(run_name + 'minpres_rho', minpres_rho)
numpy.save(run_name + 'maxpres_rho', maxpres_rho)
numpy.save(run_name + 'minpres', minpres)
numpy.save(run_name + 'maxpres', maxpres)