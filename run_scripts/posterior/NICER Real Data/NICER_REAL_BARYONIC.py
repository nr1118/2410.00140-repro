#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import neost
from neost.eos import polytropes
from neost.Prior import Prior
from neost.Star import Star
from neost.Likelihood import Likelihood
from neost import PosteriorAnalysis
from scipy.stats import gaussian_kde
import numpy as np
from pymultinest.solve import solve
import time
import seaborn as sns
import corner as corner



# In[ ]:


import neost.global_imports as global_imports

c = global_imports._c
G = global_imports._G
Msun = global_imports._M_s
pi = global_imports._pi
rho_ns = global_imports._rhons


eos_name = 'polytropes'


EOS = polytropes.PolytropicEoS(crust = 'ceft-Hebeler',rho_t = 2e14, adm_type = 'None')

# Here we implement old NICER data on J0740 and J0030 from Riley et al.


# PSR J0740+6620 (arXiv: 2105.06980)
muM = 2.07  
muR = 12.39
sigM = 0.07  #published uncertainty
sigR_plus = 1.30
sigR_minus = 0.98
J0740_mr_samples = np.load('Riley00740_mr_posterior_samples.npy').T
J0740_kde = gaussian_kde(J0740_mr_samples)

# PSR J0030+0451(arXiv: 1912.05702)
muM2 = 1.34 
muR2 = 12.71
sigM2 = 0.15 #published uncertainty
sigR2_plus = 1.13
sigR2_minus = 1.18
J0030_mr_samples = np.load('Riley0030_mr_posterior_samples.npy').T
J0030_kde = gaussian_kde(J0030_mr_samples)


likelihood_functions = [J0740_kde.pdf,J0030_kde.pdf]
likelihood_params = [['Mass', 'Radius'],['Mass','Radius']]

# This is not a GW event so we set chirp mass to None
chirp_mass = [None,None]
number_stars = len(chirp_mass)

run_name = "NICER_REAL_BARYONIC_"


variable_params = {'gamma1':[1., 4.5], 'gamma2':[0., 8.], 'gamma3':[0.5, 8.], 'rho_t1':[1.5, 8.3], 'rho_t2':[1.5, 8.3],'ceft':[EOS.min_norm, EOS.max_norm]}

for i in range(number_stars):
	variable_params.update({'rhoc_' + str(i+1):[14.6, 16]})




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
cube = np.random.rand(50, len(variable_params))
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


PosteriorAnalysis.compute_minimal_auxiliary_data_Baryonic(run_name, EOS,
                                         variable_params, static_params, chirp_mass)