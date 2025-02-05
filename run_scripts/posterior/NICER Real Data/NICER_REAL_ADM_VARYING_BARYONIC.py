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

run_name = "NICER_REAL_ADM_VARYING_BARYONIC_"

#posterior on mchi is influences such that we don't same mchi < 100 (10^2) see prior corner plots
variable_params = {'gamma1':[1., 4.5], 'gamma2':[0., 8.], 'gamma3':[0.5, 8.], 'rho_t1':[1.5, 8.3], 'rho_t2':[1.5, 8.3],
                  'mchi':[2, 9],'gchi_over_mphi': [-5,3],'adm_fraction':[0., 1.7],'ceft':[EOS.min_norm, EOS.max_norm]}
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
result = solve(LogLikelihood=likelihood.call, Prior=prior.inverse_sample, n_live_points=3000, evidence_tolerance=0.1,
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

    epscent_dm = EOS.find_epsdm_cent(par['adm_fraction'],rhocpar)
    star = Star(rhocpar,epscent_dm)
    # star.solve_structure(tabulated_example.energydensities, tabulated_example.pressures)
    star.solve_structure(EOS.energydensities, EOS.pressures,
                             EOS.energydensities_dm, EOS.pressures_dm) 
                             #ADM_fraction=par['adm_fraction'], mchi=par['mchi'], gchi_over_mphi=par['gchi_over_mphi'])
        
    scattered.append([rhocpar, EOS.eos(rhocpar), star.Mgrav, star.Req,star.tidal,star.Rdm_halo]) ## if Rdm_halo == 999 then we have a halo if Rdm_halo ==0 then we have a core



       
 

scattered = numpy.array(scattered)
    


#print(scattered)
numpy.savetxt(run_name + 'scattered.txt', scattered)
print("Testing done")

print('Generating the posterior corner plot')
C = numpy.loadtxt(run_name + 'post_equal_weights.dat')
mchi = C[:,5]
gchi_over_mphi = C[:,6]
Fchi = C[:,7]


Matrix = numpy.zeros((len(mchi),3))
###### Everything below this needs to be fixed!####
for i in range(len(mchi)):
    Matrix[i] = numpy.log10(mchi[i]),numpy.log10(gchi_over_mphi[i]),Fchi[i]


scattered_prior = numpy.loadtxt('PRIOR_FERMIONIC_REAL_07nsaturation_scattered.txt')

mchi_prior = numpy.log10(scattered_prior[:,0])
gchi_over_mphi_prior= numpy.log10(scattered_prior[:,1])
Fchi_prior= scattered_prior[:,2]

Matrix_prior = numpy.zeros((len(scattered_prior),3)) 

for i in range(len(scattered_prior)):
    if scattered_prior[i][-1] == 0 and scattered_prior[i][5]>=1.:
        Matrix_prior[i] = mchi_prior[i], gchi_over_mphi_prior[i], Fchi_prior[i]       
    
ell = corner.corner(Matrix_prior,smooth = 1.0,color = '#377eb8',group = 'prior',range = [(2,9),(-5,3),(0,1.7)],
                   plot_datapoints = False,plot_density = False,plot_contours = False,divergences = False,
                    hist_kwargs = {'linestyle': '--'})

figure = corner.corner(Matrix,smooth = 1.0,fig = ell,labels = [r"$\mathrm{log_{10}}(m_\chi/MeV)$",r"$\mathrm{log_{10}}(\frac{g_\chi}{m_\phi/MeV})$",r"$F_\chi$"]
                       ,color = 'black',range = [(2,9),(-5,3),(0,1.7)],
                       show_titles = True,label_kwargs = {"fontsize":20},title_kwargs = {"fontsize":15},
                       hist_kwargs = {'linewidth':1.0})

#quantiles =(0.001,0.999)
figure.subplots_adjust(right=1.15,top=1.15)
for ax in figure.get_axes():
    ax.tick_params(axis='both', labelsize=15) 
    
figure.legend(handles =[matplotlib.lines.Line2D([],[],color = 'black',label = 'Posterior'),
                        matplotlib.lines.Line2D([],[],color = '#377eb8',label = 'Prior',linestyle = '--')],
          fontsize = 20,frameon = False,loc = "upper right")
    
figure.savefig(run_name + 'Corner.png',dpi = 300,bbox_inches='tight')
figure.savefig(run_name + 'Corner.pdf',dpi = 300,bbox_inches='tight')
figure.savefig(run_name + 'Corner.jpeg',bbox_inches='tight')


