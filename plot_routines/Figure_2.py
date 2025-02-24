#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib as plt
from matplotlib import pyplot


# In[2]:


from matplotlib.patches import Ellipse 
import pathlib


# In[3]:


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


# In[4]:


data_directory = '../data/'

plots_directory = '../plots/' 


# In[5]:


noadm = np.load(data_directory + 'NO_ADM_ground_truth.npy')
M_B = noadm[:,0]
R_B = noadm[:,1]

adm = np.load(data_directory + 'ADM_Core_ground_truth.npy')
M_ADM = adm[:,0]
R_ADM = adm[:,1]


# In[ ]:





# In[6]:


# central energy density: 15.25140196
muM_strobex = 2.232  
muR_strobex = 10.999
sigM_strobex = muM_strobex * 0.02   #2% uncertainty
sigR_strobex = muR_strobex * 0.02 #2% uncertainty in radius


# central energy density: 15.16132928
muM2_strobex = 2.001
muR2_strobex = 11.418
sigM2_strobex = muM2_strobex * 0.02    #2% uncertainty in mass
sigR2_strobex = muR2_strobex * 0.02   #2% uncertainty in radius


# central energy density: 15.10621903
muM3_strobex = 1.886   
muR3_strobex = 11.474
sigM3_strobex = muM3_strobex * 0.02   #2% uncertainty in mass
sigR3_strobex = muR3_strobex * 0.02  # 2% uncertainty in radius 


# central energy density: 15.07244177
muM4_strobex = 1.636   
muR4_strobex = 11.524
sigM4_strobex = muM4_strobex * 0.02 # 2 % uncerainty in mass 
sigR4_strobex = muR4_strobex * 0.02  # 2 % uncertainty in radius 



# central energy density: 15.02562768
muM5_strobex = 1.457   
muR5_strobex = 11.522
sigM5_strobex = muM5_strobex * 0.02   # 2% uncertainty in mass
sigR5_strobex = muR5_strobex * 0.02  # 2 % uncertainty in radius 


# central energy density: 14.94562892
muM6_strobex = 1.263    
muR6_strobex = 11.499
sigM6_strobex = muM6_strobex * 0.02   # 2% uncertainty in mass
sigR6_strobex = muR6_strobex * 0.02  # 2 % uncertainty in radius 


# In[7]:


# central energy density: 15.12103361
muM_noADMstrobex = 2.232  
muR_noADMstrobex = 11.285
sigM_noADMstrobex = muM_noADMstrobex * 0.02   #2% uncertainty
sigR_noADMstrobex = muR_noADMstrobex * 0.02   #2% uncertainty in radius


# central energy density: 15.08251569
muM2_noADMstrobex = 2.001
muR2_noADMstrobex = 11.557
sigM2_noADMstrobex =muM2_noADMstrobex * 0.02      #2% uncertainty
sigR2_noADMstrobex = muR2_noADMstrobex * 0.02     #2% uncertainty in radius


# central energy density: 15.05940494
muM3_noADMstrobex = 1.886   
muR3_noADMstrobex = 11.606
sigM3_noADMstrobex = muM3_noADMstrobex * 0.02   #2% uncertainty
sigR3_noADMstrobex = muR3_noADMstrobex * 0.02  #2% uncertainty in radius 


# central energy density: 15.03747935
muM4_noADMstrobex = 1.636   
muR4_noADMstrobex = 11.647
sigM4_noADMstrobex = muM4_noADMstrobex * 0.02# 2 % uncerainty in mass 
sigR4_noADMstrobex = muR4_noADMstrobex * 0.02  #2% uncertainty in radius 



# central energy density: 14.99896143
muM5_noADMstrobex = 1.457    
muR5_noADMstrobex = 11.641
sigM5_noADMstrobex = muM5_noADMstrobex * 0.02       #2% uncertainty in mass
sigR5_noADMstrobex = muR5_noADMstrobex * 0.02      #2% uncertainty in radius 



# central energy density: 14.928444
muM6_noADMstrobex = 1.263    
muR6_noADMstrobex = 11.616
sigM6_noADMstrobex = muM6_noADMstrobex * 0.02        #2% uncertainty
sigR6_noADMstrobex = muR6_noADMstrobex * 0.02        #2% uncertainty in radius 


# In[ ]:





# In[ ]:





# In[ ]:





# In[8]:



ells_strobex = [Ellipse((muR_strobex,muM_strobex),width = 2*sigR_strobex,height = 2*sigM_strobex
                        ,angle = 0,linestyle = '--',fill = False,color = 'C1',lw = 2.0),
        Ellipse((muR_noADMstrobex,muM_noADMstrobex),width = 2*sigR_noADMstrobex,height = 2*sigM_noADMstrobex
                ,angle = 0,linestyle = '-',fill = False,color = '#377eb8',lw = 2.0)]
fig, ax = pyplot.subplots(1,1, figsize=(12, 10))
#STROBE-X stars
stars_strobex = [Ellipse((muR2_strobex,muM2_strobex),width = 2*sigR2_strobex,height = 2*sigM2_strobex
                         ,angle = 0,linestyle = '--',fill = False,color = 'C1',lw = 2.0),
                 Ellipse((muR3_strobex,muM3_strobex),width = 2*sigR3_strobex,height = 2*sigM3_strobex
                         ,angle = 0,linestyle = '--',fill = False,color = 'C1',lw = 2.0),
        Ellipse((muR4_strobex,muM4_strobex),width = 2*sigR4_strobex,height = 2*sigM4_strobex
                ,angle = 0,linestyle = '--',fill = False,color = 'C1',lw = 2.0),
                 Ellipse((muR5_strobex,muM5_strobex),width = 2*sigR5_strobex,height = 2*sigM5_strobex
                         ,angle = 0,linestyle = '--',fill = False,color = 'C1',lw = 2.0)
                 ,Ellipse((muR6_strobex,muM6_strobex),width = 2*sigR6_strobex,height = 2*sigM6_strobex,angle = 0
                          ,linestyle = '--',fill = False,color = 'C1',lw = 2.0) ]

stars_noADMstrobex = [Ellipse((muR2_noADMstrobex,muM2_noADMstrobex),width = 2*sigR2_noADMstrobex,
                              height = 2*sigM2_noADMstrobex,angle = 0,linestyle = '-',fill = False
                      ,color = '#377eb8',lw = 2.0),Ellipse((muR3_noADMstrobex,muM3_noADMstrobex),width = 2*sigR3_noADMstrobex,
                                                  height = 2*sigM3_noADMstrobex,angle = 0
                                                  ,linestyle = '-',fill = False,color = '#377eb8',lw = 2.0),
        Ellipse((muR4_noADMstrobex,muM4_noADMstrobex),width = 2*sigR4_noADMstrobex,height = 2*sigM4_noADMstrobex
                ,angle = 0,linestyle = '-',fill = False,color = '#377eb8',lw = 2.0),Ellipse((muR5_noADMstrobex,muM5_noADMstrobex),
                width = 2*sigR5_noADMstrobex,height = 2*sigM5_noADMstrobex,angle = 0,linestyle = '-',fill = False
                      ,color = '#377eb8',lw = 2.0),Ellipse((muR6_noADMstrobex,muM6_noADMstrobex),width = 2*sigR6_noADMstrobex
                        ,height = 2*sigM6_noADMstrobex,angle = 0,linestyle = '-',fill = False,color = '#377eb8',lw = 2.0)]



for i in range(len(stars_strobex)):
    ax.add_patch(stars_strobex[i])
    ax.add_patch(stars_noADMstrobex[i])
for i, e in enumerate(ells_strobex):
    ax.add_patch(e)
    
    
ax.legend(ells_strobex, ['ADM Core Model','No ADM Model'],fontsize = 26)
ax.set_title(r'FUTURE-X',fontsize = 28)
ax.tick_params(top = 1, right = 1, axis='both', which='both',direction = 'in', labelsize=22)
ax.plot(R_ADM,M_ADM, linestyle = '--',color = 'C1',lw = 2.0)
ax.plot(R_B,M_B, linestyle = '-',color = '#377eb8',lw = 2.0)
xlim = (10, 14)
ylim = (1.0, 2.4)
ax.set_xlim(xlim)
ax.set_ylim(ylim)

# Ticks
xticks = [ 11,12,13, 14,15, 16]
yticks = [1.0, 1.4, 1.8, 2.2]
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.minorticks_on()
ax.set_ylabel(r'$\mathrm{Mass \, [M_\odot]}$',fontsize = 28)
ax.set_xlabel(r'$\mathrm{Radius \, [km]}$',fontsize = 28)


fig.savefig(plots_directory + 'Ellipses_FUTUREX.pdf',bbox_inches='tight')
#plt.show()
    


# In[9]:


J0030_R19_STPST_MR_95 = np.loadtxt(data_directory + f'J0030_Riley19_NICERonly_STPST_95.txt')
J0030_R19_STPST_M_95 = J0030_R19_STPST_MR_95[:,1]
J0030_R19_STPST_R_95 = J0030_R19_STPST_MR_95[:,0]

J0030_R19_STPST_MR_68 = np.loadtxt(data_directory + f'J0030_Riley19_NICERonly_STPST_68.txt')
J0030_R19_STPST_M_68 = J0030_R19_STPST_MR_68[:,1]
J0030_R19_STPST_R_68 = J0030_R19_STPST_MR_68[:,0]


# In[10]:


J0740_R21_STU_MR_68 = np.loadtxt(data_directory + f'J0740_Riley21_68.txt')
J0740_R21_STU_MR_95 = np.loadtxt(data_directory + f'J0740_Riley21_95.txt')

J0740_R21_STU_M_68 = J0740_R21_STU_MR_68[:,1]
J0740_R21_STU_R_68 = J0740_R21_STU_MR_68[:,0]
J0740_R21_STU_M_95 = J0740_R21_STU_MR_95[:,1]
J0740_R21_STU_R_95 = J0740_R21_STU_MR_95[:,0]


# In[12]:


fig, axes = pyplot.subplots(1, 1, figsize=(12,10))

ax_baseline = axes
ax_baseline.annotate('J0030', xy=(11.0, 1.4), color='xkcd:aqua',fontsize = 24)
ax_baseline.annotate('J0740', xy=(11.0, 2.25), color='xkcd:light violet',fontsize = 24)
ax_baseline.tick_params(top = 1, right = 1, which="both", axis="both", direction="in",labelsize=22 )

xlim = (10, 16)
ylim = (1.0, 2.4)
ax_baseline.set_xlim(xlim)
ax_baseline.set_ylim(ylim)

# Ticks
xticks = [11, 12,13, 14,15, 16]
yticks = [1.0, 1.4, 1.8, 2.2]
ax_baseline.set_xticks(xticks)
ax_baseline.set_yticks(yticks)
ax_baseline.minorticks_on()


axes.set_title(r'REAL DATA',fontsize = 28)
axes.tick_params(axis='both', which='major', labelsize=22)
axes.set_ylabel(r'$\mathrm{Mass \, [M_\odot]}$',fontsize = 28)
axes.set_xlabel(r'$\mathrm{Radius \, [km]}$',fontsize = 28)

ax_baseline.fill(J0030_R19_STPST_R_95, J0030_R19_STPST_M_95, linewidth=2.0, color='xkcd:aqua', alpha=0.3)
ax_baseline.fill(J0030_R19_STPST_R_68, J0030_R19_STPST_M_68, linewidth=2.0, color='xkcd:aqua', alpha=0.6)


ax_baseline.fill(J0740_R21_STU_R_95, J0740_R21_STU_M_95, linewidth=2.0, color='xkcd:light violet', alpha=0.4)
ax_baseline.fill(J0740_R21_STU_R_68, J0740_R21_STU_M_68, linewidth=2.0, color='xkcd:light violet', alpha=0.8)
fig.savefig(plots_directory + 'MRdata.png',bbox_inches='tight')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




