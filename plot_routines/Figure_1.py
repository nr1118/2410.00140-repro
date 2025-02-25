#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import corner as corner
import seaborn as sns
import matplotlib as plt
from matplotlib import pyplot
import matplotlib.patches as mpatches

import os
import pathlib
import argparse
import sys


# In[2]:


plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.family'] ='serif'


# In[3]:

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--repro', action='store_true')
args = parser.parse_args()


data_directory = f'../results/prior/' if not args.repro else f'../repro/prior/'

plots_directory = '../plots/' 

run_name = 'FERMIONIC_REAL_DATA_PRIOR_'


# In[4]:


tmp = np.loadtxt(data_directory + 'FERMIONIC_REAL_DATA_PRIOR_post_equal_weights.dat')
print('Generating the prior corner plot')


Matrix_prior = np.zeros((len(tmp),3))

#m_chi = tmp[:,5]
#g_chi = tmp[:,6]
#F_chi = tmp[:,7]

for i in range(len(tmp)):
    Matrix_prior[i] =np.log10(tmp[:,5][i]),np.log10(tmp[:,6][i]),tmp[:,7][i]

    


figure = corner.corner(Matrix_prior,smooth = 1.0,labels = [r"log$_{10}$(m$_\chi$/MeV)",r"log$_{10}$($\frac{\mathdefault{g}_\chi}{\mathdefault{m}_\phi/\mathdefault{MeV}})$",r"F$_\chi$ [$\%$]"],
                      range = [(1,9),(-5,3),(0,1.7)], show_titles = True,label_kwargs = {"fontsize":18,"font":'serif'},title_kwargs = {"fontsize":15},
                      color = '#377eb8',hist_kwargs = {'linestyle': '--','linewidth': 2.0}, contour_kwargs = {'linestyles':'dashed','linewidths': 2.0} )



figure.subplots_adjust(right=1.15,top=1.15)

for ax in figure.get_axes():
    ax.tick_params(axis='both', labelsize=15)
    
figure.savefig(plots_directory + run_name + 'Corner.png',bbox_inches='tight')


# In[ ]:




