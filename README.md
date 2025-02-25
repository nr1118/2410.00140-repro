# 2410.00140-repro
Reproduction package for arXiv:2410.00140


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


TODO: Write the README file


This repository contains all scripts necessary to reproduce all results from the paper "Probing fermionic asymmetric dark matter cores using global neutron star properties" by Rutherford et al., PRD. (2025) [preprint:  https://arxiv.org/abs/2410.00140 | doi: insert doi from PRD when the time comes].

REQUIREMENTS
============
For reproducing plots using the supplied data in the results/ directory, the following 3rd-party python libraries are required:

  - matplotlib
  - numpy
  - seaborn
  - scipy
  - corner

To reproduce our data (see Sampling below), the following additional requirements apply (see also NEoST installation instructions):

  - neost (https://github.com/xpsi-group/neost, version 2.1.0)
  - A modified PosteriorAnalysis.py script (see Modified_neost_file/ directory)
  - cython
  - pymultinest
  - gsl
  - kalepy
  - numba (optional)

numba is only needed to use the Python TOV solvers; normally the much faster Cython TOV solvers should be used. While the modified PosteriorAnalysis.py file is nearly identical to the one used in NEoST v2.1.0, but containts two additional functions, namely compute_minimal_auxliarly_data_ADM() and compute_minimal_auxliarly_data_Baryonic(). These functions are simply minimal versions of the compute_auxliarly_data() function for ADM and Baryonic matter, which only compute the necessary data used in the manuscript. As opposed to additoinal data that would be computed in the original compute_auxliarly_data() function in NEoST v2.1.0. Please see the README within Modified_neost_file/ directory for installation instructions on the modified PosteriorAnalysis.py script. 

REPRODUCING PLOTS
=================
All figures in the paper can be reproduced by simply executing the generate_figs.sh script or by going to the plot_routines folder and running each python script within the directory. Refer to this script and the called plot scripts in plot_routines/ for a complete account of all options available.

The most important option is the -r (--repro) flag, which almost all scripts recognize. By default, generate_figs.sh uses data supplied in the results/ directory---which contains the results published in the paper---to produce figures. The -r flag tells the plot scripts to instead use user-generated data in the repro/ directory. If you wish to use your own generated data, please ensure it is in the repro/ directory, then the generate_figs.sh can recognize the following options, which are recognized by all plotting scripts in plot_routines/ directory:

  - name_prior: run name for the reproduced prior results instead of published
  - name_posterior_incladm_real: run name for the reproduced real data including ADM results instead of published
  - name_posterior_negladm_real: run name for the reproduced real data neglecting ADM results instead of published
  - name_posterior_incladm_adm: run name for the reproduced Future-X ADM core model including ADM results instead of published
  - name_posterior_negladm_adm: run name for the reproduced Future-X ADM core model neglecting ADM results instead of published
  - name_posterior_incladm_noadm: run name for the reproduced Future-X No ADM model including ADM results instead of published
  - name_posterior_negladm_noadm: run name for the reproduced Future-X No ADM model neglecting ADM results instead of published
  - name_appendix: run name for the reproduced Appendix plots of approximating zero self-repulsion

To generate figures using the results from the paper:

  ./generate_figs.sh

To instead use your own data:

  ./generate_figs.sh -r

Note that you may be required to set the executable permission for generate_figs.sh with

  chmod +x generate_figs.sh


RUN SCRIPTS
========
The run scripts for all prior, posterior, and Appendix B calculations can be found in the run_scripts/ directory and can be used to reproduce all of the results in the paper. You need to adapt these scripts if you want to use them with custom run names and/or output directories. Furthermore, if you wish to use the plotting scripts with these runs, they must be in the repro/ directory. The overall structure of the run_scripts folder is as follows:

- run_scripts/posterior: All posterior script files
      - /NICER_Real_Data: The real data posterior scripts which use the Riley et al. 2019 & 2021 MR inferences. Here, the NICER_REAL_ADM_VARYING.py script is the script which varies both the Baryonic and ADM equation of state parameters, whereas the NICER_REAL_BARYONIC.py script is the one which neglects the ADM EoS and only considers the Baryonic matter EoS model.
      -  /Future-X: The Future-X data posterior scripts which consider the synthetically generated MR values corresponding to the ADM core model and No ADM models. This directory is further split into the ADM Core Model/ and No ADM Model/ directories. Within the ADM Core Model/ the FUTUREX_ADM_VARYING_BARYONIC.py is the script which samples both the ADM and baryonic matter EoS parameters using the ADM core model sources, while the FUTUREX_ADM_BARYONIC_ONLY.py script is the oone which neglected the ADM EoS parameters and only samples the baryonic matter EoS parameters using the ADM core model sources. For the No ADM Model/ directory it is the same as the ADM Core Model/ directory, except now for the No ADM Model sources. Thus the scripts for the including ADM and neglecting ADM scenarios are FUTUREX_NO_ADM_VARYING_BARYONIC.py and FUTUREX_NO_ADM_BARYONIC_ONLY.py, respectively

-  run_scripts/prior: All prior script files
      -  Containts the prior script which varies the baryonic matter and ADM EoS parameters named FERMIONIC_PRIOR.py.
      -  /Appendix B: Contains the run scripts to generate the data used in Appendix B, which is the approximation of zero self-repulsion.
