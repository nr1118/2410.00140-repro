# 2410.00140-repro
Reproduction package for arXiv:2410.00140


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


This repository contains all scripts necessary to reproduce all results from the paper "Probing fermionic asymmetric dark matter cores using global neutron star properties" by Rutherford et al., PRD. (2025) [preprint:  https://arxiv.org/abs/2410.00140 | doi: insert doi from PRD when the time comes].

NOTICE
======================
Run times range from about day for runs with posteriors which only varying the baryonic matter EoS to about 10 days for the priors/posteriors which additionally include the ADM EoS on a supercomputing cluster. Therefore, it is NOT recommended that the prior and posterior scripts be run on a personal computer or laptop due to the computational demand of the runs. However, small test runs with a about 100 live points can run on a personal computer or laptop to ensure the script is running as desired. 

REQUIREMENTS
============
To reproduce our data and plots, the following additional requirements apply (see also NEoST installation instructions):

  - neost (https://github.com/xpsi-group/neost, version 2.1.0)
  - Modified PosteriorAnalysis.py with minimum auxiliary data function script (see below and Modified_neost_file/ directory for further detail)
  - cython
  - pymultinest
  - gsl
  - kalepy
  - scipy
  - matplotlib
  - seaborn
  - numpy
  - corner
  - numba (optional)

The numba is only needed to use the Python TOV solvers. However, normally the much faster Cython TOV solvers should be used (see https://xpsi-group.github.io/neost/index.html for details). While the modified PosteriorAnalysis.py file is nearly identical to the one used in NEoST v2.1.0, but containts two additional functions, namely compute_minimal_auxliarly_data_ADM() and compute_minimal_auxliarly_data_Baryonic(). These functions are simply minimal versions of the compute_auxliarly_data() function for ADM and Baryonic matter, which only compute the necessary data used in the manuscript. As opposed to additoinal data that would be computed in the original compute_auxliarly_data() function in NEoST v2.1.0. Please see the README within Modified_neost_file/ directory for installation instructions on the modified PosteriorAnalysis.py script. 

REPRODUCING PLOTS
=================
All figures in the paper can be reproduced by simply executing the generate_figs.sh script or by going to the plot_routines/ directory and running each python script within the directory. Note, if reproducing the plots via the plot_routines/ directory, the neost conda enviroment must be activated. Refer to this script and the called plot scripts in plot_routines/ for a complete account of all options available. The output of the plot_routines/ directory goes directory to the plots/ directory. Moreover, some of the scripts in the plots_routines/ directory pull from the data/ directory, which containts the mass-radius contours for J0030 and J0740 as well as their posterior samples, which are used in the run_scripts. The data/ directory also containts the mass radius curves for the ADM core and No ADM models as well, which are clearly labeled.

The most important option is the -r (--repro) flag, which all scripts recognize. By default, generate_figs.sh uses data supplied in the results/ directory---which contains the results published in the paper---to produce figures. The -r flag tells the plot scripts to instead use user-generated data in the repro/ directory. If you wish to use your own generated data, please ensure it is in the repro/ directory, then the generate_figs.sh can recognize the following options, which are recognized by all plotting scripts in plot_routines/ directory:

  - name_prior: run name for the reproduced prior results instead of published
  - name_posterior_incladm_real: run name for the reproduced real data including ADM results instead of published
  - name_posterior_negladm_real: run name for the reproduced real data neglecting ADM results instead of published
  - name_posterior_incladm_adm: run name for the reproduced Future-X ADM core model including ADM results instead of published
  - name_posterior_negladm_adm: run name for the reproduced Future-X ADM core model neglecting ADM results instead of published
  - name_posterior_incladm_noadm: run name for the reproduced Future-X No ADM model including ADM results instead of published
  - name_posterior_negladm_noadm: run name for the reproduced Future-X No ADM model neglecting ADM results instead of published
  - name_appendix: run name for the reproduced Appendix plots of approximating zero self-repulsion

To generate figures using the results from the paper using generate_figs.sh:

  ./generate_figs.sh

To instead use your own data:

  ./generate_figs.sh -r -name_prior -name_posterior_incladm_real ... -name_posterior_negladm_noadm

  Note, the other name tags to repro data are only there if you wish to run your own sampling for all cases used in the manuscript. Less name tags can be used depending on what scenarios you wish to do your own inference on.
  
 You may be required to set the executable permission for generate_figs.sh with

  chmod +x generate_figs.sh


RUN SCRIPTS
========
The run scripts for all prior, posterior, and Appendix B calculations can be found in the run_scripts/ directory and can be used to reproduce all of the results in the paper. You need to adapt these scripts if you want to use them with custom run names and/or output directories. Furthermore, if you wish to use the plotting scripts with these runs, they must be in the repro/ directory. The overall structure of the run_scripts folder is as follows:

- run_scripts/posterior: All posterior script files
      - /NICER_Real_Data: The real data posterior scripts which use the Riley et al. 2019 & 2021 MR inferences. Here, the NICER_REAL_ADM_VARYING.py script is the script which varies both the Baryonic and ADM equation of state parameters, whereas the NICER_REAL_BARYONIC.py script is the one which neglects the ADM EoS and only considers the Baryonic matter EoS model.
      -  /Future-X: The Future-X data posterior scripts which consider the synthetically generated MR values corresponding to the ADM core model and No ADM models. This directory is further split into the ADM Core Model/ and No ADM Model/ directories. Within the ADM Core Model/ the FUTUREX_ADM_VARYING_BARYONIC.py is the script which samples both the ADM and baryonic matter EoS parameters using the ADM core model sources, while the FUTUREX_ADM_BARYONIC_ONLY.py script is the oone which neglected the ADM EoS parameters and only samples the baryonic matter EoS parameters using the ADM core model sources. For the No ADM Model/ directory it is the same as the ADM Core Model/ directory, except now for the No ADM Model sources. Thus the scripts for the including ADM and neglecting ADM scenarios are FUTUREX_NO_ADM_VARYING_BARYONIC.py and FUTUREX_NO_ADM_BARYONIC_ONLY.py, respectively

-  run_scripts/prior: All prior script files
      -  Containts the prior scripts which varies the baryonic matter and ADM EoS parameters named FERMIONIC_PRIOR.py.
      -  /Appendix B: Contains the run scripts to generate the data used in Appendix B, which is the approximation of zero self-repulsion using 10^{-5} MeV^{-1}.
 

OUTPUT FILE STRUCTURES
======================
- results/: The first level is either "prior" or "posterior". For the prior, this is followed by the output files of the runs which vary both the baryonic matter and ADM EoS, namely the output files ending with post_equal_weights.dat, MR_prpr.txt, and pressures.npy.
      - post_equal_weights.dat: Standard Multinest output file, which containts the admixed (baryonic + ADM) EoS model parameters of the PP model and fermionic ADM EoS, the sampled central density, and log-likelihood evaluation of each sample.
      - MR_prpr.txt: Standard NEoST output file, which containts the corresponding mass and radius samples of each sampled admixed EoS and central density.
      - pressures.npy: Stnadard NEoST output file, which containts the pressures of the baryonic EoS in which ADM was considered during prior sampling.
  Moreover, the prior directory is also followed by two other directories:
      - Appendix_B/: Relative radial percent differences comparing zero ADM self-repulsion to 10^{-5} MeV^{-1}. The tail end of the file indicates which choice of baryonic EoS is used and if a different step of ADM particle mass or mass-fraction were used. For example, the file Relcent_diff_intermediate_stiff_baryonic_newfchistep.npy are the results which used the intermediately stiff baryonic EoS in the manuscript with a smaller step in fchi compared to the Relcent_diff_intermediate_stiff_baryonic.npy file.
       - Fchi_prior_calculation/: Directory containing the scripts that print the estimated ADM mass-fraction and distances to the Galactic center for J0437 and J0030. Note, these are run scripts because the results are simply printed results and not stored in a file structure.

  For the posterior, the directory is followed by two more directories, Future-X/ and NICER_Real_Data/, which are the posteriors corresponding to the Future-X and Real data inferences in the manuscript, respectfully.
    - Future-X/: Followed by two more directories for the ADM_core_model/ and No_ADM_Model/, which of course correspond to the ADM core model and No ADM model, respectfully. Both of these directories then have the output directories for the posteriors in which ADM is included (FUTUREX_ADM_VARYING_BARYONIC for the ADM core model and FUTUREX_NO_ADM_VARYING_BARYONIC for the No ADM model) and neglected from sampling (FUTUREX_ADM_BARYONIC_ONLY for the ADM core model and FUTUREX_NO_ADM_BARYONIC_ONLY for the No ADM model). These output direcories have the post_equal_weights.dat, MR_prpr.txt, and pressures.npy files. Furthermore, they also have the minpres.npy and maxpres.npy file, which are simply the upper (maxpres) and lower (minpres) pressure bounds on the 68% and 95% confidence regions, which are derived from the pressures.npy files.

    - NICER_Real_Data/: Followed by the output directories NICER_REAL_ADM_VARYING_BARYONIC and NICER_REAL_BARYONIC, which are the directories which including ADM and neglect it during sampling, respectively. The output files contain the same ouput files as that of the Future-X/ posteriors. 



