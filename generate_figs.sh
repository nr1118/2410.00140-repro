#!/bin/bash

# -r: use reproduced results instead of published
# -name_prior: run name for the reproduced prior results instead of published
# -name_posterior_incladm_real: run name for the reproduced real data including ADM results instead of published
# -name_posterior_negladm_real: run name for the reproduced real data neglecting ADM results instead of published
# -name_posterior_incladm_adm: run name for the reproduced Future-X ADM core model including ADM results instead of published
# -name_posterior_negladm_adm: run name for the reproduced Future-X ADM core model neglecting ADM results instead of published
# -name_posterior_incladm_noadm: run name for the reproduced Future-X No ADM model including ADM results instead of published
# -name_posterior_negladm_noadm: run name for the reproduced Future-X No ADM model neglecting ADM results instead of published
#-name_appendix: run name for the reproduced Appendix plots of approximating zero self-repulsion
optional_args=$1

if [ "$optional_args" == "-r" ]
then
	printf "Using reproduced results.\n"
else
	printf "Using published results.\n"
	printf "Use flag -r with this script to plot reproduced (as opposed to published) results.\n"
fi

mkdir -p plots

printf "\nFig. 1\n"
cmd="python plot_routines/Figure_1.py" # -r / -name_prior
printf "$cmd\n"
$cmd

printf "\nFig. 2\n"
cmd="python plot_routines/Figure_2.py"
printf "$cmd\n"
$cmd

printf "\nFig. 3\n"
cmd="python plot_routines/Figure_3.py $optional_args" # -r / -name_prior / -name_posterior_incladm_real
printf "$cmd\n"
$cmd

printf "\nFig. 4\n"
cmd="python plot_routines/Figure_4.py $optional_args" # -r / -name_prior / -name_posterior_incladm_real / -name_posterior_negladm_real 
printf "$cmd\n"
$cmd

printf "\nFig. 5\n"
cmd="python plot_routines/Figure_5.py $optional_args" # -r / -name_prior / -name_posterior_incladm_adm / -name_posterior_incladm_noadm
$cmd

printf "\nFig. 6\n"
cmd="python plot_routines/Figure_6.py $optional_args" # -r / -name_prior / -name_posterior_incladm_adm / -name_posterior_negladm_adm / -name_posterior_incladm_noadm / -name_posterior_negladm_noadm 
printf "$cmd\n"
$cmd

printf "\nFig. 7\n"
cmd="python plot_routines/Figures_7_and_8.py $optional_args" # -r / -name_prior / -name_appendix
printf "$cmd\n"
$cmd



