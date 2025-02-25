#!/bin/bash

# -r: use reproduced results instead of published

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
cmd="python plot_routines/cEFT_band.py" # -t
printf "$cmd\n"
$cmd

printf "\nFig. 2\n"
cmd="python plot_routines/MR_priors.py $optional_args" # -r
printf "$cmd\n"
$cmd

printf "\nFig. 3\n"
cmd="python plot_routines/pressure_energydensity_priors.py $optional_args" # -r
printf "$cmd\n"
$cmd

printf "\nFig. 4\n"
cmd="python plot_routines/data_likelihood.py" 
printf "$cmd\n"
$cmd

printf "\nFig. 5\n"
cmd="python plot_routines/MR_baseline_new.py $optional_args" # -r 
$cmd

printf "\nFig. 6\n"
cmd="python plot_routines/pressure_energydensity_posteriors.py $optional_args" # -r 
printf "$cmd\n"
$cmd

printf "\nFig. 7\n"
cmd="python plot_routines/MR_new_heatmap.py $optional_args" # -r 
printf "$cmd\n"
$cmd

printf "\nFig. 8\n"
cmd="python plot_routines/pressure_histograms.py $optional_args" # -r 
printf "$cmd\n"
$cmd

