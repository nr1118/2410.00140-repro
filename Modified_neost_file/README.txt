The file in this folder is nearly identical to the PosteriorAnalysis.py file used in NEoST v2.1.0. However, since Rutherford et al. 2024b (arXiv:2410.00140) didn't use the extra data associated with compute_auxiliarly_data() and compute_prior_auxliarly_data(), we made a minimal function that only computes the necessary files needed to generate the figures in the manuscript. This function is called compute_minimal_auxliarly_data_ADM() and compute_minimal_auxliarly_data_Baryonic().

Instructions to use the above functions:

1. Replace the PosteriorAnalysis.py file in your NEoST v2.1.0 download with the one in this folder. 
2. Then activate the neost conda environment and make sure you are in the main directory
3. run pip install .
