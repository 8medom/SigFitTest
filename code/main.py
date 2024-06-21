#!/usr/bin/env python3


from MS_functions import *


"""
This is the main function of the SigFitTest package.
It contains descriptions of the main functions and examples for how they can be used.
Uncomment a line of your choice to reproduce results shown in the article
``A comprehensive comparison of tools for fitting mutational signatures'' or
write your own function calls.

matus.medo@unifr.ch, 2024
"""


# # generate mutational catalogs for the provided cancer types (see directory 'cosmic tissue data' for further cancer types)
# # output: mutational catalogs saved in the directory 'data'
# # example calls:
# generate_synthetic_catalogs(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'])
# generate_synthetic_catalogs(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'], out_of_reference_weights = [0.15, 0.05])


# # generate simple mutational catalogs with provided signature weights (all samples have the same signature composition),
# # use the fitting tool set in the variable 'tool' in MS_config.py, and evaluate the estimated signature weights
# # output: estimated signature weights ('signature_results-*.zip) and evaluation results (results-*.dat); these files are
# # saved in the main directory, one level up from the directory code where main.py is located
# # example calls:
# fit_with_cosmic3_synthetic_simple(sig_weights = {'SBS3': 1}, code_name = 'SET1')
# fit_with_cosmic3_synthetic_simple(sig_weights = {'SBS1': 0.7, 'SBS5': 0.3}, code_name = 'SET2')


# # generate mutational catalogs for the provided cancer types, use the fitting tool set in the variable 'tool' in
# # MS_config.py, and evaluate the estimated signature weights
# # when out_of_reference_weights are specified, randomly chosen COSMICv3.3.1 signatures that are absent in COSMICv3
# # are assigned these weights
# # output: estimated signature weights ('signature_results-*.zip) and evaluation results (results-*.dat); these files
# # are saved in the main directory, one level up from the directory code where main.py is located
# # example calls:
fit_with_cosmic3_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'], code_name = 'SET3')
# fit_with_cosmic3_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA', 'Lung-AdenoCA', 'Skin-Melanoma', 'CNS-GBM', 'Stomach-AdenoCA', 'Liver-HCC', 'Lymph-BNHL'], out_of_reference_weights = [0.1, 0.1], code_name = 'SET6')


# # generate mutational catalogs for the provided cancer types, load the previously computed results of
# # fit_with_cosmic3_synthetic ("set6"), use them to find which signatures are sufficiently active, fit the catalogs
# # using the active signatures as a reference, and evaluate the estimated signature weights
# # output: estimated signature weights ('signature_results-*.zip) and evaluation results (results-*.dat); these files
# # are saved in the main directory, one level up from the directory code where main.py is located
# # example calls:
# prune_reference_and_fit_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'], code_name = 'SET4')


# # generate mutational catalogs for the provided cancer types, introduce systematic differences in signature
# # weights between odd and even samples, use the fitting tool set in the variable 'tool' in MS_config.py,
# # and test whether the estimated signature weights differ significantly or not
# # when out_of_reference_weights are specified, randomly chosen COSMICv3.3.1 signatures that are absent in COSMICv3
# # are assigned these weights
# # output: same as fit_with_cosmic3_synthetic() plus comparison_results-*.dat which contains the results of
# # the Wilcoxon rank-sum test for true as well as estimated signature weights
# # example calls:
# fit_with_cosmic3_synthetic_compare_groups(cancer_types = ['CNS-GBM'], which_sig = 'SBS40', difference_magnitude = 0.3, code_name = 'SET5A')
# fit_with_cosmic3_synthetic_compare_groups(cancer_types = ['CNS-GBM'], which_sig = 'SBS1', difference_magnitude = 0.2, code_name = 'SET5B')


# # generate mutational catalogs by subsampling from real mutational catalogs, use the fitting tool set in the
# # variable 'tool' in MS_config.py, evaluate the estimated signature weights by comparing with a pre-computed ground
# # truth file
# # output: estimated signature weights ('signature_results-*.zip), evaluation results (results-*.dat), and evaluation
# # results for individual samples (individual_samples-*.dat); these files are saved in the main directory, one level
# # up from the directory code where main.py is located
# # example calls:
# fit_with_cosmic3_subsampled_real_catalogs(code_name = 'SET7')
