#!/usr/bin/env python3


from MS_functions import *
from MS_run_and_evaluate import *


"""
This is the main file of the SigFitTest package.
It contains descriptions of the main functions and examples for how they can be used.
Uncomment a line of your choice to reproduce results shown in the article
``A comprehensive comparison of tools for fitting mutational signatures'' or
write your own function calls.

matus.medo@unifr.ch, 2025
"""


# # generate mutational catalogs for the provided cancer types (see directory 'cosmic tissue data' for further cancer types)
# # output:
# # mutational catalogs saved in the directory 'data'
# # sample calls:
# generate_synthetic_catalogs(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'])
# generate_synthetic_catalogs(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'], out_of_reference_weights = [0.15, 0.05])


# # generate simple mutational catalogs with provided signature weights (all samples have the same signature composition),
# # use the fitting tool set in the variable 'tool' in MS_config.py, and evaluate the estimated signature weights
# # output:
# # estimated signature weights ('signature_results-*.zip) and evaluation results (results-*.dat); these files are
# # saved in the main directory, one level up from main.py
# # sample calls:
# fit_with_cosmic3_synthetic_simple(sig_weights = {'SBS3': 1}, code_name = 'SET1')
# fit_with_cosmic3_synthetic_simple(sig_weights = {'SBS1': 0.7, 'SBS5': 0.3}, code_name = 'SET2')


# # generate mutational catalogs for the provided cancer types, use the fitting tool set in the variable 'tool' in
# # MS_config.py, and evaluate the estimated signature weights
# # when out_of_reference_weights are specified, randomly chosen COSMICv3.3.1 signatures that are absent in COSMICv3
# # are assigned these weights
# # output:
# # estimated signature weights ('signature_results-*.zip), evaluation results (results-*.dat), and the running
# # time (running_times-*.dat); these files are saved in the main directory, one level up from main.py
# # sample calls:
# fit_with_cosmic3_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'], code_name = 'SET3')
# fit_with_cosmic3_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA', 'Lung-AdenoCA', 'Skin-Melanoma', 'CNS-GBM', 'Stomach-AdenoCA', 'Liver-HCC', 'Lymph-BNHL'], out_of_reference_weights = [0.1, 0.1], code_name = 'SET6')
# fit_with_cosmic3_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'], code_name = 'SET3', evaluate_fit_quality = True)
# fit_with_cosmic3_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA', 'Lung-AdenoCA', 'Skin-Melanoma', 'CNS-GBM', 'Stomach-AdenoCA', 'Liver-HCC', 'Lymph-BNHL'], out_of_reference_weights = [0.1, 0.1], code_name = 'SET6', evaluate_fit_quality = True)


generate_mixed_synthetic_catalogs(num_samples = 100, cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'], rng_seed = 0, code_name = 'MIXED')
for ct in ['Head-SCC', 'ColoRect-AdenoCA']:
    fit_and_assess(input_catalog = '../generated_data/data_for_deconstructSigs_MIXED_{}.dat'.format(ct),
                   code_name = 'MIXED_{}'.format(ct))


# # generate mutational catalogs for the provided cancer types, load the previously computed results of
# # fit_with_cosmic3_synthetic ("SET3"), use them to find which signatures are sufficiently active, fit
# # the catalogs using the active signatures as a reference, and evaluate the estimated signature weights
# # output:
# # estimated signature weights ('signature_results-*.zip), evaluation results (results-*.dat), and
# # information about the chosen signatures (chosen_signatures-*.dat); these files are saved in the
# # main directory, one level up from main.py
# # sample calls:
# prune_reference_and_fit_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'], code_name = 'SET4')


# # generate mutational catalogs for the provided cancer types, introduce systematic differences in signature
# # weights between odd and even samples, use the fitting tool set in the variable 'tool' in MS_config.py,
# # and test whether the estimated signature weights differ significantly or not
# # when out_of_reference_weights are specified, randomly chosen COSMICv3.3.1 signatures that are absent in
# # COSMICv3 are assigned these weights
# # output:
# # same as fit_with_cosmic3_synthetic() plus comparison_results-*.dat which contains the results of the
# # Wilcoxon rank-sum test for true as well as estimated signature weights
# # sample calls:
# fit_with_cosmic3_synthetic_compare_groups(cancer_types = ['CNS-GBM'], which_sig = 'SBS40', difference_magnitude = 0.3, code_name = 'SET5A')
# fit_with_cosmic3_synthetic_compare_groups(cancer_types = ['CNS-GBM'], which_sig = 'SBS1', difference_magnitude = 0.2, code_name = 'SET5B')


# # generate mutational catalogs by subsampling from real mutational catalogs, use the fitting tool set in the
# # variable 'tool' in MS_config.py, evaluate the estimated signature weights by comparing with a pre-computed ground
# # truth file
# # output:
# # estimated signature weights ('signature_results-*.zip), evaluation results (results-*.dat), and evaluation
# # results for individual samples (individual_samples-*.dat); these files are saved in the main directory, one level
# # up from main.py
# # sample calls:
# fit_with_cosmic3_subsampled_real_catalogs(code_name = 'SET7')


# # use the fitting tool set in the variable 'tool' in MS_config.py to fit a provided input catalog, evaluate
# # the estimated signature weights against the provided catalog ground truth (catalog_GT); a subset of COSMICv3.3
# # can be provided to constrain the reference signature catalog
# # output:
# # estimated signature activities (contribution-*.dat.lzma) and performance metrics (results-*.dat)
# # sample calls:
# fit_external(input_catalog = '../sample_data.dat', catalog_GT = '../sample_true_weights.dat', code_name = 'EXT1')
# fit_external(input_catalog = '../sample_data.dat', catalog_GT = '../sample_true_weights.dat', code_name = 'EXT2', reference_signatures = ['SBS1', 'SBS2', 'SBS5', 'SBS13', 'SBS16', 'SBS18'])


# # analyze real samples and evaluate the differences between results obtained by different tools
# differences_real_samples()


# # assess the fit quality for a given set of signature estimates and the corresponding mutational catalog
# import glob
# result_files = glob.glob('../*-contribution.dat')
# for fname in result_files:
#     assess_fit_quality(input_catalog = '../samples_Sasha_Blay-num_muts_threshold_100.dat', sig_estimates = fname, ref_genome = 'GRCh37')
# assess_fit_quality(info_label = 'SPA-COSMICv34', input_catalog = '../mutation_catalogs_CBTN_bulk_and_clones-num_muts_threshold_100.dat', sig_estimates = '../SPA-COSMICv34-contribution.dat', ref_sigs = 'COSMIC_v3.4', ref_genome = 'GRCh38')


# fit synthetic samples & assess the fit quality
# fit_with_cosmic3_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA', 'Lung-AdenoCA', 'Skin-Melanoma', 'CNS-GBM', 'Stomach-AdenoCA', 'Liver-HCC', 'Lymph-BNHL'], out_of_reference_weights = {0.5: [0.2]}, code_name = 'FIND1', evaluate_fit_quality = True)


# # clones & bulk synthetic catalogs, modeled after a specific patient sample
# for tot_out_of_reference in [0, 0.1, 0.2, 0.3]:
#     fit_synthetic_clones('like_P48B', clone_sizes = {1: 2846, 2: 2622, 3: 11710},
#                          sig_weights = {1: {'SBS16': 0.21, 'SBS18': 0.37, 'SBS21': 0.11, 'SBS46': 0.31},
#                                         2: {'SBS5': 0.16, 'SBS16': 0.12, 'SBS18': 0.44, 'SBS21': 0.10, 'SBS46': 0.18},
#                                         3: {'SBS3': 0.59, 'SBS5': 0.19, 'SBS18': 0.22}},
#                          tot_out_of_reference = tot_out_of_reference)


#! generalize fit_synthetic_clones to take sig_weights values as integers---the corresponding number of random signatures will be chosen for each sample then
# clones & bulk synthetic catalogs, simple signature activities, gradually increasing the number of mutations
# for muts in [500, 2000, 8000]:
#     for tot_out_of_reference in [0, 0.1, 0.2, 0.3]:
#         fit_synthetic_clones('simple_clones2', clone_sizes = {1: muts, 2: muts},
#                              sig_weights = {1: {'SBS1': 0.5, 'SBS7a': 0.3, 'SBS7c': 0.2},
#                                             2: {'SBS1': 0.3, 'SBS2': 0.4, 'SBS13': 0.3}},
#                              tot_out_of_reference = tot_out_of_reference)


# # check what happens when the clones are not determined perfectly
# for muts in [500, 2000, 8000]:
#     for tot_out_of_reference in np.linspace(0, 0.4, 9):
#         for clone_error_mag in np.linspace(0, 0.4, 9):
#             fit_synthetic_clones('imperfect_clones_1~error~{:.0f}%'.format(100 * clone_error_mag),
#                                  clone_sizes = {1: muts, 2: muts},
#                                  sig_weights = {1: {'SBS7a': 0.3, 'SBS7c': 0.7},
#                                                 2: {'SBS2': 0.6, 'SBS13': 0.4}},
#                                  # sig_weights = {1: {'SBS1': 0.5, 'SBS7a': 0.3, 'SBS7c': 0.2},
#                                  #                2: {'SBS1': 0.3, 'SBS2': 0.4, 'SBS13': 0.3}},
#                                  # sig_weights = {1: {'SBS1': 0.48, 'SBS2': 0.04, 'SBS13': 0.03, 'SBS7a': 0.27, 'SBS7c': 0.18},
#                                  #                2: {'SBS1': 0.32, 'SBS2': 0.36, 'SBS13': 0.27, 'SBS7a': 0.03, 'SBS7c': 0.02}},
#                                  tot_out_of_reference = tot_out_of_reference,
#                                  clone_error_mag = clone_error_mag)
