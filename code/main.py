#!/usr/bin/env python3


from MS_functions import *
from MS_run_and_evaluate import *
import glob


"""
This is the main file of the SigFitTest package.
It contains descriptions of the main functions and examples for how they can be used.
You can uncomment a line of your choice to reproduce results shown in articles
``A comprehensive comparison of tools for fitting mutational signatures'' and
``Inconsistency of signature estimates in real samples with clonal structure''.
You can also use the provided examples to write your own function calls.

matus.medo@unifr.ch, 2025
"""


print('main.py: uncomment a line with a function call or add your own function call')



#####################################################################################################################
# calls related to the benchmarking paper ``A comprehensive comparison of tools for fitting mutational signatures'' #
#####################################################################################################################


# # generate mutational catalogs for the provided cancer types (see folder 'cosmic tissue data' for other cancer types)
# # when out_of_reference_weights are specified, randomly chosen COSMICv3.3.1 signatures that are absent in COSMICv3
# # are assigned these weights
# # output:
# # mutational catalogs saved in the folder '../generated_data'
# # sample calls:
# generate_synthetic_catalogs(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'])
#
# generate_synthetic_catalogs(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'], out_of_reference_weights = [0.15, 0.05])


# # generate simple mutational catalogs with provided signature weights (all samples have the same signature composition),
# # use the fitting tool set in the variable 'tool' in MS_config.py, and evaluate the estimated signature weights
# # output:
# # estimated signature weights ('signature_results-*.zip) and evaluation results (results-*.dat); these files are
# # saved in the main directory, one level up from main.py
# # sample calls:
# fit_with_cosmic3_synthetic_simple(sig_weights = {'SBS3': 1}, code_name = 'SET1')
#
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
#
# fit_with_cosmic3_synthetic(cancer_types = ['Head-SCC', 'ColoRect-AdenoCA', 'Lung-AdenoCA', 'Skin-Melanoma', 'CNS-GBM', 'Stomach-AdenoCA', 'Liver-HCC', 'Lymph-BNHL'], out_of_reference_weights = [0.1, 0.1], code_name = 'SET6')


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
#
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
#
# fit_external(input_catalog = '../sample_data.dat', catalog_GT = '../sample_true_weights.dat', code_name = 'EXT2',
#              reference_signatures = ['SBS1', 'SBS2', 'SBS5', 'SBS13', 'SBS16', 'SBS18'])


# # analyze real samples and evaluate the differences between results obtained by different tools
# differences_real_samples()



#############################################################################################################
# calls related to the paper ``Inconsistency of signature estimates in real samples with clonal structure'' #
#############################################################################################################


# # create synthetic mutational catalogs with clonal samples and bulk samples that are their sums
# # the synthetic catalogs are then fitted and the inconsistency between bulk and clonal samples
# # is computed; when error_mag > 0, the clones are not determined perfectly
# # output:
# # inconsistency values for all analyzed samples (inconsistency_values-%tool%-%code_name%.dat)
# # and summary metrics for signature fitting (results-%tool%-%code_name%.dat)
# # sample calls:
# for muts in [500, 2000, 8000]:
#     for tot_out_of_reference in [0, 0.1, 0.2, 0.3]:
#         fit_synthetic_clones(code_name = '2_simple_clones', clone_sizes = {1: muts, 2: muts},
#                              sig_weights = {1: {'SBS1': 0.5, 'SBS7a': 0.3, 'SBS7c': 0.2},
#                                             2: {'SBS1': 0.3, 'SBS2': 0.4, 'SBS13': 0.3}},
#                              tot_out_of_reference = tot_out_of_reference)
#
# for muts in [500, 2000, 8000]:
#     for tot_out_of_reference in np.linspace(0, 0.4, 9):
#         for clone_error_mag in np.linspace(0, 0.4, 9):
#             fit_synthetic_clones('imperfect_clones_1~error~{:.0f}%'.format(100 * clone_error_mag),
#                                  clone_sizes = {1: muts, 2: muts},
#                                  sig_weights = {1: {'SBS7a': 0.3, 'SBS7c': 0.7},
#                                                 2: {'SBS2': 0.6, 'SBS13': 0.4}},
#                                  # sig_weights = {1: {'SBS1': 0.5, 'SBS7a': 0.3, 'SBS7c': 0.2},
#                                  #                2: {'SBS1': 0.3, 'SBS2': 0.4, 'SBS13': 0.3}},
#                                  tot_out_of_reference = tot_out_of_reference,
#                                  clone_error_mag = clone_error_mag)
#
# for tot_out_of_reference in [0, 0.1, 0.2, 0.3]:
#     fit_synthetic_clones('like_real', clone_sizes = {1: 2846, 2: 2622, 3: 11710},
#                          sig_weights = {1: {'SBS16': 0.21, 'SBS18': 0.37, 'SBS21': 0.11, 'SBS46': 0.31},
#                                         2: {'SBS5': 0.16, 'SBS16': 0.12, 'SBS18': 0.44, 'SBS21': 0.10, 'SBS46': 0.18},
#                                         3: {'SBS3': 0.59, 'SBS5': 0.19, 'SBS18': 0.22}},
#                          tot_out_of_reference = tot_out_of_reference)


# # generate synthetic catalogs where samples differ in their mutation burden and the activity of
# # out-of-reference signatures; unless specified otherwise using parameters muts and weights_OOR,
# # the mutation burdens are 500, 1000, 2000, 4000, 8000 and the out-of-reference signature
# # activities are 0, 0.1, 0.2, 0.3 (20 different combinations in total)
# # output:
# # mutational catalogs data_for_*.dat and true signature activities true_weights_*.dat in folder
# # ../generated_data
# # sample calls:
# generate_mixed_synthetic_catalogs(num_samples = 100, cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'],
#                                   code_name = 'MIXED', rng_seed = 0)
#
generate_mixed_synthetic_catalogs(num_samples = 60, cancer_types = ['Head-SCC', 'ColoRect-AdenoCA'],
                                  code_name = 'MIXED', muts = [1000, 10000], weights_OOR = [0, 0.2])


# # fit samples in the provided mutatational catalog using the fitting tool set in the variable 'tool'
# # in MS_config.py; in the second step, the quality of the fits is assessed by fitting synthetic data
# # using the same tool; the catalog at the provided path (input_catalog) can be prepared using functions
# # generate_synthetic_catalogs() or generate_mixed_synthetic_catalogs() mentioned above, or they can have
# # an external origin
# # output:
# # 1) estimated signature weights (*_contribution.dat) and fit assessment results (fit_assessment-*.dat)
# #    are saved in the main folder
# # 2) if true signature activities (ground truth, GT) are provided, then fitting performance metrics are
# #    saved for all samples (results_individual_samples-*.dat) and for the cohort (results-*.dat)
# # sample calls:
for ct in ['Head-SCC', 'ColoRect-AdenoCA']:
    fit_and_assess(input_catalog = '../generated_data/data_for_deconstructSigs_MIXED_{}.dat'.format(ct),
                   GT = '../generated_data/true_weights_MIXED_{}.dat'.format(ct),
                   code_name = 'MIXED_{}'.format(ct))
