#!/usr/bin/env python3


import MS_config as cfg
from MS_create_synthetic import *
from MS_run_and_evaluate import *
from os.path import isfile, isdir                       # OS-level utilities
from os import mkdir                                    # OS-level utilities
import shutil                                           # recursive directory removal


ttt = open('running_times-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool), 'a')
xxx = open('stdout-{}.txt'.format(cfg.tool), 'a')
yyy = open('stderr-{}.txt'.format(cfg.tool), 'a')


code_name = 'set6'                                      # fitting using empirical signatures weights, use all COSMICv3 as a reference
for cancer_type in ['Head-SCC', 'ColoRect-AdenoCA', 'Lung-AdenoCA', 'Skin-Melanoma', 'CNS-GBM', 'Stomach-AdenoCA', 'Liver-HCC', 'Lymph-BNHL']:
    rng = np.random.default_rng(0)                      # reset the RNG
    empirical = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')                                # load empirical signature distribution for this cancer type
    for rep in range(cfg.num_realizations):
        print('\n{}: starting run {} for {} & {}'.format(cfg.tool, rep, code_name, cancer_type))
        who = rng.choice(empirical.shape[0], size = cfg.N_samples, replace = True)  # choose N_samples from all samples in the empirical signature distribution (with repetition)
        for num_muts in cfg.num_muts_list_short:        # gradually increase the number of signatures
            contribs = generate_weights_empirical(num_muts, empirical.iloc[who].astype('double').reset_index(drop = True))
            if contribs is not None:
                info_label = '{}\t{}\tw_{}\t{}'.format(cfg.tool, code_name, rep, num_muts)
                prepare_data(rng = rng, num_muts = num_muts, contribs = contribs)   # prepare synthetic mutational catalogs
                true_res = num_muts * contribs                                      # data frame with true signature contributions
                timeout_run(info_label, ttt, xxx, yyy, extra_col = cancer_type)     # run the fitting tool defined in variable tool in MS_config.py
                evaluate_main(info_label, true_res.T, num_muts, extra_col = cancer_type)    # evaluate the estimated signature weights
    system('zip signature_results-{}-{}-{}-{}.zip signature_results/contribution-*.lzma'.format(cfg.WGS_or_WES, cfg.tool, code_name, cancer_type))  # prepare a zip file with compressed (lzma) estimated signature weights for all cohorts
    shutil.rmtree('signature_results')  # remove all result files


code_name = 'set12'                                     # fitting using empirical signatures weights, use only the signatures that appear enough in the results of set6
stats = []                                              # to save the number of reference signatures in individual runs
for cancer_type in ['Head-SCC', 'ColoRect-AdenoCA', 'Lung-AdenoCA', 'Skin-Melanoma', 'CNS-GBM', 'Stomach-AdenoCA', 'Liver-HCC', 'Lymph-BNHL']:
    rng = np.random.default_rng(0)                      # reset the RNG
    empirical = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')
    system('unzip signature_results-{}-{}-set6-{}.zip'.format(cfg.WGS_or_WES, cfg.tool, cancer_type))   # unpack previously computed results where all COSMICv3 signatures were used as a reference
    rename('signature_results', 'results_set6')         # the directory with previously computed is renamed to results_set6
    system('lzma -d results_set6/*.lzma')               # decompress all individual result files
    mkdir('signature_results')                          # prepare new directory for results
    for rep in range(cfg.num_realizations):
        print('\n{}: starting run {} for {} & {}'.format(cfg.tool, rep, code_name, cancer_type))
        who = rng.choice(empirical.shape[0], size = cfg.N_samples, replace = True)
        for num_muts in cfg.num_muts_list_short:        # gradually increase the number of signatures
            if num_muts == 100: rel_weight_thr = 0.05   # the threshold for how many samples need to have a signature to keep this signature in the second round of the fitting process; this threshold is chosen differently for different total numbers of mutations (threshold of 0.01 for 100 mutations would mean that one single mutation would be enough to label a signature as active -> we use a higher relative threshold when the number of mutations is small)
            elif num_muts == 2000: rel_weight_thr = 0.03
            elif num_muts == 50000: rel_weight_thr = 0.01
            else: rel_weight_thr = 0.05
            contribs = generate_weights_empirical(num_muts, empirical.iloc[who].astype('double').reset_index(drop = True))
            try: res_set6 = pd.read_csv('results_set6/contribution-{}-{}-set6-w_{}-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool, rep, num_muts, num_muts), sep = ',', index_col = 0)  # try to load the previously computed results with COSMICv3 as a reference (set6)
            except:
                print('set6 results missing for {}, {}, {}, {}'.format(cfg.WGS_or_WES, cfg.tool, rep, num_muts))
            else:
                sample_sums = res_set6.sum()
                res_set6 = res_set6 / sample_sums       # normalize the results
                frac_sig_active = (res_set6 >= rel_weight_thr).sum(axis = 1) / res_set6.shape[1]    # for each signature, compute the fraction of samples where this signature has the relative weight at least rel_weight_thr
                active_signatures = []                  # prepare the list of active signatures
                for sig in frac_sig_active.index:
                    if frac_sig_active[sig] >= 0.05: active_signatures.append(sig)
                print('{}, realization {}, {} muts: {} active signatures'.format(cancer_type, rep, num_muts, len(active_signatures)))
                stats.append({'cancer_type': cancer_type, 'weights': 'w_{}'.format(rep), '#muts': num_muts, 'sigs_active': len(active_signatures)})
                prepare_relevant_COSMIC(cancer_type, active_sigs = active_signatures)   # prepare the reference list that includes only the sufficiently active signatures
                if contribs is not None:
                    info_label = '{}\t{}\tw_{}\t{}'.format(cfg.tool, code_name, rep, num_muts)
                    prepare_data(rng = rng, num_muts = num_muts, contribs = contribs)
                    true_res = num_muts * contribs
                    timeout_run(info_label, ttt, xxx, yyy, extra_col = cancer_type, which_setup = 'Y')
                    evaluate_main(info_label, true_res.T, num_muts, extra_col = cancer_type)
    system('zip signature_results-{}-{}-{}-{}.zip signature_results/contribution-*.lzma'.format(cfg.WGS_or_WES, cfg.tool, code_name, cancer_type))
    shutil.rmtree('signature_results')
    shutil.rmtree('results_set6')
to_save = pd.DataFrame(stats)   # prepare the data frame with the number of active signatures for each realization and save it
to_save.to_csv('active_signatures-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), sep = '\t', index = False)


ttt.close()
xxx.close()
yyy.close()
