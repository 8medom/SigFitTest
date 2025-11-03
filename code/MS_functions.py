import MS_config as cfg
from MS_create_synthetic import *
from MS_run_and_evaluate import *
from os.path import isfile, isdir                       # OS-level utilities
from os import mkdir, system                            # OS-level utilities
import shutil                                           # recursive directory removal
import sys                                              # emergency stop
from time import sleep                                  # to be able to pause when needed


# check if the weights are okay
def check_provided_weights(weights):
    if len(weights) > len(cfg.out_sigs):
        print('cannot introduce {} out-of-reference signatures (only {} signatures are in COSMICv3.3.1 but not in COSMICv3)'.format(len(weights), len(cfg.out_sigs)))
        sys.exit(1)
    if sum(weights) > 1 + cfg.EPSILON:
        print('sum of out-of-reference weights cannot exceed 1 (now {:.2f})'.format(sum(weights)))
        sys.exit(1)
    for weight in weights:
        if weight < 0:
            print('out-of-reference weights cannot be negative ({})'.format(weight))
            sys.exit(1)


# generate and save synthetic mutational catalogs with signature weights driven by COSMIC results on real tissue data
# out-of-reference signatures are assigned the weights given in out_of_reference_weights (default: no out-of-reference signatures)
def generate_synthetic_catalogs(cancer_types, out_of_reference_weights = []):
    tot_out_of_reference = sum(out_of_reference_weights)
    if tot_out_of_reference > 0: print('=== generating synthetic mutational catalogs based on real tissue data; out-of-reference signatures have joint weight {:.4f} ==='.format(tot_out_of_reference))
    else: print('=== generating synthetic mutational catalogs based on real tissue data ===')
    check_provided_weights(out_of_reference_weights)
    if not isdir('data'): mkdir('data')
    for cancer_type in cancer_types:
        rng = np.random.default_rng(0)                  # initialize the RNG
        empirical = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')                                # load empirical signature distribution for this cancer type
        for rep in range(cfg.num_realizations):         # run the chosen number of realizations
            print('\ngenerating catalog #{} for {}'.format(rep, cancer_type))
            who = rng.choice(empirical.shape[0], size = cfg.N_samples, replace = True)  # choose N_samples from all samples in the empirical signature distribution (with repetition)
            if tot_out_of_reference > 0:                # if needed, choose which out-of-reference signatures to add
                new_active = rng.choice(cfg.out_sigs, size = len(out_of_reference_weights), replace = False)
            for num_muts in cfg.num_muts_list_short:    # gradually increase the number of mutations
                print('{} mutations...'.format(num_muts), end = ' ', flush = True)
                contribs = generate_weights_empirical(num_muts, empirical.iloc[who].astype('double').reset_index(drop = True))
                if contribs is not None:
                    if tot_out_of_reference > 0:
                        contribs *= (1 - tot_out_of_reference)
                        for n, weight in enumerate(out_of_reference_weights):
                            contribs[new_active[n]] = weight
                    info_label = '{}-w_{}-{}'.format(cancer_type, rep, num_muts)
                    # generate and save synthetic mutational catalogs
                    muts = pd.Series(num_muts, index = contribs.index)
                    counts = prepare_data_from_signature_activity(rng = rng, muts = muts, contribs = contribs)
                    save_catalogs(counts = counts, info_label = info_label)
                    # keep only the active signatures
                    tmp = contribs[contribs.columns[(contribs.sum(axis = 0) != 0)]]
                    # save true signature contributions
                    tmp.to_csv('data/true_weights-{}.dat'.format(info_label), sep = '\t', float_format = '%.6f')
    shutil.move('data', '../generated_data')
    print('\n\ndone, datasets saved in folder generated_data')


# fitting syntetic mutational catalogs with identical signatures activity for all samples (specified by sig_weights)
# examples: sig_weights = {'SBS3': 1}; sig_weights = {'SBS1': 0.7, 'SBS5': 0.3}
# all COSMICv3 signatures are used as a reference
def fit_with_cosmic3_synthetic_simple(sig_weights, code_name):
    print('=== fitting simple synthetic mutational catalogs using COSMICv3 ===')
    print('active signatures: {}'.format(', '.join(['{} ({})'.format(sig, sig_weights[sig]) for sig in sig_weights.keys()])))
    ttt = open('../running_times-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool), 'a')
    xxx = open('../stdout-{}.txt'.format(cfg.tool), 'a')
    yyy = open('../stderr-{}.txt'.format(cfg.tool), 'a')
    if not isdir('data'): mkdir('data')
    if not isdir('signature_results'): mkdir('signature_results')
    rng = np.random.default_rng(0)                              # initialize the RNG
    print(cfg.header_line)
    for num_muts in cfg.num_muts_list:                          # to gradually increase the number of mutations
        contribs = pd.DataFrame(0, index = ['S{}'.format(x) for x in range(cfg.N_samples)], columns = cfg.input_sigs.columns, dtype = float)
        for sig in sig_weights: contribs[sig] = sig_weights[sig]
        sig_info = '~'.join(['{}({})'.format(sig, sig_weights[sig]) for sig in sig_weights.keys()])
        info_label = '{}\t{}\t{}\t{}'.format(cfg.tool, code_name, sig_info, num_muts)
        muts = pd.Series(num_muts, index = contribs.index)      # prepare synthetic mutational catalogs
        counts = prepare_data_from_signature_activity(rng = rng, muts = muts, contribs = contribs)
        save_catalogs(counts = counts)
        true_res = (contribs.T * muts).T                        # data frame with true signature contributions
        timeout_run(info_label, ttt, xxx, yyy)                  # run the fitting tool defined in variable tool in MS_config.py
        evaluate_main(info_label, true_res.T, muts)             # evaluate the estimated signature weights
    zip_oname = '../signature_results-{}-{}-{}.zip'.format(cfg.WGS_or_WES, cfg.tool, code_name)
    system('zip {} signature_results/contribution-*.lzma'.format(zip_oname))    # zip the estimated signature weights for all cohorts
    shutil.rmtree('signature_results')                          # remove all result files
    shutil.rmtree('data')                                       # remove the directory with synthetic data
    ttt.close()
    xxx.close()
    yyy.close()


# fitting synthetic mutational catalogs with empirical signatures weights, using all COSMICv3 signatures as a reference
# out-of-reference signatures are assigned the weights given in out_of_reference_weights (default: no out-of-reference signatures)
# examples: out_of_reference_weights = [0.2, 0.1] means that one out of reference signature has weight 20% and another has weight 10%
# out-of-reference signatures are the same for all samples in a cohort
def fit_with_cosmic3_synthetic(cancer_types, code_name, out_of_reference_weights = [], save_true_weights = False, evaluate_fit_quality = False):
    tot_out_of_reference = sum(out_of_reference_weights)
    if tot_out_of_reference > 0: print('=== fitting synthetic mutational catalogs using COSMICv3; out-of-reference signatures have joint weight {:.4f} ==='.format(tot_out_of_reference))
    else: print('=== fitting synthetic mutational catalogs using COSMICv3 ===')
    check_provided_weights(out_of_reference_weights)
    ttt = open('../running_times-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool), 'a')
    xxx = open('../stdout-{}.txt'.format(cfg.tool), 'a')
    yyy = open('../stderr-{}.txt'.format(cfg.tool), 'a')
    for cancer_type in cancer_types:
        if not isdir('data'): mkdir('data')
        if not isdir('signature_results'): mkdir('signature_results')
        rng = np.random.default_rng(0)                              # initialize the RNG
        empirical = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')                                # load empirical signature distribution for this cancer type
        for rep in range(cfg.num_realizations):
            print('\n{}: starting run {} for {} & {}'.format(cfg.tool, rep, code_name, cancer_type))
            print('cancer_type\t{}'.format(cfg.header_line))
            # choose N_samples from all samples in the empirical signature distribution (with repetition)
            who = rng.choice(empirical.shape[0], size = cfg.N_samples, replace = True)
            for num_muts in cfg.num_muts_list_short:                # to gradually increase the number of mutations
                contribs = generate_weights_empirical(num_muts, empirical.iloc[who].astype('double').reset_index(drop = True))
                if contribs is not None:
                    info_label = '{}\t{}\tw_{}\t{}'.format(cfg.tool, code_name, rep, num_muts)
                    if tot_out_of_reference > 0:
                        contribs *= (1 - tot_out_of_reference)
                        new_active = rng.choice(cfg.out_sigs, size = len(out_of_reference_weights), replace = False)
                        for n, weight in enumerate(out_of_reference_weights):
                            contribs[new_active[n]] = weight
                    if save_true_weights:                           # save true signature contributions
                        tmp = contribs[contribs.columns[(contribs.sum(axis = 0) != 0)]]
                        tmp.to_csv('data/true_weights_{}_{}_w{}_{}.dat'.format(cancer_type, code_name, rep, num_muts), sep = '\t', float_format = '%.6f')
                    # prepare synthetic mutational catalogs
                    muts = pd.Series(num_muts, index = contribs.index)
                    counts = prepare_data_from_signature_activity(rng = rng, muts = muts, contribs = contribs)
                    save_catalogs(counts = counts)
                    true_res = (contribs.T * muts).T                # data frame with true signature contributions
                    timeout_run(info_label, ttt, xxx, yyy, extra_col = cancer_type)         # fit using the tool defined in MS_config.py
                    if evaluate_fit_quality: evaluate_fits(info_label, counts, extra_col = cancer_type)
                    evaluate_main(info_label, true_res.T, muts, extra_col = cancer_type)
        zip_oname = '../signature_results-{}-{}-{}-{}.zip'.format(cfg.WGS_or_WES, cfg.tool, code_name, cancer_type)
        system('zip {} signature_results/contribution-*.lzma'.format(zip_oname))            # zip signature activity estimates
        if cfg.tool == 'sigfit':                                    # if running sifgit, save also its lower activity estimates
            system('zip ../lower90-{}-{}-{}-{}.zip signature_results/lower90-*.lzma'.format(cfg.WGS_or_WES, cfg.tool, code_name, cancer_type))
        shutil.rmtree('signature_results')                          # remove all result files
        if save_true_weights:
            system('zip ../true_signature_weights-{}-{}.zip data/true_weights*.dat'.format(cancer_type, code_name))
        shutil.rmtree('data')                                       # remove the directory with synthetic data
    ttt.close()
    xxx.close()
    yyy.close()


# fitting synthetic mutational catalogs with empirical signatures weights, using only the signatures that are sufficiently active
# results of SET3 (fitting against all COSMICv3) for the corresponding cancer type-number of mutations need to be in the code directory
def prune_reference_and_fit_synthetic(cancer_types, code_name):
    print('=== using previously-computed fitting results to select reference signatures and fit again ===')
    ttt = open('../running_times-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool), 'a')
    xxx = open('../stdout-{}.txt'.format(cfg.tool), 'a')
    yyy = open('../stderr-{}.txt'.format(cfg.tool), 'a')
    stats = []                                              # to save the number of reference signatures in individual runs
    for cancer_type in cancer_types:
        if not isdir('data'): mkdir('data')
        if not isdir('signature_results'): mkdir('signature_results')
        rng = np.random.default_rng(0)                      # initialize the RNG
        empirical = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')
        # unpack previously computed results where all COSMICv3 signatures were used as a reference
        system('unzip signature_results-{}-{}-SET3-{}.zip'.format(cfg.WGS_or_WES, cfg.tool, cancer_type))
        # the directory with previously computed is renamed to results_SET3
        rename('signature_results', 'results_SET3')
        # decompress all individual result files
        system('lzma -d results_SET3/*.lzma')
        # prepare new directory for results
        mkdir('signature_results')
        for rep in range(cfg.num_realizations):
            print('\n{}: starting run {} for {} & {}'.format(cfg.tool, rep, code_name, cancer_type))
            who = rng.choice(empirical.shape[0], size = cfg.N_samples, replace = True)
            for num_muts in cfg.num_muts_list_short:        # gradually increase the number of mutations
                if num_muts == 100: rel_weight_thr = 0.05   # the threshold for how many samples need to have a signature to keep this signature in the second round of the fitting process; this threshold is chosen differently for different total numbers of mutations (threshold of 0.01 for 100 mutations would mean that one single mutation would be enough to label a signature as active -> we use a higher relative threshold when the number of mutations is small)
                elif num_muts == 2000: rel_weight_thr = 0.03
                elif num_muts == 50000: rel_weight_thr = 0.01
                else: rel_weight_thr = 0.05
                contribs = generate_weights_empirical(num_muts, empirical.iloc[who].astype('double').reset_index(drop = True))
                try:    # try to load the previously computed results with COSMICv3 as a reference (SET3)
                    res_set3 = pd.read_csv('results_SET3/contribution-{}-{}-{}-SET3-w_{}-{}.dat'.format(cfg.WGS_or_WES, cancer_type, cfg.tool, rep, num_muts), sep = ',', index_col = 0)
                except:
                    print('SET3 results missing for {}, {}, {}, {}, {}'.format(cfg.WGS_or_WES, cancer_type, cfg.tool, rep, num_muts))
                else:
                    sample_sums = res_set3.sum()
                    res_set3 = res_set3 / sample_sums   # normalize the results
                    # compute the fraction of samples where each signature has the relative weight at least rel_weight_thr
                    frac_sig_active = (res_set3 >= rel_weight_thr).sum(axis = 1) / res_set3.shape[1]
                    active_signatures = []              # prepare the list of active signatures
                    for sig in frac_sig_active.index:
                        if frac_sig_active[sig] >= 0.05: active_signatures.append(sig)
                    print('{}, realization {}, {} muts: {} active signatures'.format(cancer_type, rep, num_muts, len(active_signatures)))
                    stats.append({'cancer_type': cancer_type, 'weights': 'w_{}'.format(rep), '#muts': num_muts, '#sigs_chosen': len(active_signatures), 'sigs_chosen': ','.join(active_signatures)})
                    # prepare the reference list that includes only the sufficiently active signatures
                    prepare_relevant_COSMIC(cancer_type, chosen_sigs = active_signatures)
                    if contribs is not None:
                        info_label = '{}\t{}\tw_{}\t{}'.format(cfg.tool, code_name, rep, num_muts)
                        muts = pd.Series(num_muts, index = contribs.index)
                        counts = prepare_data_from_signature_activity(rng = rng, muts = muts, contribs = contribs)
                        save_catalogs(counts = counts)
                        true_res = (contribs.T * muts).T
                        timeout_run(info_label, ttt, xxx, yyy, extra_col = cancer_type, which_setup = 'Y')
                        print('cancer_type\t{}'.format(cfg.header_line))
                        evaluate_main(info_label, true_res.T, muts, extra_col = cancer_type)
        system('zip ../signature_results-{}-{}-{}-{}.zip signature_results/contribution-*.lzma'.format(cfg.WGS_or_WES, cfg.tool, code_name, cancer_type))
        shutil.rmtree('signature_results')  # remove all result files
        shutil.rmtree('data')               # remove the directory with synthetic data
        shutil.rmtree('results_SET3')       # remove the directory with previously computed results
    to_save = pd.DataFrame(stats)   # prepare the data frame with the number of active signatures for each realization and save it
    to_save.to_csv('../chosen_signatures-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), sep = '\t', index = False)
    ttt.close()
    xxx.close()
    yyy.close()


# fitting synthetic mutational catalogs with empirical signatures weights, using all COSMICv3 signatures as a reference
# out-of-reference signatures are assigned the weights given in out_of_reference_weights (default: no out-of-reference signatures)
# systematic differences are introduced between odd and even samples
# Wilcoxon rank-sum test is used to compare true as well as estimated signature weights between the two groups of samples
def fit_with_cosmic3_synthetic_compare_groups(cancer_types, code_name, which_sig, difference_magnitude, out_of_reference_weights = [], save_true_weights = False):
    tot_out_of_reference = sum(out_of_reference_weights)
    if tot_out_of_reference > 0: print('=== fitting synthetic mutational catalogs using COSMICv3; out-of-reference signatures have joint weight {:.4f} ==='.format(tot_out_of_reference))
    else: print('=== fitting synthetic mutational catalogs using COSMICv3 ===')
    check_provided_weights(out_of_reference_weights)
    ttt = open('../running_times-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool), 'a')
    xxx = open('../stdout-{}.txt'.format(cfg.tool), 'a')
    yyy = open('../stderr-{}.txt'.format(cfg.tool), 'a')
    for cancer_type in cancer_types:
        if not isdir('data'): mkdir('data')
        if not isdir('signature_results'): mkdir('signature_results')
        rng = np.random.default_rng(0)                              # initialize the RNG
        empirical = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')                                # load empirical signature distribution for this cancer type
        for rep in range(cfg.num_realizations):
            for cohort_size in [50, 100, 150, 200]:
                print('\n{}: starting run {} for {} with {} samples'.format(cfg.tool, rep, cancer_type, cohort_size))
                # choose cohort_size samples from all samples in the empirical signature distribution (with repetition)
                who = rng.choice(empirical.shape[0], size = cohort_size, replace = True)
                contribs = generate_weights_empirical(200, empirical.iloc[who].astype('double').reset_index(drop = True), cohort_size = cohort_size)    # all signatures with weight < 0.05 are set to zero (they are weak signatures contributing less than 10 mutations; 200 * 0.05 = 10)
                if tot_out_of_reference > 0:                            # introduce out-of-reference signatures if necessary
                    contribs *= (1 - tot_out_of_reference)
                    new_active = rng.choice(cfg.out_sigs, size = len(out_of_reference_weights), replace = False)
                    for n, weight in enumerate(out_of_reference_weights):
                        contribs[new_active[n]] = weight
                if save_true_weights:                                   # save true signature contributions
                    tmp = contribs[contribs.columns[(contribs.sum(axis = 0) != 0)]]
                    tmp.to_csv('data/true_weights_{}_{}_w{}.dat'.format(cancer_type, code_name, rep), sep = '\t', float_format = '%.6f')
                contribs = introduce_weight_differences(contribs, which_sig, difference_magnitude)
                info_label = '{}\t{}\t{}\tw_{}\t{}'.format('actual weights', code_name, cohort_size, rep, 0)
                significance_value = compare_groups(info_label, which_sig, 0, difference_magnitude, true_res = contribs, extra_col = cancer_type)
                if significance_value < 0.05:                           # if true weights differ significantly, fit mutational catalogs
                    for n, num_muts in enumerate([100, 1000, 10000]):   # gradually increase the number of mutations
                        info_label = '{}\t{}\t{}\tw_{}\t{}'.format(cfg.tool, code_name, cohort_size, rep, num_muts)
                        # prepare synthetic mutational catalogs
                        muts = pd.Series(num_muts, index = contribs.index)
                        counts = prepare_data_from_signature_activity(rng = rng, muts = muts, contribs = contribs)
                        save_catalogs(counts = counts)
                        # data frame with true signature contributions
                        true_res = (contribs.T * muts).T
                        # run the fitting tool defined in variable tool in MS_config.py
                        timeout_run(info_label, ttt, xxx, yyy, extra_col = cancer_type)
                        # evaluate the estimated signature weights
                        evaluate_main(info_label, true_res.T, muts, extra_col = cancer_type, compress_result_file = False)
                        compare_groups(info_label, which_sig, muts, difference_magnitude, extra_col = cancer_type)
                else:
                    print('true weights of {} do not differ significantly -> skipping this realization'.format(which_sig))
        # prepare a zip file with compressed (lzma) estimated signature weights for all cohorts
        system('zip ../signature_results-{}-{}-{}-{}.zip signature_results/contribution-*.lzma'.format(cfg.WGS_or_WES, cfg.tool, code_name, cancer_type))
        if cfg.tool == 'sigfit':
            system('zip ../lower90-{}-{}-{}-{}.zip signature_results/lower90-*.lzma'.format(cfg.WGS_or_WES, cfg.tool, code_name, cancer_type))
        shutil.rmtree('signature_results')  # remove all result files
        if save_true_weights:
            system('zip ../true_signature_weights-{}-{}.zip data/true_weights*.dat'.format(cancer_type, code_name))
        shutil.rmtree('data')               # remove the directory with synthetic data
    ttt.close()
    xxx.close()
    yyy.close()


# fitting subsampled real mutational catalogs
# real catalogs and their respective ground truth solutions have to be in folder 'real mutational catalogs'
# ground truth solutions are obtained by analyzing complete real mutational catalogs with 50,000+ mutations
# GT_averaged_consensus2 is the ground truth that relies on the results obtained with SigProfilerAssignment (SPA),
# sigLASSO, and MuSiCal; only the signatures that are found by at least two tools are kept and their weights
# are averaged over those tools
def fit_with_cosmic3_subsampled_real_catalogs(code_name, real_catalogs = 'chosen_samples_WGS_PCAWG-input_profiles.dat', real_GT = 'chosen_samples_WGS_PCAWG-GT_averaged_consensus2.dat', GT_normalized = True):
    print('=== fitting subsampled real mutational catalogs using COSMICv3 ===')
    ttt = open('../running_times-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool), 'a')
    oname = '../individual_samples-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool)
    if not isfile(oname):                               # output file with results for individual samples
        aaa = open(oname, 'w')
        aaa.write('muts\trun\tsample\tTAE\twT\tn_FP\twT_FP\tn_FN\twT_FN\tP\tR\tF1\tMCC\n')
    else: aaa = open(oname, 'a')
    xxx = open('../stdout-{}.txt'.format(cfg.tool), 'a')
    yyy = open('../stderr-{}.txt'.format(cfg.tool), 'a')
    if not isdir('data'): mkdir('data')
    if not isdir('signature_results'): mkdir('signature_results')
    rng = np.random.default_rng(0)                      # initialize the RNG
    real_data_full = pd.read_csv('../real mutational catalogs/{}'.format(real_catalogs), sep = '\t', index_col = 'Type')
    GT = pd.read_csv('../real mutational catalogs/{}'.format(real_GT), sep = '\t', index_col = 'signature', comment = '#')
    if not GT_normalized:                               # if the ground truth is not normalized, normalize it now
        for col in GT.columns:
            GT[col] /= real_data_full.sum()[col]
    counter, new_col_name = 0, {}                       # simplify sample names; keep only the samples that are the input data and the GT
    for col in real_data_full.columns:
        if col in GT.columns:
            new_col_name[col] = 'S{}'.format(counter)
            counter += 1
    real_data_full = real_data_full[real_data_full.columns[real_data_full.columns.isin(new_col_name.keys())]]
    real_data_full = real_data_full.rename(columns = new_col_name)
    muts_full = real_data_full.sum()
    smallest_num_muts = muts_full.min()
    GT = GT[GT.columns[GT.columns.isin(new_col_name.keys())]]
    GT = GT.rename(columns = new_col_name)
    print('{} samples, smallest number of mutations in a sample is {}'.format(GT.shape[1], smallest_num_muts))
    which_num_muts = []                                 # for subsampling, keep all num_muts that are smaller than smallest_num_muts
    for num_muts in cfg.num_muts_list:
        if num_muts < smallest_num_muts: which_num_muts.append(num_muts)
    if smallest_num_muts / which_num_muts[-1] > 1.2:    # if smallest_num_muts is substantially higher than the last element, add it too
        which_num_muts.append(smallest_num_muts)
    for rep in range(cfg.num_realizations):
        print('\n{}: starting run {} for {}'.format(cfg.tool, rep, code_name))
        print('run\tsamples\tmuts\tMAE\twT\tn_FP\twT_FP\tn_FN\twT_FN\tP\tR\tF1\tMCC')
        for num_muts in which_num_muts:                 # gradually increase the number of mutations
            info_label = '{}\t{}\trun_{}\t{}'.format(cfg.tool, code_name, rep, num_muts)
            # prepare synthetic mutational catalogs
            muts = pd.Series(num_muts, index = real_data_full.columns)
            counts = prepare_data_from_real_data_by_subsampling(rng = rng, muts = muts, real_data = real_data_full)
            save_catalogs(counts = counts)
            # run the fitting tool defined in variable tool in MS_config.py
            timeout_run(info_label, ttt, xxx, yyy)
            # evaluate the estimated signature weights
            evaluate_real_catalogs(info_label, GT, muts, aaa)
        aaa.flush()
    # prepare a zip file with compressed (lzma) estimated signature weights for all cohorts
    system('zip ../signature_results-{}-{}-{}.zip signature_results/contribution-*.lzma'.format(cfg.WGS_or_WES, cfg.tool, code_name))
    shutil.rmtree('signature_results')  # remove all result files
    shutil.rmtree('data')               # remove the directory with synthetic data
    ttt.close()
    aaa.close()
    xxx.close()
    yyy.close()


# fitting synthetic mutational catalogs with empirical signatures weights, using all COSMICv3 signatures as a reference
# out-of-reference signatures are assigned the weights given in out_of_reference_weights (default: no out-of-reference signatures)
# examples: out_of_reference_weights = [0.2, 0.1] means that one out of reference signature has weight 20% and another has weight 10%
# out-of-reference signatures are the same for all samples in a cohort
def fit_external(input_catalog, catalog_GT, code_name, reference_signatures = None):
    ttt = open('../running_times-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool), 'a')
    xxx = open('../stdout-{}.txt'.format(cfg.tool), 'a')
    yyy = open('../stderr-{}.txt'.format(cfg.tool), 'a')
    if not isdir('data'): mkdir('data')
    if not isdir('signature_results'): mkdir('signature_results')
    counts = pd.read_csv(input_catalog, sep = None, index_col = 0, engine = 'python')   # load the mutational catalog
    counts = counts.reindex(index = cfg.input_sigs.index)                               # the desired order of SBS mutation types
    muts = counts.sum()                                                                 # number of mutations for each sample
    print('read mutational catalogs from {} samples'.format(muts.size))
    print('number of mutations varies from {} to {} (median is {:.0f})'.format(muts.min(), muts.max(), muts.median()))
    contribs = pd.read_csv(catalog_GT, sep = None, index_col = 0, engine = 'python')    # load the ground truth
    for sig in cfg.input_sigs.columns:
        if sig not in contribs.columns: contribs[sig] = 0
    contribs = contribs[cfg.input_sigs.columns]
    contribs = contribs.reindex(index = counts.columns)                                 # make the sample order same as the counts
    true_res = (contribs.T * muts).T                                                    # true absolute signature contributions
    save_catalogs(counts = counts)                                                      # save the catalogs for various fitting tools
    info_label = '{}\t{}'.format(cfg.tool, code_name)                                   # ID string for this run
    print('\n{}: starting the analysis of {} with the ground truth from {}'.format(cfg.tool, input_catalog, catalog_GT))
    if reference_signatures == None:
        print('all COSMICv3 signatures are used as a reference')
        which_setup = 'X'   # tool scripts starting with X use all COSMICv3 signatures as a reference
    else:
        print('a specified list of {} reference signatures is used as a reference'.format(len(reference_signatures)))
        which_setup = 'Y'   # tool scripts starting with Y use selected COSMICv3 signatures as a reference
        prepare_relevant_COSMIC(found_sigs = reference_signatures)
    print('external_id\t{}'.format(cfg.header_line))
    # run the fitting tool defined in variable tool in MS_config.py
    timeout_run(info_label, ttt, xxx, yyy, which_setup = which_setup)
    # evaluate the estimated signature weights
    evaluate_main(info_label, true_res.T, muts)
    oname = 'contribution-{}-{}-{}.dat.lzma'.format(cfg.WGS_or_WES, cfg.tool, code_name)
    rename('signature_results/{}'.format(oname), '../{}'.format(oname))
    shutil.rmtree('signature_results')  # remove all result files
    shutil.rmtree('data')               # remove the directory with synthetic data
    ttt.close()
    xxx.close()
    yyy.close()


# fitting real mutational catalogs using all COSMICv3 signatures as a reference
# three tools (SigProfilerAssignment, MuSiCal, and sigLASSO) that performed well on synthetic mutational catalogs are used
# differences (total absolute difference) between relative signature activities estimated by these tools are reported
def differences_real_samples():
    # run the analysis with three chosen tools
    for tool in ['SPA', 'MuSiCal', 'sigLASSO']:
        if tool in cfg.Python_tools:
            print('+++ started the analysis with {} (script PCAWG-{}.py)'.format(tool, tool))
            system('python3 PCAWG-{}.py'.format(tool))
        else:
            print('started the analysis with {} (script PCAWG-{}.R)'.format(tool, tool))
            system('Rscript PCAWG-{}.R'.format(tool))
        print('--- finished the analysis with {}\n\n'.format(tool))
    # load and pre-process the results
    catalogs = pd.read_csv('../real mutational catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat', sep = '\t', index_col = 'Type', comment = '#')
    muts = catalogs.sum()
    res1 = pd.read_csv('../WGS_PCAWG-146_samples-SPA-contribution.dat', index_col = 0)
    res1[res1 < 10] = 0
    res1 /= muts
    res2 = pd.read_csv('../WGS_PCAWG-146_samples-MuSiCal-contribution.dat', index_col = 0)
    res2[res2 < 10] = 0
    res2 /= muts
    res3 = pd.read_csv('../WGS_PCAWG-146_samples-sigLASSO-contribution.dat', index_col = 0)
    new_col_names = []
    for col in res3.columns:
        new_col_names.append(col.replace('..', '::').replace('.', '-'))
    res3.columns = new_col_names
    res3 *= muts
    res3[res3 < 10] = 0
    res3 /= muts
    # compute the total absolute difference between the three pairs of results
    rows = []
    for sample in muts.index:
        cancer_type = sample.split(':')[0]
        if cancer_type not in ['Skin-Melanoma', 'Stomach-AdenoCA', 'Lung-SCC']: cancer_type = 'Other'
        diff1 = (res1[sample] - res2[sample]).abs().sum()
        diff2 = (res1[sample] - res3[sample]).abs().sum()
        diff3 = (res2[sample] - res3[sample]).abs().sum()
        rows.append({'Sample': sample, 'Cancer type': cancer_type, 'Tot difference': (diff1 + diff2 + diff3) / 3})
    tmp = pd.DataFrame(rows)
    print('Number of samples of each cancer type:')
    print(tmp.groupby('Cancer type')['Tot difference'].count())
    print('\n\nTotal difference between the results stratified by cancer type:')
    print(tmp.groupby('Cancer type')['Tot difference'].median())
