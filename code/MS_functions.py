import MS_config as cfg
from MS_create_synthetic import *
from MS_run_and_evaluate import *
from scipy.spatial.distance import cosine as cosine_d
from os.path import isfile, isdir                       # OS-level utilities
from os import mkdir, rename, system, remove            # OS-level utilities
import shutil                                           # recursive directory removal
import sys                                              # emergency stop
from time import sleep                                  # to be able to pause when needed
import glob


# check if the provided dictoniary with weights of out of reference signatures is okay
def check_provided_weights(weights):
    if len(weights) > len(cfg.out_sigs):
        print('cannot introduce {} out-of-reference signatures (only {} signatures are in COSMICv3.3.1 but not in COSMICv3)'.format(len(weights), len(cfg.out_sigs)))
        sys.exit(1)
    if sum(weights) > 1 + cfg.EPSILON:
        print('sum of out-of-reference weights cannot exceed 1 (now {:.4f})'.format(sum(out_of_ref)))
        sys.exit(1)
    for weight in weights:
        if weight < 0:
            print('out-of-reference weights cannot be negative ({})'.format(weight))
            sys.exit(1)


# generate and save synthetic mutational catalogs with signature weights driven by COSMIC results on real tissue data
# out-of-reference signatures are assigned the weights given in out_of_reference_weights (default: no out-of-reference signatures)
def generate_synthetic_catalogs(cancer_types, out_of_reference_weights = []):
    tot_out_of_reference = sum(out_of_reference_weights)
    if tot_out_of_reference > 0:
        print('=== generating synthetic mutational catalogs based on real tissue data; out-of-reference signatures have joint weight {:.4f} ==='.format(tot_out_of_reference))
    else:
        print('=== generating synthetic mutational catalogs based on real tissue data ===')
    check_provided_weights(out_of_reference_weights)
    if not isdir('data'): mkdir('data')
    for cancer_type in cancer_types:
        rng = np.random.default_rng(0)                  # initialize the RNG
        empirical = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')                                # load empirical signature distribution for this cancer type
        for rep in range(cfg.num_realizations):         # run the chosen number of realizations
            print('\ngenerating catalog #{} for {}'.format(rep, cancer_type))
            who = rng.choice(empirical.shape[0], size = cfg.N_samples, replace = True)  # choose N_samples from all samples in the empirical signature distribution (with repetition)
            for num_muts in cfg.num_muts_list_short:    # gradually increase the number of mutations
                print('{} mutations...'.format(num_muts), end = ' ', flush = True)
                contribs = generate_weights_empirical(num_muts, empirical.iloc[who].astype('double').reset_index(drop = True))
                if tot_out_of_reference > 0:            # add out-of-reference signature for each sample (if needed)
                    for ix in contribs.index:
                        contribs.loc[ix] *= (1 - tot_out_of_reference)
                        new_active = rng.choice(cfg.out_sigs, size = len(out_of_reference_weights), replace = False)
                        for n, weight in enumerate(out_of_reference_weights):
                            contribs.loc[ix, new_active[n]] = weight
                info_label = '{}-w_{}-{}'.format(cancer_type, rep, num_muts)
                # generate and save synthetic mutational catalogs
                muts = pd.Series(num_muts, index = contribs.index)
                counts = prepare_data_from_signature_activity(rng = rng, muts = muts, contribs = contribs)
                save_catalogs(counts = counts, info_label = info_label)
                # keep only the active signatures
                tmp = contribs[contribs.columns[(contribs.sum(axis = 0) != 0)]]
                # save true signature contributions
                tmp.to_csv('data/true_weights-{}.dat'.format(info_label), sep = '\t', float_format = '%.6f')
    if not isdir('../generated_data'): mkdir('../generated_data')
    for f in glob.glob('data/*.dat'):
        rename(f, '../generated_data/{}'.format(f.split('/')[1]))
    print('\n\ndone, datasets saved in folder generated_data')


# generate and save synthetic mutational catalogs with signature weights driven by COSMIC results on real tissue data
# the samples differ in their numbers of mutations and weights of an out-of-reference signature (1 OOR signature per sample)
# ideally, num_samples should produce the same number of samples of each kind (kind == [muts, w_out])
# to achieve that, num_samples should be an integer multiple of len(muts) * len(weights_OOR)
def generate_mixed_synthetic_catalogs(num_samples, cancer_types, rng_seed = 0, code_name = 'MIXED', muts = [500, 1000, 2000, 4000, 8000], weights_OOR = [0, 0.1, 0.2, 0.3]):
    print('~~~ generating mixed synthetic mutational catalogs based on real tissue data ~~~')
    print('number of mutations ({}) and out-of-reference signature weights ({}) vary across the samples'.format(muts, weights_OOR))
    print('with {} samples in each catalog, there will be {:.0f} samples for each condition (number of mutations and out-of-reference weight)\n'.format(num_samples, num_samples / (len(muts) * len(weights_OOR))))
    if not isdir('data'):                                   # directory does not exist -> create it
        mkdir('data')
    else:                                                   # directory exists -> remove already existing files
        for f in glob.glob('data/*.dat'): remove(f)
    for cancer_type in cancer_types:
        print('generating catalog for {}'.format(cancer_type))
        rng = np.random.default_rng(rng_seed)               # initialize the RNG
        empirical = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')                         # load empirical signature distribution for this cancer type
        who = rng.choice(empirical.shape[0], size = num_samples, replace = True)    # choose from empirical signature activities
        contribs = generate_weights_empirical(500, empirical.iloc[who].astype('double').reset_index(drop = True), cohort_size = num_samples)
        num_muts = []
        for i, ix in enumerate(contribs.index):
            m = muts[(i // len(weights_OOR)) % len(muts)]   # number of mutations in this sample
            num_muts.append(m)
            w_OOR = weights_OOR[i % len(weights_OOR)]       # activities of out of reference signatures in this sample
            new_active = rng.choice(cfg.out_sigs)           # one out of reference signature is chosen for each sample
            if w_OOR > 0:
                contribs.loc[ix] *= (1 - w_OOR)
                contribs.loc[ix, new_active] = w_OOR
        info_label = '{}_{}'.format(code_name, cancer_type)
        # generate and save synthetic mutational catalogs
        sample_muts = pd.Series(num_muts, index = contribs.index)
        counts = prepare_data_from_signature_activity(rng = rng, muts = sample_muts, contribs = contribs)
        save_catalogs(counts = counts, info_label = info_label)
        # keep only the active signatures
        tmp = contribs[contribs.columns[(contribs.sum(axis = 0) != 0)]]
        tmp.index.name = 'Sample'
        # save true signature contributions
        tmp.to_csv('data/true_weights_{}.dat'.format(info_label), sep = '\t', float_format = '%.6f')
    if not isdir('../generated_data'): mkdir('../generated_data')
    for f in glob.glob('data/*.dat'):
        rename(f, '../generated_data/{}'.format(f.split('/')[1]))
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
# for examples, out_of_reference_weights = [0.2, 0.1] means that one out of reference signature has weight 20% and another has weight 10%
def fit_with_cosmic3_synthetic(cancer_types, code_name, out_of_reference_weights = [], save_true_weights = False):
    print('=== fitting synthetic mutational catalogs ===')
    tot_out_of_reference = sum(out_of_reference_weights)
    if tot_out_of_reference > 0:                                        # some samples have out of reference signatures
        print('generated samples have out of reference signatures ({})'.format(out_of_reference_weights))
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
                # prepare true signature weights for synthetic mutational catalogs
                contribs = generate_weights_empirical(num_muts, empirical.iloc[who].astype('double').reset_index(drop = True))
                if tot_out_of_reference > 0:                        # add out-of-reference signature for each sample (if needed)
                    for ix in contribs.index:
                        contribs.loc[ix] *= (1 - tot_out_of_reference)
                        new_active = rng.choice(cfg.out_sigs, size = len(out_of_reference_weights), replace = False)
                        for n, weight in enumerate(out_of_reference_weights):
                            contribs.loc[ix, new_active[n]] = weight
                if save_true_weights:                               # save true signature contributions
                    tmp = contribs[contribs.columns[(contribs.sum(axis = 0) != 0)]]
                    tmp.to_csv('data/true_weights_{}_{}_w{}_{}.dat'.format(cancer_type, code_name, rep, num_muts), sep = '\t', float_format = '%.6f')
                # prepare synthetic mutational catalogs
                info_label = '{}\t{}\tw_{}\t{}'.format(cfg.tool, code_name, rep, num_muts)
                muts = pd.Series(num_muts, index = contribs.index)
                counts = prepare_data_from_signature_activity(rng = rng, muts = muts, contribs = contribs)
                save_catalogs(counts = counts)
                true_res = (contribs.T * muts).T                    # data frame with true signature contributions
                timeout_run(info_label, ttt, xxx, yyy, extra_col = cancer_type)         # fit using the tool defined in MS_config.py
                evaluate_main(info_label, true_res.T, muts, extra_col = cancer_type)
        # zip signature activity estimates
        zip_oname = '../signature_results-{}-{}-{}-{}.zip'.format(cfg.WGS_or_WES, cfg.tool, code_name, cancer_type)
        system('zip {} signature_results/contribution-*.lzma signature_results/fit_quality_results-*.dat'.format(zip_oname))
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
    oname = '../results_individual_samples-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool)
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
    counter, new_col_name = 0, {}                       # keep only samples that are the input data and the GT, simplify sample names
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


# fit given mutational catalogs and compare the obtained signature estimates with the provided ground truth activities
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
        prepare_relevant_COSMIC(chosen_sigs = reference_signatures)
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


# fitting synthetic mutational catalogs that correspond to clones and their corresponding bulk (sum over clones)
# clone sizes are specified by the clone_sizes; signature weights are specified for each clone by sig_weights
# for fitting, all COSMICv3 signatures are used as a reference
# when tot_out_of_reference > 0, the specified fraction of mutations in each clone are due to an out-of-reference signature
# the out-of-reference signature is the same for all clones from a given patient (this assumption can be modified)
# when clone_error_mag > 0, some fraction of mutations from clone 1 are considered as mutations from clone 2 and vice versa
# (this feature is implemented only for samples with two clones)
def fit_synthetic_clones(code_name, clone_sizes, sig_weights, tot_out_of_reference = 0, clone_error_mag = 0, seed = 0):
    print('\n=== fitting synthetic mutational catalogs with simple clonal structures using COSMICv3 ===')
    print('clone sizes: {}'.format(', '.join([str(clone_sizes[c]) for c in clone_sizes])))
    print('out of reference weight: {:.2f}'.format(tot_out_of_reference))
    n_clones = len(clone_sizes.keys())
    n_bulk_samples = cfg.N_samples // (1 + n_clones)                # number of bulk samples whose catalogs are sums of the clones
    muts_bulk = sum(clone_sizes.values())
    if cfg.N_samples % (n_clones + 1) != 0:
        print('change N_samples in MS_config.py so that it can be divided by {} (number of clones + 1 for bulk); now it is {}'.format(n_clones + 1, cfg.N_samples))
        sys.exit(1)
    for c in clone_sizes:
        if c not in sig_weights:
            print('signature composition of clone {} is missing'.format(c))
            sys.exit(1)
        sig_weights_one_clone = sig_weights[c]
        sum_of_weights = sum(sig_weights_one_clone.values())
        if abs(sum_of_weights - 1) > cfg.EPSILON:
            print('signature weights in clone {} sum to {:4f}, not to 1'.format(c, sum_of_weights))
            sys.exit(1)
    ttt = open('../running_times-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool), 'a')
    xxx = open('../stdout-{}.txt'.format(cfg.tool), 'a')
    yyy = open('../stderr-{}.txt'.format(cfg.tool), 'a')
    if not isdir('data'): mkdir('data')
    if not isdir('signature_results'): mkdir('signature_results')
    rng = np.random.default_rng(seed)                               # initialize the RNG
    for rep in range(cfg.num_realizations):
        print('\n{}: starting run {} for {}'.format(cfg.tool, rep, code_name))
        print(cfg.header_line)
        contribs = pd.DataFrame(0, index = ['S{}'.format(x) for x in range(cfg.N_samples)], columns = cfg.input_sigs.columns, dtype = float)
        num_muts_list = []
        for n in range(cfg.N_samples):
            which_clone = n % (n_clones + 1)
            which_bulk = n // (n_clones + 1)
            if which_clone == 0:                                    # placeholder for future bulk samples
                num_muts_list.append(muts_bulk)
                contribs.loc['S{}'.format(n), 'SBS1'] = 1
            else:
                num_muts_list.append(clone_sizes[which_clone])
                sigs_here = sig_weights[which_clone]
                for sig in sigs_here:                           # scale the signature weights to make room for out of reference signatures
                    contribs.loc['S{}'.format(n), sig] = (1 - tot_out_of_reference) * sigs_here[sig]
                out_of_ref_sig = rng.choice(cfg.out_sigs)       # one out of reference signature is chosen for each sample
                contribs.loc['S{}'.format(n), out_of_ref_sig] = tot_out_of_reference
        # prepare synthetic mutational catalogs
        muts = pd.Series(num_muts_list, index = contribs.index)
        num_muts = muts.mean()                                      # mean number of mutations in the analyzed samples
        info_label = '{}\t{}~{:.0f}%\t{}\tw_{}\t{:.0f}'.format(cfg.tool, code_name, 100 * tot_out_of_reference, n_bulk_samples, rep, num_muts)
        counts = prepare_data_from_signature_activity(rng = rng, muts = muts, contribs = contribs)
        for n in range(n_bulk_samples):
            bulk_sample = n * (n_clones + 1)
            contribs.loc['S{}'.format(bulk_sample)] = 0             # signature contributions are computed as weighted average of clones
            counts['S{}'.format(bulk_sample)] = 0                   # drop generated counts and compute them by summing clones instead
            for c in range(1, n_clones + 1):
                contribs.loc['S{}'.format(bulk_sample)] += contribs.loc['S{}'.format(bulk_sample + c)] * clone_sizes[c] / muts_bulk
                counts['S{}'.format(bulk_sample)] += counts['S{}'.format(bulk_sample + c)]
        if clone_error_mag > 0:
            if n_clones != 2:
                print('clone_error_mag is so far implemented only for samples with two clones!')
                sys.exit(1)
            for n in range(n_bulk_samples):
                bulk_sample = n * (n_clones + 1)
                c1 = counts['S{}'.format(bulk_sample + 1)].copy()               # profile of clone C1
                c1_in_c2 = round(c1.sum() * clone_error_mag * rng.random())     # number of C1 mutations that will be counted as C2
                from_c1 = pd.Series(0, index = cfg.input_sigs.index)
                for _ in range(c1_in_c2):
                    ix = rng.choice(c1.index, p = c1 / c1.sum())
                    c1[ix] -= 1
                    from_c1[ix] += 1
                c2 = counts['S{}'.format(bulk_sample + 2)].copy()               # profile of clone C2
                c2_in_c1 = round(c2.sum() * clone_error_mag * rng.random())     # number of C2 mutations that will be counted as C1
                from_c2 = pd.Series(0, index = cfg.input_sigs.index)
                for _ in range(c2_in_c1):
                    ix = rng.choice(c2.index, p = c2 / c2.sum())
                    c2[ix] -= 1
                    from_c2[ix] += 1
                counts['S{}'.format(bulk_sample + 1)] += from_c2 - from_c1      # move the chosen mutations between the clones
                muts['S{}'.format(bulk_sample + 1)] += from_c2.sum() - from_c1.sum()
                counts['S{}'.format(bulk_sample + 2)] += from_c1 - from_c2
                muts['S{}'.format(bulk_sample + 2)] += from_c1.sum() - from_c2.sum()
                contribs_1 = (muts['S1'] - c1_in_c2) * contribs.loc['S{}'.format(bulk_sample + 1)] + c2_in_c1 * contribs.loc['S{}'.format(bulk_sample + 2)]
                contribs_1 /= contribs_1.sum()
                contribs_2 = (muts['S2'] - c2_in_c1) * contribs.loc['S{}'.format(bulk_sample + 2)] + c1_in_c2 * contribs.loc['S{}'.format(bulk_sample + 1)]
                contribs_2 /= contribs_2.sum()
                contribs.loc['S{}'.format(bulk_sample + 1)] = contribs_1
                contribs.loc['S{}'.format(bulk_sample + 2)] = contribs_2
        save_catalogs(counts = counts)
        # data frame with true signature contributions
        true_res = (contribs.T * muts).T
        # run the fitting tool defined in variable tool in MS_config.py
        timeout_run(info_label, ttt, xxx, yyy)
        # evaluate the consistency between the bulk and the clones
        evaluate_inconsistency(info_label, n_clones, muts)
        # evaluate the estimated signature weights
        oname = 'results-{}-{}'.format(cfg.tool, code_name)
        evaluate_main(info_label, true_res.T, muts, compress_result_file = False, oname = oname)
        rename('signature_results/{}-contribution.dat'.format(cfg.tool), 'signature_results/contribution-{}-{}-{}-w_{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool, code_name, rep, num_muts))
    # # prepare a zip file with compressed (lzma) estimated signature weights for all cohorts
    # system('zip ../signature_results-{}-{}-{}-.zip signature_results/contribution-*.lzma'.format(cfg.WGS_or_WES, cfg.tool, code_name))
    # shutil.rmtree('signature_results')  # remove all result files
    # shutil.rmtree('data')               # remove the directory with synthetic data
    ttt.close()
    xxx.close()
    yyy.close()


# fit a given mutational catalog, then go over all samples and generate num_bootstrap synthetic samples with signature
# activities given by the initially obtained estimates, fit the generated samples and compare their reconstruction
# metrics with the initially computed reconstruction metrics for the input sample
def fit_and_assess(input_catalog, code_name, num_bootstrap = 100, GT = None):
    rng = np.random.default_rng(0)                  # initialize the RNG
    cosmic3 = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38.txt', sep = '\t', index_col = 0)   # reference catalog used for fitting
    if not isdir('data'): mkdir('data')
    if not isdir('signature_results'): mkdir('signature_results')
    ttt = open('../running_times-{}-{}.dat'.format(cfg.WGS_or_WES, cfg.tool), 'a')
    xxx = open('../stdout-{}.txt'.format(cfg.tool), 'a')
    yyy = open('../stderr-{}.txt'.format(cfg.tool), 'a')
    oname = '../fit_assessment-{}-{}.dat'.format(cfg.tool, code_name)                   # fit quality evaluation results
    if isfile(oname): Qeval = open(oname, 'a')
    else:
        Qeval = open(oname, 'w')
        Qeval.write('sample\tL1\tL2\tcos\tL1_B\tL2_B\tcos_B\tz-L1\tz-L2\tz-cos\n')
    counts = pd.read_csv(input_catalog, sep = None, index_col = 0, engine = 'python')   # load the mutational catalog
    counts = counts.reindex(index = cfg.input_sigs.index)                               # the desired order of SBS mutation types
    muts = counts.sum()                                                                 # number of mutations for each sample
    save_catalogs(counts = counts)                                                      # save the catalogs for various fitting tools
    info_label = '{}\t{}'.format(cfg.tool, code_name)                                   # ID string for this run
    print('\nstarting the analysis of {} using {}'.format(input_catalog, cfg.tool))
    print('there are {} samples in the input data'.format(counts.shape[1]))
    which_setup = 'X'       # tool scripts starting with X use all COSMICv3 signatures as a reference
    # run the fitting tool defined in variable tool in MS_config.py
    timeout_run(info_label, ttt, xxx, yyy, which_setup = which_setup)
    print('the input data ({} samples) have been fitted with COSMICv3'.format(counts.shape[1]))
    if GT != None:
        oname = '../results_individual_samples-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool)
        if not isfile(oname):                                                           # fitting metrics for individual samples
            aaa = open(oname, 'w')
            aaa.write('muts\trun\tsample\tTAE\twT\tn_FP\twT_FP\tn_FN\twT_FN\tP\tR\tF1\tMCC\n')
        else: aaa = open(oname, 'a')
        GT = pd.read_csv(GT, sep = '\t', index_col = 0, comment = '#').T
        print('\nrun\tsamples\tmuts\tMAE\twT\tn_FP\twT_FP\tn_FN\twT_FN\tP\tR\tF1\tMCC')
        evaluate_real_catalogs(info_label, GT, muts, aaa, compress_result_file = False)
        aaa.close()
    res = pd.read_csv('signature_results/{}-contribution.dat'.format(cfg.tool), sep = ',', index_col = 0)
    rename('signature_results/{}-contribution.dat'.format(cfg.tool), '../{}_{}_contribution.dat'.format(cfg.tool, code_name))
    rows = []
    for n, col in enumerate(res.columns):           # assess the quality of fit for each sample
        normed_sample = counts[col] / muts[col]     # normalized sample
        res[col] /= res[col].sum()                  # force normalization of the estimated signature activities
        reconstruction = cosmic3.dot(res[col])      # reconstructed profile
        L1 = (normed_sample - reconstruction).abs().sum()
        L2 = np.sqrt(np.power(normed_sample - reconstruction, 2).sum())
        if reconstruction.sum() > 0: cos = 1 - cosine_d(normed_sample, reconstruction)
        else: cos = 0
        print('=== sample {} (muts = {}) ==='.format(col, muts[col]))
        print('original sample: L1 = {:.2f}, L2 = {:.2f}, cos = {:.2f}'.format(L1, L2, cos))
        Xmuts = pd.Series(muts[col], index = range(num_bootstrap))
        Xcontribs = res[col].copy()
        for sig in cfg.input_sigs.columns:
            if sig not in Xcontribs.index: Xcontribs[sig] = 0
        Xcontribs = pd.concat([Xcontribs] * num_bootstrap, axis = 1).T
        Xcounts = prepare_data_from_signature_activity(rng = rng, muts = Xmuts, contribs = Xcontribs)
        save_catalogs(counts = Xcounts)
        timeout_run(info_label, ttt, xxx, yyy, which_setup = which_setup)
        Xres = pd.read_csv('signature_results/{}-contribution.dat'.format(cfg.tool), sep = ',', index_col = 0)
        Xres /= Xres.sum()                          # bootsrap: force normalization of the estimated signature activities
        Xreconstruction = cosmic3.dot(Xres)         # bootsrap: reconstructed profiles
        Xnormed_sample = Xcounts / Xcounts.sum()    # bootsrap: normalized samples
        XL1 = (Xnormed_sample - Xreconstruction).abs().sum(axis = 0)
        XL2 = np.sqrt(np.power(Xnormed_sample - Xreconstruction, 2).sum(axis = 0))
        Xcos = []
        for sample in Xres.columns:
            Xcos.append(1 - cosine_d(Xnormed_sample[sample], Xreconstruction[sample]))
        Xcos = np.array(Xcos)
        print('      bootstrap: L1 = {:.2f}, L2 = {:.2f}, cos = {:.2f}'.format(XL1.mean(), XL2.mean(), Xcos.mean()))
        zL1 = (L1 - XL1.mean()) / XL1.std()
        zL2 = (L2 - XL2.mean()) / XL2.std()
        zcos = (Xcos.mean() - cos) / Xcos.std()
        print('      bootstrap: zL1 = {:.2f}, zL2 = {:.2f}, zcos = {:.2f}'.format(zL1, zL2, zcos))
        Qeval.write('{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(col, L1, L2, cos, XL1.mean(), XL2.mean(), Xcos.mean(), zL1, zL2, zcos))
        if n % 10 == 0: Qeval.flush()   # regularly flush the data to the output file
    shutil.rmtree('data')               # remove the directory with synthetic data
    Qeval.close()
    ttt.close()
    xxx.close()
    yyy.close()
