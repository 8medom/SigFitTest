import numpy as np                                      # for numerics
import pandas as pd                                     # for data structures
from os.path import isfile                              # OS-level utilities
from os import rename, system                           # OS-level utilities
import time                                             # to keep track of elapsed time
from subprocess import Popen                            # to execute external scripts for fitting mutational signatures
from threading import Timer                             # to time-out a process after a pre-defined time (cfg.timeout_time)
import shlex                                            # to split a string using shell-like syntax
import lzma                                             # to compress result files for future analysis
from tqdm import tqdm                                   # for progress bar
from scipy.stats import ranksums, pearsonr, norm        # for the evaluation of results
from scipy.spatial.distance import cosine               # for sample reconstruction evaluation
import MS_config as cfg                                 # all global stuff


# shorten a string for text output on the screen
def shorten_string(s, max_length = 15):
    if len(s) <= max_length: return(s)
    else: return(s[:max_length - 3] + '...')


# prepare subsets of the default COSMIC signatures; these subsets can be either defined by which signatures are active for a given cancer type (when chosen_sigs == None) or by which signatures are found in previously computed fitting results (those signatures are then passed as a list in chosen_sigs)
# the generated subsets are formatted for all evaluated tools and their names have suffix RELEVANT
def prepare_relevant_COSMIC(cancer_type = None, also_artefacts = False, chosen_sigs = None):
    if cancer_type != None:                             # if cancer_type is given, read the list of signatures active in the tissue data
        tmp = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')
    df = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38.txt', sep = '\t', index_col = 'Type')
    cols_to_save = []
    for col in df.columns:                              # to get the default order of signatures
        if chosen_sigs != None:                         # if a list of signatures is provided, use only them
            if col in chosen_sigs: cols_to_save.append(col)
        else:                                           # otherwise, use the COSMIC tissue signatures
            if also_artefacts:
                if col in tmp.columns or col in cfg.artifact_sigs: cols_to_save.append(col)
            else:
                if col in tmp.columns: cols_to_save.append(col)
    df_to_save = df[cols_to_save]
    df_to_save.to_csv('../input/COSMIC_v3_SBS_GRCh38-RELEVANT.txt', sep = '\t')
    df = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38_for_MutationalPatterns.txt', sep = ',')
    df_to_save = df[cols_to_save]
    df_to_save.to_csv('../input/COSMIC_v3_SBS_GRCh38_for_MutationalPatterns-RELEVANT.txt', sep = ',', index = False)
    df = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38_for_YAPSA.txt', sep = ',', index_col = 0)
    df_to_save = df[cols_to_save]
    df_to_save.to_csv('../input/COSMIC_v3_SBS_GRCh38_for_YAPSA-RELEVANT.txt', sep = ',')
    df = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38_for_mmsig.txt', sep = ',')
    xcols = ['sub', 'tri']
    for col in cols_to_save: xcols.append(col)
    df_to_save = df[xcols]
    df_to_save.to_csv('../input/COSMIC_v3_SBS_GRCh38_for_mmsig-RELEVANT.txt', sep = ',', index = False)
    df = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38_for_sigLASSO.txt', sep = ',', index_col = 0)
    df_to_save = df[cols_to_save]
    df_to_save.to_csv('../input/COSMIC_v3_SBS_GRCh38_for_sigLASSO-RELEVANT.txt', sep = ',')
    df = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38_for_spss.txt', sep = ',')
    xcols = ['Type', 'SubType']
    for col in cols_to_save: xcols.append(col)
    df_to_save = df[xcols]
    df_to_save.to_csv('../input/COSMIC_v3_SBS_GRCh38_for_spss-RELEVANT.txt', sep = ',', index = False)
    df = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38_for_sigfit.txt', sep = ',', index_col = 0)
    df_to_save = df.loc[cols_to_save]
    df_to_save.to_csv('../input/COSMIC_v3_SBS_GRCh38_for_sigfit-RELEVANT.txt', sep = ',')


# to make sure that the file with results starts with a header
def check_open(oname, extra_col):
    if not isfile(oname):
        obj = open(oname, 'w')
        if extra_col == 'real_data':
            obj.write('run\tsamples\tmuts\tTAE\tTAE_std\twT\tn_FP\twT_FP\tn_FN\twT_FN\tP\tP_std\tR\tR_std\tF1\tF1_std\tMCC\tMCC_std\n')
        else:
            base = '\t'.join(cfg.header_line_full.split('\t')[1:])
            if extra_col == None: obj.write('sig_weights\t' + base + '\n')
            else: obj.write('cancer_type\tweights\t' + base + '\tr_S1\tr_S2\tr_S3\tr_S4\tr_S5\tr_S6\n')
    else: obj = open(oname, 'a')
    return obj


# to call commands with timeout (https://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout)
def timeout_run(info_label, ttt, xxx, yyy, timeout_sec = cfg.timeout_time, extra_col = None, which_setup = 'X'):
    if cfg.WGS_or_WES == 'WGS':
        if cfg.tool in cfg.Python_tools: cmd = 'python3 {}-{}.py'.format(which_setup, cfg.tool)
        else: cmd = 'Rscript {}-{}.R'.format(which_setup, cfg.tool)
    else:
        if cfg.tool in cfg.Python_tools: cmd = 'python3 {}-{}-WES.py'.format(which_setup, cfg.tool)
        else: cmd = 'Rscript {}-{}-WES.R'.format(which_setup, cfg.tool)
    t_start = time.time()
    proc = Popen(shlex.split(cmd), stdout = xxx, stderr = yyy)
    timer = Timer(timeout_sec, proc.kill)
    try:
        timer.start()
        stdout, stderr = proc.communicate()
    finally:
        timer.cancel()
    t_end = time.time()
    parts = info_label.split('\t')
    if extra_col == None:
        if info_label.count('\t') == 1:     # for the analysis of external data where the realization (parts[2]) and mean_num_muts (parts[3]) are undefined
            ttt.write('{}\t{:.1f}\n'.format(parts[1], t_end - t_start))
        else:                               # for simple cohorts
            ttt.write('{}\t{}\t{}\t{:.1f}\n'.format(parts[1], parts[2], parts[3], t_end - t_start))
    else:                                   # for heterogeneous cohorts where extra_col stores the cancer type
        ttt.write('{}_{}\t{}\t{}\t{:.1f}\n'.format(parts[1], extra_col, parts[2], parts[3], t_end - t_start))
    ttt.flush()


# load the result file, return absolute signature contributions
def load_results(info_label, muts, extra_col, fname):
    try:
        df = pd.read_csv(fname, sep = ',', index_col = 0)               # first column contains the fitted signatures
    except:
        if extra_col == None: out_string = 'no results for {}, {}\n'.format(extra_col, info_label)
        else: out_string = 'no results for {} & {}\n'.format(extra_col, info_label)
        print(out_string)
        ppp = open('../problems-{}.txt'.format(cfg.tool), 'a')             # if there results cannot be loaded, save this information
        ppp.write(out_string)
        ppp.close()
        return None
    else:
        if cfg.tool in cfg.tools_that_produce_relative_contributions:   # switch from relative to absolute signature contributions
            df = df * muts                                              # number of mutations is a series (different for each sample)
        df[df < 10] = 0                                                 # signatures contributing less than 10 mutations are set to zero
        return df


# compute the evaluation metrics for df when true weights are true_res
def eval_results(info_label, true_res, muts, df, extra_col, recommended = False, ref_sigs = 'COSMIC_v3'):
    if true_res.ndim == 1:                      # possible legacy issue: true_res is a vector (i.e., all samples have the same composition)
        print('true_res should be a DataFrame specifying the true signature weight for each sample, not a vector')
        sys.exit(1)
    if info_label.count('\t') == 1:             # when fitting external data
        code_name, line_start = info_label.split('\t')[1], info_label.split('\t')[1]
    else:                                       # for other cases, there is useful info (e.g., cancer type) in info_label.split('\t')[2]
        code_name, line_start = info_label.split('\t')[1], info_label.split('\t')[2]
    if recommended:
        output_file = check_open('../results-{}-{}-{}_recommended.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), extra_col)
    else:
        output_file = check_open('../results-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), extra_col)
    try:                                        # reference signatures used for fitting
        tmp = pd.read_csv('../input/{}_SBS_GRCh38.txt'.format(ref_sigs), sep = '\t', index_col = 'Type')
    except:
        print('cannot find file {} with the reference signatures used for fitting'.format('../input/{}_SBS_GRCh38.txt'.format(ref_sigs)))
        sys.exit(1)
    for sig in tmp.columns:
        if sig not in true_res.index:           # make sure that all reference signatures are included in the true result data frame
            true_res.loc[sig] = 0
        if sig not in df.index:                 # make sure that all reference signatures are included in the results data frame
            df.loc[sig] = 0
    for ix in df.index:                         # to make sure that signatures from the results are not missing in the true result
        if ix not in true_res.index:
            true_res.loc[ix] = 0
    for ix in true_res.index:                   # to make sure that signatures from the true result are not missing in the results
        if ix not in df.index:
            df.loc[ix] = 0
    err = df.sub(true_res, axis = 0).abs().sum() / (2 * muts)   # total absolute error (aka "fitting error") for each sample
    TAE = err.mean()                                            # mean of the fitting error over all samples
    TAE_std = err.std()                                         # std of the fitting error between the samples
    nRMSE = np.sqrt(np.power(err, 2).mean())                    # root mean square witting error
    wtot_FP, wtot_FP_squared, num_FP_sigs, wtot_FN, num_FN_sigs, MAE_active, pearson_vals, P_vals, R_vals, S_vals, F_vals, MCC_vals = 0, 0, 0, 0, 0, 0, [], [], [], [], [], []
    for sample in true_res.columns:
        either_pos = (true_res[sample] > 0) | (df[sample] > 0)      # see how many signatures have positive true or estimated weight
        if either_pos.sum() >= 3:                                   # if they are at least three, compute the Pearson correlation between true and estimated weights (taking only those chosen signatures into account)
            if np.std(true_res[sample][either_pos]) > cfg.EPSILON and np.std(df[sample][either_pos]) > cfg.EPSILON:   # to avoid one set of results to be all identical values (Pearson correlation then cannot be computed)
                pearson_vals.append(pearsonr(true_res[sample][either_pos], df[sample][either_pos])[0])
        wtot_FP_one_sample, num_TP, num_TN, num_FP, num_FN = 0, 0, 0, 0, 0
        for sig in true_res.index:
            if true_res.loc[sig, sample] > 0:                       # active signature
                MAE_active += np.abs(df.loc[sig, sample] - true_res.loc[sig, sample]) / muts[sample]    # error for active signatures
                if df.loc[sig, sample] == 0:                        # false negative
                    wtot_FN += true_res.loc[sig, sample] / muts[sample]
                    num_FN_sigs += 1
                    num_FN += 1
                else: num_TP += 1
            else:                                                   # inactive signature
                if df.loc[sig, sample] > 0:                         # false positive
                    wtot_FP += df.loc[sig, sample] / muts[sample]   # increment the average weight assigned to inactive signatures
                    wtot_FP_one_sample += df.loc[sig, sample]       # the same but for this one sample only
                    num_FP_sigs += 1
                    num_FP += 1
                else: num_TN += 1                                   # true negative
        wtot_FP_one_sample /= muts[sample]                          # total weight assigned to false positive signatures in the sample
        wtot_FP_squared += wtot_FP_one_sample * wtot_FP_one_sample
        if num_TP + num_FP > 0:                                     # precision & recall/sensitivity
            Psample = num_TP / (num_TP + num_FP)
            Rsample = num_TP / (num_TP + num_FN)
        else:
            Psample = 0
            Rsample = 0
        if num_TN + num_FP > 0:                                     # specificity (true negatives to all negatives)
            Ssample = num_TN / (num_TN + num_FP)
        else: Ssample = 0
        if Psample + Rsample > 0:                                   # F1 score
            Fsample = 2 * Psample * Rsample / (Psample + Rsample)
        else: Fsample = 0
        if Psample + Rsample > 0 and num_TN + num_FP > 0 and num_TN + num_FN > 0:   # Matthews correlation coefficient
            MCCsample = (num_TP * num_TN - num_FP * num_FN) / np.sqrt((num_TP + num_FP) * (num_TP + num_FN) * (num_TN + num_FP) * (num_TN + num_FN))
        else: MCCsample = 0
        P_vals.append(Psample)
        R_vals.append(Rsample)
        S_vals.append(Ssample)
        F_vals.append(Fsample)
        MCC_vals.append(MCCsample)
    MAE_active /= (2 * df.shape[1])                                 # normalize all metrics
    wtot_FP /= df.shape[1]
    wtot_FP_squared /= df.shape[1]
    wtot_FP_squared -= wtot_FP * wtot_FP
    num_FP_sigs /= df.shape[1]                                      # average number of false positive signatures per sample
    wtot_FN /= df.shape[1]
    num_FN_sigs /= df.shape[1]                                      # average number of false negative signatures per sample
    weight_tot = (df.sum() / muts).mean()                           # sum of the assigned weights, averaged over samples, normalized
    w_norm = df / df.sum()
    n_eff = np.ma.masked_invalid(1 / np.power(w_norm, 2).sum()).mean()  # effective number of estimated signatures per sample
    if np.ma.is_masked(n_eff): n_eff = np.nan                       # if n_eff is masked for all samples, mean cannot be computed
    if len(pearson_vals) >= 3: pearson = np.nanmean(pearson_vals)
    else: pearson = np.nan
    corr_out = ''
    if extra_col != None and extra_col in cfg.top_sigs.keys():      # for empirical weights, cohort-wide correlations for top signatures
        for sig in cfg.top_sigs[extra_col]:
            if sig not in true_res.index: true_res.loc[sig] = 0     # top signature missing in the true result (can happen rarely)
            if sig not in df.index: df.loc[sig] = 0                 # top signature missing in the results
            if true_res.loc[sig].std() > 0 and df.loc[sig].std() > 0:
                pearson_value = pearsonr(true_res.loc[sig], df.loc[sig])[0]
            else: pearson_value = np.nan
            corr_out = corr_out + '\t{:.4f}'.format(pearson_value)
    mean_num_muts = '{:.0f}'.format(muts.mean())
    full_string = '{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}{}\n'.format(line_start, df.shape[1], mean_num_muts, TAE, TAE_std, nRMSE, weight_tot, n_eff, MAE_active, wtot_FP, np.sqrt(wtot_FP_squared), num_FP_sigs, wtot_FN, num_FN_sigs, np.mean(P_vals), np.std(P_vals), np.mean(R_vals), np.std(R_vals), np.mean(S_vals), np.std(S_vals), np.mean(F_vals), np.std(F_vals), np.mean(MCC_vals), np.std(MCC_vals), pearson, corr_out)
    short_string = '{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}'.format(line_start, df.shape[1], mean_num_muts, TAE, MAE_active, n_eff, weight_tot, wtot_FP, num_FP_sigs, wtot_FN, num_FN_sigs, np.mean(P_vals), np.mean(R_vals), np.mean(S_vals), np.mean(F_vals), np.mean(MCC_vals), pearson)
    if extra_col == None:
        if recommended: print(short_string + '(recommended settings)')
        else: print(short_string)
        output_file.write(full_string)
    else:
        if recommended: print('{:15s}\t'.format(shorten_string(extra_col)) + short_string + '(recommended settings)')
        else: print('{:15s}\t'.format(shorten_string(extra_col)) + short_string)
        output_file.write(extra_col + '\t' + full_string)
    output_file.close()


# main function for evaluating the estimated signature weights
def evaluate_main(info_label, true_res, muts, extra_col = None, compress_result_file = True):
    df = load_results(info_label, muts, extra_col, 'signature_results/{}-contribution.dat'.format(cfg.tool))
    if df is not None:
        eval_results(info_label, true_res, muts, df, extra_col)
        if cfg.tool.startswith('sigfit'):    # use recommended options for sigfit
            dfx = pd.read_csv('signature_results/{}-contribution_lower90.dat'.format(cfg.tool), sep = ',') # lower estimate
            dfx = dfx.rename(columns={'Unnamed: 0': 'signature'})
            dfx = dfx.set_index('signature', drop = True)
            df2 = pd.read_csv('signature_results/{}-contribution.dat'.format(cfg.tool), sep = ',') # load the results
            df2 = df2.rename(columns={'Unnamed: 0': 'signature'})
            df2 = df2.set_index('signature', drop = True)
            df2[dfx < 0.01] = 0     # sigfit vignette: " In practice, ‘sufficiently non-zero’ means that the lower end of the Bayesian HPD interval (see the previous section) is above a threshold value close to zero (by default 0.01, and adjustable via the thresh argument)."
            df2 = (df2.T * muts).T  # number of mutations is a series (different for each sample)
            df2[df2 < 10] = 0
            eval_results(info_label, true_res, muts, df2, extra_col, recommended = True)
        elif cfg.tool == 'deconstructSigs': # use recommended options for deconstructSigs
            df2 = pd.read_csv('signature_results/{}-contribution.dat'.format(cfg.tool), sep = ',') # load the results
            df2 = df2.rename(columns={'Unnamed: 0': 'signature'})
            df2 = df2.set_index('signature', drop = True)
            df2[df2 < 0.06] = 0     # deconstructSigs vignette: "...by default, deconstructSigs uses "signature.cutoff = 0.06""
            df2 = (df2.T * muts).T
            df2[df2 < 10] = 0
            eval_results(info_label, true_res, muts, df2, extra_col, recommended = True)
        if compress_result_file:
            simplified_label = info_label.replace('\t', '-').replace('(', '_').replace(')', '')
            if extra_col != None: simplified_label = extra_col + '-' + simplified_label
            rename('signature_results/{}-contribution.dat'.format(cfg.tool), 'signature_results/contribution-{}-{}.dat'.format(cfg.WGS_or_WES, simplified_label))
            system('lzma signature_results/contribution-{}-{}.dat'.format(cfg.WGS_or_WES, simplified_label))


# evaluate the signature weights estimated for subsampled real mutational catalogs
def evaluate_real_catalogs(info_label, true_res, muts, aaa_file, compress_result_file = True, extra_col = 'real_data'):
    code_name, weights = info_label.split('\t')[1], info_label.split('\t')[2]
    df = load_results(info_label, muts, extra_col, 'signature_results/{}-contribution.dat'.format(cfg.tool))
    for sample in df.columns:                   # normalize the loaded results to weighted signature contributions
        df[sample] /= muts[sample]
    recommended = False
    if cfg.tool == 'deconstructSigs':           # by default, deconstructSigs uses "signature.cutoff = 0.06"
        df[df < 0.06] = 0
        recommended = True
    elif cfg.tool == 'sigfit':
        dfx = pd.read_csv('signature_results/{}-contribution_lower90.dat'.format(cfg.tool), sep = ',') # lower estimate
        dfx = dfx.rename(columns={'Unnamed: 0': 'signature'})
        dfx = dfx.set_index('signature', drop = True)
        df[dfx < 0.01] = 0                      # sigfit vignette: " In practice, ‘sufficiently non-zero’ means that the lower end of the Bayesian HPD interval (see the previous section) is above a threshold value close to zero (by default 0.01, and adjustable via the thresh argument)."
        recommended = True
    if recommended:
        output_file = check_open('../results-{}-{}-{}_recommended.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), extra_col)
    else: output_file = check_open('../results-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), extra_col)
    err = np.abs(df.sub(true_res, axis = 0)).sum() / 2
    TAE = err.mean()                            # mean of the fitting error over all samples
    TAE_std = err.std()                         # std of the fitting error between the samples
    num_FP_values, wtot_FP_values, num_FN_values, wtot_FN_values, P_vals, R_vals, F_vals, MCC_vals = [], [], [], [], [], [], [], []
    for sample in err.index:
        num_FP, num_FN, num_TP, num_TN, wtot_FP, wtot_FN = 0, 0, 0, 0, 0, 0     # initialize error counters for this sample
        for sig in true_res.index:
            if true_res.loc[sig, sample] > 0:   # active signature
                if df.loc[sig, sample] == 0:    # false negative
                    wtot_FN += true_res.loc[sig, sample]
                    num_FN += 1
                else: num_TP += 1
            else:                               # inactive signature
                if df.loc[sig, sample] > 0:     # false positive
                    wtot_FP += df.loc[sig, sample]
                    num_FP += 1
                else: num_TN += 1
        if num_TP + num_FP > 0:                 # precision & recall
            Psample = num_TP / (num_TP + num_FP)
            Rsample = num_TP / (num_TP + num_FN)
        else:
            Psample = 0
            Rsample = 0
        if Psample + Rsample > 0:               # F1 metric & Matthews correlation coefficient
            Fsample = 2 * Psample * Rsample / (Psample + Rsample)
            MCCsample = (num_TP * num_TN - num_FP * num_FN) / np.sqrt((num_TP + num_FP) * (num_TP + num_FN) * (num_TN + num_FP) * (num_TN + num_FN))
        else:
            Fsample = 0
            MCCsamples = 0
        P_vals.append(Psample)
        R_vals.append(Rsample)
        F_vals.append(Fsample)
        MCC_vals.append(MCCsample)
        aaa_file.write('{}\t{}\t{}\t{:.4f}\t{:.4f}\t{}\t{:.4f}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(muts[sample], weights, sample, err[sample], df[sample].sum(), num_FP, wtot_FP, num_FN, wtot_FN, Psample, Rsample, Fsample, MCCsample))
        num_FP_values.append(num_FP)
        wtot_FP_values.append(wtot_FP)
        num_FN_values.append(num_FN)
        wtot_FN_values.append(wtot_FN)
    mean_num_muts = '{:.0f}'.format(muts.mean())
    full_string = '{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(weights, df.shape[1], mean_num_muts, TAE, TAE_std, df.sum().mean(), np.mean(num_FP_values), np.mean(wtot_FP_values), np.mean(num_FN_values), np.mean(wtot_FN_values), np.mean(P_vals), np.std(P_vals), np.mean(R_vals), np.std(R_vals), np.mean(F_vals), np.std(F_vals), np.mean(MCC_vals), np.std(MCC_vals))
    short_string = '{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}'.format(weights, df.shape[1], mean_num_muts, TAE, df.sum().mean(), np.mean(num_FP_values), np.mean(wtot_FP_values), np.mean(num_FN_values), np.mean(wtot_FN_values), np.mean(P_vals), np.mean(R_vals), np.mean(F_vals), np.mean(MCC_vals))
    if recommended: print(short_string + ' (recommended settings)')
    else: print(short_string)
    output_file.write(full_string)
    output_file.close()
    if compress_result_file:
        rename('signature_results/{}-contribution.dat'.format(cfg.tool), 'signature_results/contribution-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-')))
        system('lzma signature_results/contribution-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-')))


# Wilcoxon rank-sum test for two groups (odd- and even-numbered samples)
def compare_weights(vals, info_label, extra_col, which_sig, difference_magnitude, output_file):
    x, y = [], []
    for n, ix in enumerate(vals.index):
        if n % 2 == 0: x.append(vals[ix])
        else: y.append(vals[ix])
    res = ranksums(x = x, y = y)
    out_string = '{}\t{}\t{}\t{}\t{:.4e}\n'.format(extra_col, which_sig, difference_magnitude, info_label, res.pvalue)
    print(out_string.strip())
    output_file.write(out_string)
    return res.pvalue


# compare the weights of signature which_sig between odd- and even-numbered samples
def compare_groups(info_label, which_sig, muts, difference_magnitude, true_res = None, extra_col = None):
    code_name, weights = info_label.split('\t')[1], info_label.split('\t')[2]
    oname = '../comparison_results-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool)
    if not isfile(oname):                           # open the output file
        output_file = open(oname, 'w')
        output_file.write('cancer_type\ttarget_signature\tdifference_magnitude\tresult_type\tscenario\tcohort_size\tweights\t#muts\tWilcoxon p-value\n')
    else: output_file = open(oname, 'a')
    if info_label.startswith('actual weights'):     # compare true weights, not the estimates
        pval = compare_weights(true_res[which_sig], info_label, extra_col, which_sig, difference_magnitude, output_file)
        output_file.close()
        return pval
    else:                                           # compare the estimated weights
        df = load_results(info_label, muts, extra_col, 'signature_results/{}-contribution.dat'.format(cfg.tool))
        if df is not None:
            compare_weights(df.loc[which_sig], info_label, extra_col, which_sig, difference_magnitude, output_file)
            mean_num_muts = '{:.0f}'.format(muts.mean())
            rename('signature_results/{}-contribution.dat'.format(cfg.tool), 'signature_results/contribution-{}-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-'), mean_num_muts))
            system('lzma signature_results/contribution-{}-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-'), mean_num_muts))
        output_file.close()


# assess the reconstruction of mutation profiles by various means (L1, L2, cosine, bootstrap,...)
# ref_sigs is the reference catalog that was used by the fitting tool to estimate signature activities
def evaluate_fits(info_label, input_profiles, extra_col = None, ref_sigs = 'COSMIC_v3'):
    ref_catalog = pd.read_csv('../input/{}_SBS_GRCh38.txt'.format(ref_sigs), sep = '\t', index_col = 0)
    code_name, weights = info_label.split('\t')[1], info_label.split('\t')[2]
    input_muts = input_profiles.sum()                                   # numbers of mutations in the samples
    try:
        fname = 'signature_results/{}-contribution.dat'.format(cfg.tool)
        estimated_sigs = pd.read_csv(fname, sep = ',', index_col = 0)
    except:
        if extra_col == None: out_string = 'fit quality evaluation: no results for {}, {}\n'.format(extra_col, info_label)
        else: out_string = 'fit quality evaluation: no results for {} & {}\n'.format(extra_col, info_label)
        print(out_string)
        return None
    if cfg.tool not in cfg.tools_that_produce_relative_contributions:   # switch from relative to absolute signature contributions
        estimated_sigs = estimated_sigs / input_muts
    worst_normalization = (estimated_sigs.sum() - 1).abs().max()
    if worst_normalization > cfg.EPSILON:
        print('there is a sample whose sum of estimated signature weights differs from 1 by {:.2e}'.format(worst_normalization))
        print('warning: current evaluation of signature fits is designed for estimated signature weights that sum to 1')
    oname = '../fit_quality_results-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool)
    if not isfile(oname):                                               # open the output file
        output_file = open(oname, 'w')
        if extra_col == None: output_file.write('cohort_size\tweights\tsample\t#muts\t\n')
        else: output_file.write('cancer_type\tcohort_size\tweights\tsample\t#muts\t\n')
    else: output_file = open(oname, 'a')
    # compute reconstruction metrics
    input_profiles_norm = input_profiles / input_muts                   # normalized input mutational catalogs
    for sample in tqdm(estimated_sigs.columns):                         # loop over samples
        profile_reconstructed = ref_catalog.dot(estimated_sigs[sample])
        cos_val = 1 - cosine(input_profiles_norm[sample], profile_reconstructed)
        L1_val = (input_profiles_norm[sample] - profile_reconstructed).abs().sum()
        L2_val = np.sqrt(np.power(input_profiles_norm[sample] - profile_reconstructed, 2).sum())
        if extra_col == None: output_file.write('{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(estimated_sigs.shape[1], weights, sample, input_muts[sample], cos_val, L1_val, L2_val))
        else: output_file.write('{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(extra_col, estimated_sigs.shape[1], weights, sample, input_muts[sample], cos_val, L1_val, L2_val))
    output_file.close()


def assess_fit_quality(info_label, input_catalog, sig_estimates, ref_sigs = 'COSMIC_v3', ref_genome = 'GRCh38', boot_realizations = 1000, extra_col = None, sample_info = None):
    # print('input mutational catalogs: {}'.format(input_catalog))
    # print('      signature estimates: {}'.format(sig_estimates))
    if 'COSMICv34' in sig_estimates: ref_sigs = 'COSMIC_v3.4'
    rng2 = np.random.default_rng(0)                                     # private RNG for this function only
    ref_catalog = pd.read_csv('../input/{}_SBS_{}.txt'.format(ref_sigs, ref_genome), sep = '\t', index_col = 0)
    if isinstance(input_catalog, str):                                  # path to the input catalog has been passed
        input_profiles = pd.read_csv(input_catalog, sep = '\t', index_col = 0)
    elif isinstance(input_catalog, pd.DataFrame):                       # DataFrame has been passed
        input_profiles = input_catalog
    else:
        print('do not know what to do with input_catalog {}'.format(type(input_catalog)))
        sys.exit(1)
    input_muts = input_profiles.sum()                                   # numbers of mutations in the samples
    try:
        estimated_sigs = pd.read_csv(sig_estimates, sep = ',', index_col = 0)
    except:
        print('cannot read signature estimates from {}'.format(sig_estimates))
        return None
    sample_sums = estimated_sigs.sum(axis = 0)
    if sample_sums.max() > 2:
        print('provided signature estimates have been inferred to be absolute')
        estimated_sigs = estimated_sigs / input_muts                    # switch from relative to absolute signature contributions
    else:
        print('provided signature estimates have been inferred to be relative')
    # file_name = sig_estimates.split('/')[-1]
    # if file_name.endswith('-contribution.dat'): file_name = file_name[:-17]
    # oname = 'signature_results/fit_quality_results-{}.dat'.format(file_name)
    # output_file = open(oname, 'w')
    simplified_label = info_label.replace('\t', '-').replace('(', '_').replace(')', '')
    if extra_col != None: simplified_label = extra_col + '-' + simplified_label
    output_file = open('signature_results/fit_quality_results-{}-{}.dat'.format(cfg.WGS_or_WES, simplified_label), 'w')
    worst_normalization = (estimated_sigs.sum() - 1).abs().max()
    if worst_normalization > cfg.EPSILON:
        warning = 'there is a sample whose sum of estimated signature weights differs from 1 by {:.2e}'.format(worst_normalization)
        print(warning)
        print('warning: current evaluation of signature fits is designed for estimated signature weights that sum to 1')
        output_file.write('# {}\n'.format(warning))
    if sample_info != None: extra = '\tsample_info'
    else: extra = ''
    output_file.write('sample\tmuts\tw_T\tcos\tL1\tL2\tE(L2_B)\tsig(L2_B)\tz-score\tp-value{}\n'.format(extra))
    input_profiles_norm = input_profiles / input_muts                   # normalized input mutational catalogs
    for i, sample in enumerate(tqdm(estimated_sigs.columns)):           # loop over samples
        profile_reconstructed = ref_catalog.dot(estimated_sigs[sample])
        profile_reconstructed_norm = profile_reconstructed / profile_reconstructed.sum()
        cos_val = 1 - cosine(input_profiles_norm[sample], profile_reconstructed)
        L1_val = (input_profiles_norm[sample] - profile_reconstructed).abs().sum()
        L2_val = np.sqrt(np.power(input_profiles_norm[sample] - profile_reconstructed, 2).sum())
        L2_bootstrap = []
        for boot_rep in range(boot_realizations):
            counts_bootstrap = np.zeros(96, dtype = int)                                                                # zero counts
            np.add.at(counts_bootstrap, rng2.choice(96, p = profile_reconstructed_norm, size = input_muts[sample]), 1)  # generate counts
            profile_bootstrap = counts_bootstrap.astype(float) / input_muts[sample]
            L2_boot = np.sqrt(np.power(profile_bootstrap - profile_reconstructed, 2).sum())
            L2_bootstrap.append(L2_boot)
        L2_boot_mean = np.mean(L2_bootstrap)
        L2_boot_std = np.std(L2_bootstrap)
        z_score = (L2_val - L2_boot_mean) / L2_boot_std
        p_value = 1 - norm.cdf(z_score)
        #! export to string differently above 0.001 and below (enginnering notation)
        # p_value = (np.array(L2_bootstrap) > L2_val).sum() / boot_realizations
        if sample_info != None: extra = '\t{}'.format(sample_info[i])
        else: extra = ''
        output_file.write('{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}{}\n'.format(sample, input_muts[sample], estimated_sigs[sample].sum(), cos_val, L1_val, L2_val, L2_boot_mean, L2_boot_std, z_score, p_value, extra))
    output_file.close()
