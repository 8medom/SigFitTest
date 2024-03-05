import numpy as np                                      # for numerics
import pandas as pd                                     # for data structures
from os.path import isfile                              # OS-level utilities
from os import rename, system                           # OS-level utilities
import time                                             # to keep track of elapsed time
from subprocess import Popen                            # to execute external scripts for fitting mutational signatures
from threading import Timer                             # to time-out a process after a pre-defined time (cfg.timeout_time)
import shlex                                            # to split a string using shell-like syntax
import lzma                                             # to compress some result files for future analysis
from scipy.stats import ranksums, pearsonr              # for evaluation of results
import MS_config as cfg                                 # all global stuff


def shorten_string(s, max_length = 15):
    if len(s) <= max_length: return(s)
    else: return(s[:max_length - 3] + '...')


# prepare subsets of the default COSMIC signatures; these subsets can be either defined by which signatures are active for a given cancer type (when active_sigs == None) or by which signatures are active in previously computed fitting results (those signatures are then passed as a list in active_sigs)
# the generated subsets are formatted for all evaluated tools and their names have suffix RELEVANT
def prepare_relevant_COSMIC(cancer_type, also_artefacts = False, active_sigs = None):
    tmp = pd.read_csv('../cosmic tissue data/signature_contributions_{}_{}.dat'.format(cfg.WGS_or_WES, cancer_type), sep = '\t', index_col = 'Sample')
    df = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38.txt', sep = '\t', index_col = 'Type')
    cols_to_save = []
    for col in df.columns:                              # to get the default order of signatures
        if active_sigs != None:                         # if a list of signatures is provided, use only them
            if col in active_sigs: cols_to_save.append(col)
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
            obj.write('run\tsamples\tmuts\tMAE\tMAE_std\twT\tn_FP\twT_FP\tn_FN\twT_FN\tP\tP_std\tR\tR_std\tF1\tF1_std\n')
        else:
            base = '\t'.join(cfg.header_line_full.split('\t')[1:])
            if extra_col == None: obj.write('sig\t' + base + '\n')
            else: obj.write('cancer_type\tweights\t' + base + '\n')
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
        ttt.write('{}\t{}\t{}\t{:.1f}\n'.format(parts[1], parts[2], parts[3], t_end - t_start))
    else:
        ttt.write('{}_{}\t{}\t{}\t{:.1f}\n'.format(parts[1], extra_col, parts[2], parts[3], t_end - t_start))
    ttt.flush()


# try to load the result file
def load_results(info_label, num_muts, extra_col, fname):
    try:
        df = pd.read_csv(fname, sep = ',')
    except:
        if extra_col == None: out_string = 'no results for {}, {}\n'.format(extra_col, info_label)
        else: out_string = 'no results for {} & {}\n'.format(extra_col, info_label)
        print(out_string)
        ppp = open('../problems-{}.txt'.format(cfg.tool), 'a')             # if there results cannot be loaded, save this information
        ppp.write(out_string)
        ppp.close()
        return None
    else:
        df = df.rename(columns={'Unnamed: 0': 'signature'})
        df = df.set_index('signature', drop = True)
        if cfg.tool in cfg.tools_that_produce_relative_contributions:   # switch from relative to absolute signature contributions
            df *= num_muts
        df[df < 10] = 0                                                 # signatures contributing less than 10 mutations are set to zero
        return df


# compute the evaluation metrics for df when true weights are true_res
def eval_results(info_label, true_res, num_muts, df, extra_col, recommended = False):
    if true_res.ndim == 1:                      # possible legacy issue: true_res is a vector (i.e., all samples have the same composition)
        print('true_res should be a DataFrame specifying the true signature weight for each sample, not a vector')
        sys.exit(1)
    code_name, which = info_label.split('\t')[1], info_label.split('\t')[2]
    if recommended:
        output_file = check_open('../results-{}-{}-{}_recommended.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), extra_col)
    else:
        output_file = check_open('../results-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), extra_col)
    for ix in df.index:                         # to make sure that signatures from the results are not missing in the true result
        if ix not in true_res.index:
            true_res.loc[ix] = 0
    for ix in true_res.index:                   # to make sure that signatures from the true result are not missing in the results
        if ix not in df.index:
            df.loc[ix] = 0
    err = df.sub(true_res, axis = 0).abs().sum() / (2 * num_muts)   # total fitting error for each sample; factor 1/2 is introduced to make it bound to the range [0, 1]
    nMAE = err.mean()                           # mean of the fitting error over all samples
    nMAE_std = err.std()                        # std of the fitting error between the samples
    nRMSE = np.sqrt(np.power(err, 2).mean())    # root mean square witting error
    wtot_FP, wtot_FP_squared, num_FP_sigs, wtot_FN, num_FN_sigs, MAE_active, pearson_vals, P_vals, R_vals, F_vals = 0, 0, 0, 0, 0, 0, [], [], [], []
    for sample in true_res.columns:
        either_pos = (true_res[sample] > 0) | (df[sample] > 0)      # see how many signatures have positive true or estimated weight
        if either_pos.sum() >= 3:                                   # if they are at least three, compute the Pearson correlation between true and estimated weights (taking only those chosen signatures into account)
            if np.std(true_res[sample][either_pos]) > 1e-8 and np.std(df[sample][either_pos]) > 1e-8:   # to avoid one set of results to be all identical values (Pearson correlation then cannot be computed)
                pearson_vals.append(pearsonr(true_res[sample][either_pos], df[sample][either_pos])[0])
        wtot_FP_one_sample, num_TP, num_TN, num_FP, num_FN = 0, 0, 0, 0, 0
        for sig in true_res.index:
            if true_res.loc[sig, sample] > 0:                       # active signature
                MAE_active += np.abs(df.loc[sig, sample] - true_res.loc[sig, sample])   # increment the error for active signatures
                if df.loc[sig, sample] == 0:                        # false negative
                    wtot_FN += true_res.loc[sig, sample]
                    num_FN_sigs += 1
                    num_FN += 1
                else: num_TP += 1
            else:                                                   # inactive signature
                if df.loc[sig, sample] > 0:                         # false positive
                    wtot_FP += df.loc[sig, sample]                  # increment the average weight assigned to inactive signatures
                    wtot_FP_one_sample += df.loc[sig, sample]       # the same but for this one sample only
                    num_FP_sigs += 1
                    num_FP += 1
                else: num_TN += 1                                   # true negative
        wtot_FP_one_sample /= num_muts
        wtot_FP_squared += wtot_FP_one_sample * wtot_FP_one_sample
        Psample = num_TP / (num_TP + num_FP)
        Rsample = num_TP / (num_TP + num_FN)
        if Psample + Rsample > 0: Fsample = 2 * Psample * Rsample / (Psample + Rsample)
        else: Fsample = 0
        P_vals.append(Psample)
        R_vals.append(Rsample)
        F_vals.append(Fsample)
    MAE_active /= (2 * num_muts * df.shape[1])                      # normalize all metrics
    wtot_FP /= (num_muts * df.shape[1])
    wtot_FP_squared /= df.shape[1]
    wtot_FP_squared -= wtot_FP * wtot_FP
    num_FP_sigs /= df.shape[1]                                      # average number of false positive signatures per sample
    wtot_FN /= (num_muts * df.shape[1])
    num_FN_sigs /= df.shape[1]                                      # average number of false negative signatures per sample
    weight_tot = df.sum().mean() / num_muts                         # sum of the assigned weights, averaged over samples, normalized
    w_norm = df / df.sum()
    n_eff = np.ma.masked_invalid(1 / np.power(w_norm, 2).sum()).mean()  # effective number of estimated signatures per sample
    if len(pearson_vals) >= 3: pearson = np.nanmean(pearson_vals)
    else: pearson = np.nan
    full_string = '{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(which, df.shape[1], num_muts, nMAE, nMAE_std, nRMSE, weight_tot, n_eff, MAE_active, wtot_FP, np.sqrt(wtot_FP_squared), num_FP_sigs, wtot_FN, num_FN_sigs, np.mean(P_vals), np.std(P_vals), np.mean(R_vals), np.std(R_vals), np.mean(F_vals), np.std(F_vals), pearson)
    short_string = '{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}'.format(which, df.shape[1], num_muts, nMAE, nRMSE, weight_tot, n_eff, MAE_active, wtot_FP, num_FP_sigs, wtot_FN, num_FN_sigs, np.mean(P_vals), np.mean(R_vals), np.mean(F_vals), pearson)
    if extra_col == None:
        if recommended: print(short_string + '(recommended settings)')
        else: print(short_string)
        output_file.write(full_string)
    else:
        if recommended: print(shorten_string(extra_col) + '\t' + short_string + '(recommended settings)')
        else: print(shorten_string(extra_col) + '\t' + short_string)
        output_file.write(extra_col + '\t' + full_string)
    output_file.close()


# main function for evaluating the estimated signature weights
def evaluate_main(info_label, true_res, num_muts, extra_col = None, compress_result_file = True):
    df = load_results(info_label, num_muts, extra_col, 'signature_results/{}-contribution.dat'.format(cfg.tool))
    if df is not None:
        eval_results(info_label, true_res, num_muts, df, extra_col)
        if cfg.tool.startswith('sigfit'):    # use recommended options for sigfit
            dfx = pd.read_csv('signature_results/{}-contribution_lower90.dat'.format(cfg.tool), sep = ',') # lower estimate
            dfx = dfx.rename(columns={'Unnamed: 0': 'signature'})
            dfx = dfx.set_index('signature', drop = True)
            df2 = pd.read_csv('signature_results/{}-contribution.dat'.format(cfg.tool), sep = ',') # load the results
            df2 = df2.rename(columns={'Unnamed: 0': 'signature'})
            df2 = df2.set_index('signature', drop = True)
            df2[dfx < 0.01] = 0    # sigfit vignette: " In practice, ‘sufficiently non-zero’ means that the lower end of the Bayesian HPD interval (see the previous section) is above a threshold value close to zero (by default 0.01, and adjustable via the thresh argument)."
            df2 *= num_muts
            df2[df2 < 10] = 0
            eval_results(info_label, true_res, num_muts, df2, extra_col, recommended = True)
        elif cfg.tool == 'deconstructSigs': # use recommended options for deconstructSigs
            df2 = pd.read_csv('signature_results/{}-contribution.dat'.format(cfg.tool), sep = ',') # load the results
            df2 = df2.rename(columns={'Unnamed: 0': 'signature'})
            df2 = df2.set_index('signature', drop = True)
            df2[df2 < 0.06] = 0    # by default, deconstructSigs uses "signature.cutoff = 0.06"
            df2 *= num_muts
            df2[df2 < 10] = 0
            eval_results(info_label, true_res, num_muts, df2, extra_col, recommended = True)
        if compress_result_file:
            simplified_label = info_label.replace('\t', '-').replace('(', '_').replace(')', '')
            rename('signature_results/{}-contribution.dat'.format(cfg.tool), 'signature_results/contribution-{}-{}.dat'.format(cfg.WGS_or_WES, simplified_label))
            system('lzma signature_results/contribution-{}-{}.dat'.format(cfg.WGS_or_WES, simplified_label))


# evaluate the signature weights estimated for subsampled real mutational catalogs
def evaluate_real_catalogs(info_label, true_res, num_muts, aaa_file, compress_result_file = True, extra_col = 'real_data'):
    code_name, which = info_label.split('\t')[1], info_label.split('\t')[2]
    df = load_results(info_label, num_muts, extra_col, 'signature_results/{}-contribution.dat'.format(cfg.tool))
    for col in df.columns:                      # normalize the loaded results to weighted signature contributions
        df[col] /= num_muts
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
    nMAE = err.mean()                           # mean of the fitting error over all samples
    nMAE_std = err.std()                        # std of the fitting error between the samples
    num_FP_values, wtot_FP_values, num_FN_values, wtot_FN_values, P_vals, R_vals, F_vals = [], [], [], [], [], [], []
    for sample in err.index:
        num_FP, num_FN, num_TP, num_TN, wtot_FP, wtot_FN = 0, 0, 0, 0, 0, 0
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
        Psample = num_TP / (num_TP + num_FP)
        Rsample = num_TP / (num_TP + num_FN)
        if Psample + Rsample > 0: Fsample = 2 * Psample * Rsample / (Psample + Rsample)
        else: Fsample = 0
        P_vals.append(Psample)
        R_vals.append(Rsample)
        F_vals.append(Fsample)
        aaa_file.write('{}\t{}\t{}\t{:.4f}\t{:.4f}\t{}\t{:.4f}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(num_muts, which, sample, err[sample], df[sample].sum(), num_FP, wtot_FP, num_FN, wtot_FN, Psample, Rsample, Fsample))
        num_FP_values.append(num_FP)
        wtot_FP_values.append(wtot_FP)
        num_FN_values.append(num_FN)
        wtot_FN_values.append(wtot_FN)
    full_string = '{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(which, df.shape[1], num_muts, nMAE, nMAE_std, df.sum().mean(), np.mean(num_FP_values), np.mean(wtot_FP_values), np.mean(num_FN_values), np.mean(wtot_FN_values), np.mean(P_vals), np.std(P_vals), np.mean(R_vals), np.std(R_vals), np.mean(F_vals), np.std(F_vals))
    short_string = '{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}'.format(which, df.shape[1], num_muts, nMAE, df.sum().mean(), np.mean(num_FP_values), np.mean(wtot_FP_values), np.mean(num_FN_values), np.mean(wtot_FN_values), np.mean(P_vals), np.mean(R_vals), np.mean(F_vals))
    if recommended: print(short_string + ' (recommended settings)')
    else: print(short_string)
    output_file.write(full_string)
    output_file.close()
    if compress_result_file:
        rename('signature_results/{}-contribution.dat'.format(cfg.tool), 'signature_results/contribution-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-')))
        system('lzma signature_results/contribution-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-')))


# Wilcoxon rank-sum test for two groups (odd- and even-numbered samples)
def compare_weights(vals, info_label, extra_col, difference_magnitude, output_file):
    x, y = [], []
    for n, ix in enumerate(vals.index):
        if n % 2 == 0: x.append(vals[ix])
        else: y.append(vals[ix])
    res = ranksums(x = x, y = y)
    out_string = '{}\t{}\t{}\t{:.4e}\n'.format(extra_col, info_label, difference_magnitude, res.pvalue)
    print(out_string.strip())
    output_file.write(out_string)
    return res.pvalue


# compare the weights of signature which_sig between odd- and even-numbered samples
def compare_groups(info_label, which_sig, num_muts, difference_magnitude, true_res = None, extra_col = None):
    code_name, which = info_label.split('\t')[1], info_label.split('\t')[2]
    oname = 'comparison_results-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool)
    if not isfile(oname):                           # open the output file
        output_file = open(oname, 'w')
        output_file.write('cancer_type\tresult_type\tscenario\tcohort_size\tweights\t#muts\tdifference_magnitude\tWilcoxon p-value\n')
    else: output_file = open(oname, 'a')
    if info_label.startswith('actual weights'):     # compare true weights, not the estimates
        pval = compare_weights(true_res[which_sig], info_label, extra_col, difference_magnitude, output_file)
        output_file.close()
        return pval
    else:                                           # compare the estimated weights
        df = load_results(info_label, num_muts, extra_col, 'signature_results/{}-contribution.dat'.format(cfg.tool))
        if df is not None:
            compare_weights(df.loc[which_sig], info_label, extra_col, difference_magnitude, output_file)
            rename('signature_results/{}-contribution.dat'.format(cfg.tool), 'signature_results/contribution-{}-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-'), num_muts))
            system('lzma signature_results/contribution-{}-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-'), num_muts))
        output_file.close()
