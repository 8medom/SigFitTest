import numpy as np                                      # for numerics
import pandas as pd                                     # for data structures
from sklearn.metrics import precision_recall_curve, auc
from scipy.stats import wilcoxon, mannwhitneyu, ranksums, pearsonr
from lifelines import KaplanMeierFitter, CoxPHFitter
from os.path import isfile                              # OS-level utilities
from os import rename, system                           # OS-level utilities
import time                                             # to keep track of elapsed time
from subprocess import Popen, PIPE                      #
from threading import Timer                             #
import shlex                                            # to split a string using shell-like syntax
import lzma                                             # to compress some result files for future analysis
import MS_config as cfg                                 # all global stuff


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
        if extra_col == None:
            obj.write('sig\t#samples\t#muts\tMAE\tMAEstd\tRMSE\twT\t#eff\tMAE_TP\twT_FP\twT_FP_std\tn_FP\tPearson_r\tAUC_PR\n')
        else:
            obj.write('cancer_type\tweights\t#samples\t#muts\tMAE\tMAEstd\tRMSE\twT\t#eff\tMAE_TP\twT_FP\twT_FP_std\tn_FP\tPearson_r\tAUC_PR\n')
    else: obj = open(oname, 'a')
    return obj


# to make sure that the file with results starts with a header
def check_open_comparison(oname, extra_col):
    if not isfile(oname):
        obj = open(oname, 'w')
        if extra_col == None:
            obj.write('result_type\tscenario\tcohort_size\tweights\t#muts\tdifference_magnitude\tWilcoxon p-value\n')
        else:
            obj.write('cancer_type\tresult_type\tscenario\tcohort_size\tweights\t#muts\tdifference_magnitude\tWilcoxon p-value\n')
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
        ppp = open('problems-{}.txt'.format(cfg.tool), 'a')             # if there results cannot be loaded, save this information
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


# compute the evaluation metrics for df when true weights are stored in true_res
def eval_results(info_label, true_res, num_muts, df, extra_col, recommended = False):
    code_name, which = info_label.split('\t')[1], info_label.split('\t')[2]
    if recommended: output_file = check_open('results-{}-{}-{}_recommended.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), extra_col)
    else: output_file = check_open('results-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), extra_col)
    for ix in df.index:                 # to make sure that signatures from the results are not missing in the true result
        if ix not in true_res.index:
            true_res.loc[ix] = 0
    for ix in true_res.index:           # to make sure that signatures from the true result are not missing in the results
        if ix not in df.index:
            df.loc[ix] = 0
    err = np.abs(df.sub(true_res, axis = 0)).sum() / (2 * num_muts)     # total fitting error for each sample; factor 1/2 is introduced to make it bound to the range [0, 1]
    nMAE = err.mean()                           # mean of the fitting error over all samples
    nMAE_std = err.std()                        # std of the fitting error between the samples
    nRMSE = np.sqrt(np.power(err, 2).mean())    # root mean square witting error
    wtot_outside, wtot_outside_one, wtot_outside_squared, num_outside, MAE_inside, pearson_vals = 0, 0, 0, 0, 0, []
    if true_res.ndim == 1:              # when true_res is a vector, i.e., all samples have the same composition
        for sig in true_res.index:
            if true_res[sig] > 0:       # active signature
                if sig in df.index: vals = df.loc[sig]
                else: vals = 0
                MAE_inside += np.abs(vals - true_res[sig]).sum()        # the total error for true positives
            elif sig in df.index:                                       # is this signature in the results?
                wtot_outside += df.loc[sig].sum()                       # what is the total weight assigned to this false positive signature
                num_outside += (df.loc[sig] > 0).sum()                  # how many samples have this false positive signature
        for sample in df.columns:                                       # to estimate the standard deviation of wtot_FP
            either_pos = (true_res > 0) | (df[sample] > 0)
            if either_pos.sum() >= 3:
                pearson_vals.append(pearsonr(true_res[either_pos], df[sample][either_pos])[0])
            wtot_outside_one = 0
            for sig in df.index:
                if true_res[sig] == 0: wtot_outside_one += df.loc[sig, sample]
            wtot_outside_one /= num_muts
            wtot_outside_squared += wtot_outside_one * wtot_outside_one
            precision, recall, thresholds = precision_recall_curve(true_res > 0, df[sample])    # precision and recall computed by assuming that the active signatures are the true positives
            aucPR = auc(recall, precision)                              # area under the precision-recall curve
    elif true_res.ndim == 2:                                            # the same as above but assuming that the true weights are different for each sample (heterogeneous cohorts)
        for sample in true_res.columns:
            either_pos = (true_res[sample] > 0) | (df[sample] > 0)      # see how many signatures have positive true or estimated weight
            if either_pos.sum() >= 3:                                   # if they are at least three, compute the Pearson correlation between true and estimated weights (taking only those chosen signatures into account)
                pearson_vals.append(pearsonr(true_res[sample][either_pos], df[sample][either_pos])[0])
            wtot_outside_one = 0
            for sig in true_res.index:
                if true_res.loc[sig, sample] > 0:
                    if sig in df.index: val = df.loc[sig, sample]
                    else: val = 0
                    MAE_inside += np.abs(val - true_res.loc[sig, sample])   # the average error for true positives
                elif sig in df.index:                                   # is this signature in the results?
                    wtot_outside += df.loc[sig, sample]                 # the average weight assigned to this false positive signature
                    wtot_outside_one += df.loc[sig, sample]             # to focus on this one sample only
                    if df.loc[sig, sample] > 0: num_outside += 1        # count false positives
            wtot_outside_one /= num_muts
            wtot_outside_squared += wtot_outside_one * wtot_outside_one
            precision, recall, thresholds = precision_recall_curve(true_res[sample] > 0, df[sample])
            aucPR = auc(recall, precision)
    MAE_inside /= (num_muts * df.shape[1])                          # normalize all metrics
    wtot_outside /= (num_muts * df.shape[1])
    wtot_outside_squared /= df.shape[1]
    wtot_outside_squared -= wtot_outside * wtot_outside
    num_outside /= df.shape[1]
    weight_tot = df.sum().mean() / num_muts                         # sum of the assigned weights, averaged over samples, normalized
    w_norm = df / df.sum()
    n_eff = np.ma.masked_invalid(1 / np.power(w_norm, 2).sum()).mean()
    if len(pearson_vals) >= 3: pearson = np.nanmean(pearson_vals)
    else: pearson = np.nan
    if extra_col == None:
        out_string = '{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(which, df.shape[1], num_muts, nMAE, nMAE_std, nRMSE, weight_tot, n_eff, MAE_inside, wtot_outside, np.sqrt(wtot_outside_squared), num_outside, pearson, aucPR)
    else:
        out_string = '{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(extra_col, which, df.shape[1], num_muts, nMAE, nMAE_std, nRMSE, weight_tot, n_eff, MAE_inside, wtot_outside, np.sqrt(wtot_outside_squared), num_outside, pearson, aucPR)
    if recommended: print(out_string.strip() + ' (recommended settings)')
    else: print(out_string.strip())
    output_file.write(out_string)
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
            rename('signature_results/{}-contribution.dat'.format(cfg.tool), 'signature_results/contribution-{}-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-'), num_muts))
            system('lzma signature_results/contribution-{}-{}-{}.dat'.format(cfg.WGS_or_WES, info_label.replace('\t', '-'), num_muts))


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
    output_file = check_open_comparison('comparison_results-{}-{}-{}.dat'.format(cfg.WGS_or_WES, code_name, cfg.tool), extra_col)
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
