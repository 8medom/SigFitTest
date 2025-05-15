import pandas as pd                     # for data structures


tool = 'sigLASSO'             # signature fitting tool to use; tool scripts named X-tool_name.R or X-tool_name.py
WGS_or_WES = 'WGS'                      # whether to use WGS or WES signatures
N_samples = 100                         # number of synthetic samples are in each cohort
num_realizations = 2                   # number of independent cohorts to generate
timeout_time = 1800                     # fitting tools processes are killed after this time (increase to 8 hrs for mmsig and signature.tools.lib)
num_muts_list = [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]     # numbers of mutations for simple cohorts
num_muts_list_short = [400, 2000, 10000]                                            # numbers of mutations for heterogeneous cohorts
EPSILON = 1e-8                                                                      # for floating-point comparisons
if WGS_or_WES == 'WGS': input_signatures = '../input/COSMIC_v3.3.1_SBS_GRCh38.txt'  # location of input WGS signatures
else: input_signatures = '../input/COSMIC_v3_SBS_GRCh38-WES.txt'                    # location of input WES signatures
out_sigs = ['SBS10c', 'SBS10d', 'SBS86', 'SBS87', 'SBS88', 'SBS89', 'SBS90', 'SBS91', 'SBS92', 'SBS93', 'SBS94', 'SBS95']   # signatures that are in COSMICv3.3.1 but they are not in COSMICv3 (they can be used as out-of-reference signatures for tools that use COSMICv3 as a reference)
tools_that_produce_relative_contributions = ['deconstructSigs', 'sigLASSO', 'sigfit', 'sigfit2', 'mmsig', 'SigsPack']       # these tools estimate relative signature weights (other tools are assumed to estimate absolute signature weights)
Python_tools = ['SPSS', 'SPA', 'MuSiCal', 'MuSiCalFull']                            # all other tools are assumed to be R-based
tools_with_recommended_settings = ['sigfit', 'deconstructSigs']                     # these tools have some recommended settings that are evaluated as well
header_line = 'weights\tsamples\tmuts\tTAE\tTAE_TP\tn_eff\twT\twT_FP\tn_FP\twT_FN\tn_FN\tP\tR\tS\tF1\tMCC\tPearson'              # header for the concise result tables (printed)
header_line_full = 'weights\tsamples\tmuts\tTAE\tTAE_std\tRMSE\twT\tn_eff\tTAE_TP\twT_FP\twT_FP_std\tn_FP\twT_FN\tn_FN\tP\tP_std\tR\tR_std\tS\tS_std\tF1\tF1_std\tMCC\tMCC_std\tPearson'      # header for the full result tables (saved)
top_sigs = {'Liver-HCC': ['SBS5', 'SBS12', 'SBS1', 'SBS29', 'SBS40', 'SBS4'], 'Stomach-AdenoCA': ['SBS5', 'SBS18', 'SBS1', 'SBS17b', 'SBS17a', 'SBS3'], 'Head-SCC': ['SBS5', 'SBS13', 'SBS2', 'SBS1', 'SBS40', 'SBS18'], 'ColoRect-AdenoCA': ['SBS5', 'SBS1', 'SBS18', 'SBS40', 'SBS44', 'SBS10a'], 'Lung-AdenoCA': ['SBS4', 'SBS5', 'SBS13', 'SBS2', 'SBS1', 'SBS40'], 'Skin-Melanoma': ['SBS7a', 'SBS7b', 'SBS5', 'SBS7c', 'SBS7d', 'SBS1'], 'Lymph-BNHL': ['SBS40', 'SBS5', 'SBS1', 'SBS9', 'SBS17b', 'SBS3'], 'CNS-GBM': ['SBS40', 'SBS1', 'SBS5', 'SBS11', 'SBS37', 'SBS15']}  # for each cancer type, five most active signatures are listed; average correlation between true and estimated signature weights in the cohort is computed for them and saved in the result files


# load the reference signatures
input_sigs = pd.read_csv(input_signatures, sep = '\t', comment = '#', index_col = 'Type')
print('\n\n~~~ starting SigFitTest for {} COSMIC signatures ~~~\n\n'.format(WGS_or_WES))


# load other auxiliary files
index_MutationalPatterns = pd.read_csv('../input/mut_matrix_order_MutationalPatterns.dat', header = None).squeeze()
index_sigLASSO = pd.read_csv('../input/sigLASSO-cosmic_v3_exo.txt', sep = ',', comment = '#').iloc[:, 0]
index_YAPSA = pd.read_csv('../input/mut_matrix_order_YAPSA.dat', header = None).squeeze()
order_spss = pd.read_csv('../input/mut_matrix_order_spss.dat', sep = ',')
