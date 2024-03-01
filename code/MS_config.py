import pandas as pd                                     # for data structures
from os import mkdir                                    # OS-level utilities
from os.path import isfile, isdir                       # OS-level utilities
from sys import exit                                    # emergency stop


tool = 'SPA'                                       # which tool to use
WGS_or_WES = 'WGS'                                      # whether to use WGS or WES signatures (we eventually do everything for WGS only)
N_samples = 100                                         # how many synthetic samples are in each cohort
num_realizations = 2                                   # how many independent cohorts do we generate
timeout_time = 1800                                     # calls to fitting methods are stopped after this time (increase to 8 hrs for mmsig)
num_muts_list = [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]         # numbers of mutations for simple cohorts
num_muts_list_short = [100, 2000, 50000]                                                # numbers of mutations for heterogeneous cohorts
if WGS_or_WES == 'WGS': input_signatures = '../input/COSMIC_v3.3.1_SBS_GRCh38.txt'      # location of input WGS signatures
else: input_signatures = '../input/COSMIC_v3_SBS_GRCh38-WES.txt'                        # location of input WES signatures
out_sigs = ['SBS10c', 'SBS10d', 'SBS86', 'SBS87', 'SBS88', 'SBS89', 'SBS90', 'SBS91', 'SBS92', 'SBS93', 'SBS94', 'SBS95']   # signatures that are in COSMICv3.3.1 but they are not in COSMICv3 (they can be used as out-of-reference signatures for tools that use COSMICv3 as a reference)
tools_that_produce_relative_contributions = ['deconstructSigs', 'sigLASSO', 'sigfit', 'sigfit2', 'mmsig', 'SigsPack']       # these tools estimate relative signature weights (all others are assumed to estimate absolute signature weights)
Python_tools = ['SPSS', 'SPA', 'MuSiCal', 'MuSiCalFull']                                # all other tools are assumed to be R-based
tools_with_recommended_settings = ['sigfit', 'deconstructSigs']                         # these tools have some recommended settings that are evaluated as well


# load the reference signatures
input_sigs = pd.read_csv(input_signatures, sep = '\t', comment = '#', index_col = 'Type')
print('\n\n~~~ starting SigFitTest for {} signatures ~~~\n\n'.format(WGS_or_WES))


# load other auxiliary files
index_MutationalPatterns = pd.read_csv('../input/mut_matrix_order_MutationalPatterns.dat', header = None).squeeze()
index_sigLASSO = pd.read_csv('../input/sigLASSO-cosmic_v3_exo.txt', sep = ',').iloc[:, 0]
index_YAPSA = pd.read_csv('../input/mut_matrix_order_YAPSA.dat', header = None).squeeze()
order_spss = pd.read_csv('../input/mut_matrix_order_spss.dat', sep = ',')
