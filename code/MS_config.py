import pandas as pd                                     # for data structures
from os import mkdir                                    # OS-level utilities
from os.path import isfile, isdir                       # OS-level utilities
from sys import exit                                    # emergency stop


tool = 'SPA'                                            # which tool to use
WGS_or_WES = 'WGS'                                      # whether to use WGS or WES signatures (we eventually do everything for WGS only)
N_samples = 100                                         # how many synthetic samples to generate & analyze
num_realizations = 50                                    # how many combinations of weights do we generate
timeout_time = 3600                                     # calls to fitting methods are stopped after this time (increase to 8 hrs for mmsig)
num_muts_list = [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]         # which numbers of mutations to test for single-signature cohorts
num_muts_list_short = [100, 2000, 50000]                # which total numbers of mutations to test for heterogeneous cohorts
num_muts_list_long = [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400, 204800, 409600, 819200, 1638400]   # which numbers of mutations to test for single-signature cohorts with 5% of an out-of-reference signature
if WGS_or_WES == 'WGS': input_signatures = '../input/COSMIC_v3_SBS_GRCh38.txt'          # location of input WGS signatures
else: input_signatures = '../input/COSMIC_v3_SBS_GRCh38-WES.txt'                        # location of input WES signatures
artifact_sigs = ['SBS27', 'SBS43', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51', 'SBS52', 'SBS53', 'SBS54',
'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60']   # list of artifact signatures (they are not used in single-signature cohorts, for example)
tools_that_produce_relative_contributions = ['deconstructSigs', 'sigLASSO', 'sigfit', 'sigfit2', 'mmsig', 'SigsPack']    # these tools estimate relative signature weights (all others are assumed to estimate absolute signature weights)
Python_tools = ['SPSS', 'SPA']                          # these tools are Python-based (all others are assumed to be R-based)
tools_with_recommended_settings = ['sigfit', 'deconstructSigs']                         # these tools have some recommended settings that are evaluated as well


# create necessary directories (if needed)
if not isdir('data'): mkdir('data')
if not isdir('signature_results'): mkdir('signature_results')

input_sigs = pd.read_csv(input_signatures, sep = '\t', comment = '#', index_col = 'Type')
print('\n\n=== evaluation of fitting {} signatures ===\n\n'.format(WGS_or_WES))
if WGS_or_WES == 'WGS':     # use overall WGS trinucleotide frequencies as noise
    noise = pd.read_csv('../input/trinucleotide_frequencies-WGS.dat', sep = '\t', index_col = 'context')
else:                       # use overall WES trinucleotide frequencies as noise
    noise = pd.read_csv('../input/trinucleotide_frequencies-WES.dat', sep = '\t', index_col = 'context')


# load other auxiliary files
index_MutationalPatterns = pd.read_csv('../input/mut_matrix_order_MutationalPatterns.dat', header = None).squeeze()
index_sigLASSO = pd.read_csv('../input/sigLASSO-cosmic_v3_exo.txt', sep = ',').iloc[:, 0]
index_YAPSA = pd.read_csv('../input/mut_matrix_order_YAPSA.dat', header = None).squeeze()
order_spss = pd.read_csv('../input/mut_matrix_order_spss.dat', sep = ',')
noise = noise.reindex(input_sigs.index)
