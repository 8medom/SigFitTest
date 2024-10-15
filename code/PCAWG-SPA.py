import pandas as pd
from os import remove
import shutil                                           # recursive directory removal
from SigProfilerAssignment import Analyzer as Analyze   # https://github.com/AlexandrovLab/SigProfilerAssignment


# prepare input file
tmp = pd.read_csv('../real mutational catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat', sep = '\t', index_col = 'Type', comment = '#')
tmp.to_csv('146_samples.dat', sep = '\t')


# fit signatures
Analyze.cosmic_fit(samples = '146_samples.dat', input_type = 'matrix', context_type="96", output = 'signature_results', genome_build = 'GRCh38', cosmic_version = 3, make_plots = False, sample_reconstruction_plots = False)
df = pd.read_csv('signature_results/Assignment_Solution/Activities/Assignment_Solution_Activities.txt', sep = '\t', index_col = 'Samples').T
df.index.name = 'signature'
df.to_csv('../WGS_PCAWG-146_samples-SPA-contribution.dat')


# clean up
shutil.rmtree('signature_results')  # remove all result files
remove('146_samples.dat')
