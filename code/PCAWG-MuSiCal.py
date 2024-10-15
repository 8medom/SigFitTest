#!/usr/bin/env python3


import musical
import pandas as pd


# prepare input data
X = pd.read_csv('../real mutational catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat', sep = '\t', index_col = 'Type', comment = '#')
MP_index = pd.read_csv('../input/mut_matrix_order_MutationalPatterns.dat', header = None).squeeze()
X = X.reindex(MP_index)
W = pd.read_csv('../input/COSMIC_v3_SBS_GRCh38.txt', sep = '\t', index_col = 0)
W = W.reindex(MP_index)


# fit signatures
H, model = musical.refit.refit(X, W, method = 'likelihood_bidirectional', thresh = 0.001)   # Likelihood-based sparse NNLS
H.index.name = 'signature'
H.to_csv('../WGS_PCAWG-146_samples-MuSiCal-contribution.dat', float_format = '%.0f')
