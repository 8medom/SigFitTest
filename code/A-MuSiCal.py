#!/usr/bin/env python3


import musical          # https://github.com/parklab/MuSiCal/tree/main
import pandas as pd


X = pd.read_csv('data/data_for_MutationalPatterns.dat', index_col = 0, sep = '\t')
W = pd.read_csv('../input/COSMIC_v3.4_SBS_GRCh37.txt', sep = '\t', index_col = 0)
MP_index = pd.read_csv('../input/mut_matrix_order_MutationalPatterns.dat', header = None).squeeze()
W = W.reindex(MP_index)
H, model = musical.refit.refit(X, W, method = 'likelihood_bidirectional', thresh = 0.001)   # Likelihood-based sparse NNLS
H.to_csv('signature_results/MuSiCal-contribution.dat')
