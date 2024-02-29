#!/usr/bin/env python3


import musical          # https://github.com/parklab/MuSiCal/tree/main
import pandas as pd


X = pd.read_csv('data/data_for_MutationalPatterns.dat', index_col = 0, sep = '\t')
catalog = musical.load_catalog('COSMIC_v3_SBS_WGS')
W = catalog.W
#H, model = musical.refit.refit(X, W, method='thresh_naive', thresh=0)                   # Naive NNLS
H, model = musical.refit.refit(X, W, method='likelihood_bidirectional', thresh=0.001)   # Likelihood-based sparse NNLS
H.to_csv('signature_results/MuSiCal-contribution.dat')
