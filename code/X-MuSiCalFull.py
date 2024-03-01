#!/usr/bin/env python3


import musical          # https://github.com/parklab/MuSiCal/tree/main
import pandas as pd
import numpy as np


X = pd.read_csv('data/data_for_MutationalPatterns.dat', index_col = 0, sep = '\t')
model = musical.DenovoSig(X,                      # parameter values copied from https://github.com/parklab/MuSiCal/blob/main/examples/example_full_pipeline.ipynb; three parameters were changed to allow for repeated runs in reasonable time (original parameter values are marked with #!)
                          min_n_components = 1,   # Minimum number of signatures to test
                          max_n_components = 20,  # Maximum number of signatures to test
                          init = 'random',        # Initialization method
                          method = 'mvnmf',       # mvnmf or nmf
                          n_replicates = 10,      # Number of mvnmf/nmf replicates to run per n_components  #! 20
                          ncpu = 1,               # Number of CPUs to use
                          max_iter = 10000,       # Maximum number of iterations for each mvnmf/nmf run   #! 100000
                          bootstrap = True,       # Whether or not to bootstrap X for each run
                          tol = 1e-6,             # Tolerance for claiming convergence of mvnmf/nmf   #! 1e-8
                          verbose = 1,            # Verbosity of output
                          normalize_X = False     # Whether or not to L1 normalize each sample in X before mvnmf/nmf
                         )
model.fit()


thresh_grid = np.array([0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2])   # four last elements of this array (0.5, 1., 2., 5.) have been removed to speed up this step; the optimal values reported in MuSiCalFull_denovo_stats.dat show that the upper bound of thresh_grid is sufficient
catalog = musical.load_catalog('COSMIC_v3_SBS_WGS')
W_catalog = catalog.W
model.assign_grid(W_catalog, 
                  method_assign = 'likelihood_bidirectional',   # Method for performing matching and refitting
                  thresh_match_grid = thresh_grid,              # Grid of threshold for matchinng
                  thresh_refit_grid = thresh_grid,              # Grid of threshold for refitting
                  thresh_new_sig = 0.0,                         # De novo signatures with reconstructed cosine similarity below this threshold will be considered novel
                  connected_sigs = False,                       # Whether or not to force connected signatures to co-occur
                  clean_W_s = False                             # An optional intermediate step to avoid overfitting to small backgrounds in de novo signatures for 96-channel SBS signatures
                 )


model.validate_grid(validate_n_replicates = 1,                  # Number of simulation replicates to perform for each grid point
                    grid_selection_method = 'pvalue',           # Method for selecting the best grid point
                    grid_selection_pvalue_thresh = 0.05         # Threshold used for selecting the best grid point
                   )
W_s = model.W_s
H_s = model.H_s
print(model.best_grid_point)
print(model.thresh_match)
print(model.thresh_refit)
H_s.to_csv('signature_results/MuSiCalFull-contribution.dat')    # save the fitting results for the optimal reference signatures
tmp = open('info_label.txt')
info_label = tmp.read()
tmp.close()
out = open('MuSiCalFull_denovo_stats.dat', 'a')                 # save information about the chosen reference signatures
out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(info_label, H_s.shape[1], model.n_components, model.thresh_match, model.thresh_refit, W_s.shape[1], ','.join(W_s.columns)))
out.close()
