#!/usr/bin/env python3


from sigproSS import spss_pcwag     # https://github.com/AlexandrovLab/SigProfilerSingleSample
import pandas as pd


spss_pcwag.single_sample_pcwag('data/data_for_spss.dat', output = 'signature_results', sigbase = '../input/COSMIC_v3_SBS_GRCh38_for_spss.txt', n_cpu = 1)  # n_cpu = 1 to allow running in parallel with other runs
df = pd.read_csv('signature_results/sig_activities.txt', sep = '\t')
df = df.set_index('Cancer Types', drop = True)
df = df.drop('Similarity', axis = 1).T
df.to_csv('signature_results/SigProfiler-contribution.dat')
