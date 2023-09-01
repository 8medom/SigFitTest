import numpy as np                                      # for numerics
import pandas as pd                                     # for data structures
import sys                                              # emergency stop
from os.path import isfile                              # OS-level utilities
import MS_config as cfg                                 # all global stuff


# generate the number of mutations in all contexts for one sample
def generate_one_sample(rng, sample_name, num_muts, contribs_given):
    if contribs_given.sum() > 1 + 1e-14:    # check if the sum of signature contributions is not too large
        print('the sum of signature contributions cannot be greater than 1, {:.4f} was given in {}'.format(contribs_given.sum(), sample_name))
        sys.exit(1)
    if np.any(contribs_given < 0):          # there must be no negative signature contributions
        print('signature contributions must be non-negative, now the smallest one is {:.1e} in {}'.format(contribs_given.min(), sample_name))
        sys.exit(1)
    context_weights = pd.Series(0, index = cfg.input_sigs.index)                        # compute weights of tri-nucleotide contexts
    for sig in cfg.input_sigs.columns:                                                  # for each context, the weights is a weighted sum over all signatures
        context_weights += contribs_given[sig] * cfg.input_sigs[sig]
    noise_mag = 1 - contribs_given.sum()
    if noise_mag > 0:                      # when the given signature contributions are less than one in total, the rest is expected to be noise; we use the overall WGS context frequencies as the noise
        for ix in context_weights.index:
            context_weights.loc[ix] += noise_mag * cfg.noise.loc[ix].values[0]
    context_weights /= sum(context_weights)                                             # make sure that the weights of all contexts sum to 1
    counts = np.zeros(96, dtype = int)                                                  # empty count array
    np.add.at(counts, rng.choice(96, p = context_weights, size = num_muts), 1)          # generate mutation positions
    return pd.DataFrame(counts, index = context_weights.index, columns = [sample_name]) # return the count array as a data frame


# generate the mutational catalogs and save them in formats that can be used by the evaluated tools
def prepare_data(rng, num_muts, contribs, cohort_size = cfg.N_samples):
    if contribs.ndim == 1:      # all samples have the same signature contributions
        counts = pd.concat([generate_one_sample(rng, 'S{}'.format(i), num_muts, contribs) for i in range(cohort_size)], axis = 1)
    elif contribs.ndim == 2:    # samples have different (empirically-driven?) signature contributions
        if contribs.shape[0] != cohort_size:
            print('mismatch between the provided signature contribution array with shape {} and the desired number of samples ({})'.format(contribs.shape, cohort_size))
            sys.exit(1)
        counts = pd.concat([generate_one_sample(rng, 'S{}'.format(i), num_muts, contribs.iloc[i]) for i in range(cohort_size)], axis = 1)
    counts.to_csv('data/data_for_deconstructSigs.dat', sep = '\t')
    counts = counts.reindex(index = cfg.index_MutationalPatterns)
    counts.index.name = 'Type'
    counts.to_csv('data/data_for_MutationalPatterns.dat', sep = '\t')
    counts3 = counts.copy()
    new_index = []
    for ix in counts3.index:
        new_index.append(ix[2:5] + ' ' + ix[0] + ix[2] + ix[6])
    counts3.index = new_index
    counts3 = counts3.reindex(index = cfg.index_YAPSA)
    counts3.to_csv('data/data_for_YAPSA.dat', sep = '\t')
    counts2 = counts.copy()
    mut_type, trinuc = [], []
    for ix in counts2.index:
        mut_type.append(ix[2:5])
        trinuc.append(ix[0] + ix[2] + ix[6])
    counts2['Mutation type']= mut_type
    counts2['Trinucleotide'] = trinuc
    cols = list(counts2.columns)
    new_cols = [cols[-2]] + [cols[-1]] + cols[:-2]
    counts2 = counts2[new_cols]
    counts2.to_csv('data/data_for_spss.dat', sep = ',', index = False)
    if cfg.tool == 'sigfit':    # data files for sigfit are large -> save them only when really needed
        ox = open('data/data_for_sigfit.dat', 'w')
        ox.write('Sample\tRef\tAlt\tTrinuc\n')
        for col in counts.columns:
            for ix in counts.index:
                counter = counts.loc[ix, col]
                s1, s2, s3 = ix[2], ix[4], ix[0] + ix[2] + ix[6]
                while counter > 0:
                    ox.write('{}\t{}\t{}\t{}\n'.format(col, s1, s2, s3))
                    counter -= 1
        ox.close()
    new_index = []
    for ix in counts.index:
        new_index.append(ix.replace('[', '').replace(']', ''))
    counts.index = new_index
    counts = counts.reindex(index = cfg.index_sigLASSO)
    counts.index.name = 'Type'
    counts.to_csv('data/data_for_sigLASSO_spectrum.dat', sep = '\t')


# generate empirically-driven signature contributions
def generate_weights_empirical(num_muts, empirical_sub, cohort_size = cfg.N_samples, noise = 0):
    for ix in empirical_sub.index:      # re-weight the signature contributions to match the required number of mutations
        empirical_sub.loc[ix] = empirical_sub.loc[ix] * (1 - noise) * num_muts / empirical_sub.loc[ix].sum()
    empirical_sub[empirical_sub < 10] = 0   # we remove weak signatures (those with the expected number of mutations below 10) because they cannot be recovered by the methods anyway (as we use the same threshold of 10 mutations for them)
    if (empirical_sub.sum(axis = 1) == 0).sum() > 0:
        print('for num_muts = {} and realization {}, there is a sample with all absolute signature weights below 10\nwe skip this iteration'.format(num_muts, rep))
        return None
    for ix in empirical_sub.index:      # re-weight the signature contributions again to compensate for possible zeros introduced above
        empirical_sub.loc[ix] *= (1 - noise) / empirical_sub.loc[ix].sum()
    contribs = pd.DataFrame(0, index = ['S{}'.format(x) for x in range(cohort_size)], columns = cfg.input_sigs.columns, dtype = float)  # move from the fitted signatures to a data frame with all COSMIC signatures (which is then used as input for prepare_data)
    for ix in empirical_sub.index:      # save the generated signature weights in a data frame
        for sig in empirical_sub.columns:
            contribs.loc['S{}'.format(ix), sig] = empirical_sub.loc[ix, sig]
    return contribs


# introduce differences in signature which_sig between odd and even samples
def introduce_weight_differences(contribs, which_sig, difference_magnitude):
    for n, ix in enumerate(contribs.index):
        if n % 2 == 1:  # odd samples have higher weights
            contribs.loc[ix, which_sig] *= (1 + difference_magnitude)
        else:
            contribs.loc[ix, which_sig] /= (1 + difference_magnitude)
        contribs.loc[ix] /= contribs.loc[ix].sum()  # re-normalize the weights
    return contribs
