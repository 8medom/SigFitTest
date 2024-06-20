import numpy as np                                      # for numerics
import pandas as pd                                     # for data structures
import sys                                              # emergency stop
import MS_config as cfg                                 # all global stuff


# generate the mutational catalog (the number of mutations in all contexts) for one sample
def generate_one_sample(rng, sample_name, num_muts, contribs_given):
    if np.abs(contribs_given.sum() - 1) > cfg.EPSILON:
        print('the sum of signature contributions must be one, now it differs by {:.2e}'.format(contribs_given.sum() - 1))
        sys.exit(1)
    if np.any(contribs_given < 0):                                                      # signature contributions must be positive
        print('signature contributions must be non-negative, now the smallest one is {:.2e} in {}'.format(contribs_given.min(), sample_name))
        sys.exit(1)
    context_weights = pd.Series(0, index = cfg.input_sigs.index)                        # compute weights of tri-nucleotide contexts
    for sig in cfg.input_sigs.columns:                                                  # for each context, the weights is a weighted sum over all signatures
        context_weights += contribs_given[sig] * cfg.input_sigs[sig]
    context_weights /= context_weights.sum()                                            # make sure that the weights sum to 1
    counts = np.zeros(96, dtype = int)                                                  # empty count array
    np.add.at(counts, rng.choice(96, p = context_weights, size = num_muts), 1)          # generate mutation positions
    return pd.DataFrame(counts, index = context_weights.index, columns = [sample_name]) # return the count array as a data frame


# generate mutational catalogs with num_muts where signature acitivity is driven by contribs
def prepare_data_from_signature_activity(rng, num_muts, contribs):
    if contribs.ndim == 1:      # all samples have the same signature contributions
        counts = pd.concat([generate_one_sample(rng, 'S{}'.format(i), num_muts, contribs) for i in range(cfg.N_samples)], axis = 1)
    elif contribs.ndim == 2:    # samples have different (e.g., empirically-driven) signature contributions
        counts = pd.concat([generate_one_sample(rng, 'S{}'.format(i), num_muts, contribs.iloc[i]) for i in range(contribs.shape[0])], axis = 1)
    return counts


# generate mutational catalogs with num_muts by subsampling from the (real) mutational catalog real_data
def prepare_data_from_real_data_by_subsampling(rng, num_muts, real_data):
    counts = pd.DataFrame(0, index = real_data.index, columns = real_data.columns)
    for col in real_data.columns:               # go over samples one by one
        tot_muts = real_data[col].sum()
        if num_muts >= tot_muts:                # if too many mutations are required, skip the downsampling and copy all mutations
            counts[col] = real_data[col]
        else:
            all_muts = [ix for ix in real_data.index for n in range(real_data.loc[ix, col])]    # list with the contexts of all mutations
            chosen_muts = rng.choice(all_muts, size = num_muts, replace = False)                # choose the needed number of mutations
            for ix in chosen_muts: counts.loc[ix, col] += 1                                     # increment the corresponding entries
    return counts


# save the created mutational catalogs in various formats that can be used by the evaluated fitting tools
def save_catalogs(counts, info_label = None):
    if info_label == None: counts.to_csv('data/data_for_deconstructSigs.dat', sep = '\t')
    else: counts.to_csv('data/data_for_deconstructSigs_{}.dat'.format(info_label), sep = '\t')
    counts = counts.reindex(index = cfg.index_MutationalPatterns)
    counts.index.name = 'Type'
    if info_label == None: counts.to_csv('data/data_for_MutationalPatterns.dat', sep = '\t')
    else: counts.to_csv('data/data_for_MutationalPatterns_{}.dat'.format(info_label), sep = '\t')
    counts3 = counts.copy()
    new_index = []
    for ix in counts3.index:
        new_index.append(ix[2:5] + ' ' + ix[0] + ix[2] + ix[6])
    counts3.index = new_index
    counts3 = counts3.reindex(index = cfg.index_YAPSA)
    if info_label == None: counts3.to_csv('data/data_for_YAPSA.dat', sep = '\t')
    else: counts3.to_csv('data/data_for_YAPSA_{}.dat'.format(info_label), sep = '\t')
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
    if info_label == None: counts2.to_csv('data/data_for_spss.dat', sep = ',', index = False)
    else: counts2.to_csv('data/data_for_spss_{}.dat'.format(info_label), sep = ',', index = False)
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
    if info_label == None: counts.to_csv('data/data_for_sigLASSO_spectrum.dat', sep = '\t')
    else: counts.to_csv('data/data_for_sigLASSO_spectrum_{}.dat'.format(info_label), sep = '\t')


# generate empirically-driven signature contributions given by empirical_sub
# the main task of this function is to drop the signatures that contribute less than 10 mutations for a given num_muts
# the reason to do that is that they cannot be recovered as we use the threshold of 10 mutations for the resulting estimates
def generate_weights_empirical(num_muts, empirical_sub, cohort_size = cfg.N_samples):
    for ix in empirical_sub.index:      # re-weight the signature contributions to match the required number of mutations
        empirical_sub.loc[ix] *= num_muts / empirical_sub.loc[ix].sum()
    empirical_sub[empirical_sub < 10] = 0   # remove weak signatures (those with the expected number of mutations below 10)
    if (empirical_sub.sum(axis = 1) == 0).sum() > 0:
        print('for num_muts = {} and realization {}, there is a sample with all absolute signature weights below 10'.format(num_muts, rep))
        print('skipping this iteration...')
        return None
    for ix in empirical_sub.index:      # re-weight the signature contributions again to compensate for possible zeros introduced above
        empirical_sub.loc[ix] /= empirical_sub.loc[ix].sum()
    # create a data frame with all COSMIC signatures and their weights (this is then used as input for prepare_data)
    contribs = pd.DataFrame(0, index = ['S{}'.format(x) for x in range(cohort_size)], columns = cfg.input_sigs.columns, dtype = float)
    for ix in empirical_sub.index:      # save the generated signature weights in a data frame
        for sig in empirical_sub.columns:
            contribs.loc['S{}'.format(ix), sig] = empirical_sub.loc[ix, sig]
    return contribs


# introduce differences (quantified by difference_magnitude) in signature which_sig between odd and even samples
def introduce_weight_differences(contribs, which_sig, difference_magnitude):
    for n, ix in enumerate(contribs.index):
        if n % 2 == 1:                              # odd samples have higher weights
            contribs.loc[ix, which_sig] *= (1 + difference_magnitude)
        else:                                       # even samples have lower weights
            contribs.loc[ix, which_sig] /= (1 + difference_magnitude)
        contribs.loc[ix] /= contribs.loc[ix].sum()  # re-normalize the weights
    return contribs
