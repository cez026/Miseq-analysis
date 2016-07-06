import os
import sys
import csv
import glob
import sqlite3
import numpy as np
import pandas as pd
from itertools import combinations, izip_longest
from utils.wildtype import WILDTYPE_PRO
from utils.mut_freqs import DB_LOCATION

#MUT_PAIRS = ['30-88', '54-82', '73-90', '46-82', '24-74', '35-36', '69-84',
#             '24-46', '24-82', '13-33', '10-93', '12-19', '33-66', '10-46',
#             '32-82', '24-64', '37-63', '33-60', '41-93', '30-35', '35-88',
#             '32-46', '20-62', '63-93']
#MUT_POSITIONS = list(set([s.split('-')[0] for s in MUT_PAIRS] + [s.split('-')[1] for s in MUT_PAIRS]))
#MUT_POSITIONS = sorted(map(int, MUT_POSITIONS))

seq2ind = dict(zip('MW.', range(3)))

"""
Given pair(s) of positions, 
1. calculate the bivariate counts in each sample
  - calculate the univariate marginals, bivariate marginals, and upper/lower bounds
2. calculate the aggregate bivariate count
  - calculate aggregate univariate marginals, bivariate marginals, and upper/lower bounds
3. calculate average univariate marginals, bivariate marginals, upper/lower bounds
"""

def get_pair_counts(pair_list, seqfile):
    # assumes we're dealing with fasta
    pair_counts = np.zeros((len(pair_list), 3, 3), dtype=float)
    with open(seqfile, 'r') as fin:
        for name, seq in izip_longest(*[fin]*2):
            seq = seq.strip()
            if not seq: break
            for k, (p1, p2) in enumerate(pair_list):
                i1, i2 = p1 - 1, p2 - 1
                pair_counts[k, seq2ind[seq[i1]], seq2ind[seq[i2]]] += 1
    return pair_counts

def get_bounds(uni):
    lower_bound = max(0, uni[0] + uni[2] - 1)
    lower_table = np.array([lower_bound, uni[0] - lower_bound, uni[2] - lower_bound,
                            1 + lower_bound - uni[0] - uni[2]])
    upper_bound = min(uni[0], uni[2])
    upper_table = np.array([upper_bound, uni[0] - upper_bound, uni[2] - upper_bound,
                            1 + upper_bound - uni[0] - uni[2]])

    assert lower_table[0] <= upper_table[0] or np.allclose(lower_table[0], upper_table[0]), (lower_table[0], upper_table[0])
    assert abs(sum(lower_table) - 1) < 1e-6
    assert abs(sum(upper_table) - 1) < 1e-6

    return lower_table, upper_table

def osjoin(iterable):
    return os.sep.join(iterable)

def sample_main(pair_list, remake=False):
    aggregate_counts, lower_tables, upper_tables = [], [], []
    for f in glob.glob('/u2/scripps/samples/*/*_pro_reads_MW.fasta'):
        basepath = osjoin(f.split(os.sep)[:-1])
        print basepath.split(os.sep)[-1]
        if os.path.exists(osjoin((basepath, 'skip'))): continue
        count_files = [osjoin((basepath, 'pair_counts_%i-%i'%pair)) for pair in pair_list]
        bimarg_files = [osjoin((basepath, 'pair_bimarg_%i-%i'%pair)) for pair in pair_list]
        unimarg_files = [osjoin((basepath, 'pair_unimarg_%i-%i'%pair)) for pair in pair_list]
        lower_files = [osjoin((basepath, 'pair_lower_%i-%i'%pair)) for pair in pair_list]
        upper_files = [osjoin((basepath, 'pair_upper_%i-%i'%pair)) for pair in pair_list]

        # get counts
        if not all(map(os.path.exists, count_files)) or remake:
            bivariate_counts = get_pair_counts(pair_list, f)
            # get rid of '.' counts
            bivariate_counts = bivariate_counts[:, :-1, :-1]
            for k, pair in enumerate(pair_list):
                np.savetxt(count_files[k], bivariate_counts[k].reshape(1,4))
        else:
            bivariate_counts = np.array(map(np.loadtxt, count_files)).reshape(-1, 2, 2)
        aggregate_counts.append(bivariate_counts)

        # get marginals
        if not all(map(os.path.exists, bimarg_files + unimarg_files)) or remake:
            bivariate_marginals = bivariate_counts / bivariate_counts.sum(axis=-1).sum(axis=-1)[:, None, None]
            univariate_marginals = np.hstack((bivariate_marginals.sum(axis=-1), bivariate_marginals.sum(axis=-2)))

            assert np.allclose(bivariate_marginals.sum(axis=-1).sum(axis=-1), np.ones(len(pair_list))), bivariate_marginals.sum(axis=-1).sum(axis=-1)
            assert np.allclose(univariate_marginals[:,:2].sum(axis=-1), np.ones(len(pair_list)))
            assert np.allclose(univariate_marginals[:,2:].sum(axis=-1), np.ones(len(pair_list)))
            for k, pair in enumerate(pair_list):
                np.savetxt(bimarg_files[k], bivariate_marginals[k].reshape(1,4))
                np.savetxt(unimarg_files[k], univariate_marginals[k].reshape(1,4))
        else:
            bivariate_marginals = np.array(map(np.loadtxt, bimarg_files)).reshape(-1, 2, 2)
            univariate_marginals = np.array(map(np.loadtxt, unimarg_files)).reshape(-1, 2, 2)
        
        # get bounds
        if not all(map(os.path.exists, lower_files + upper_files)) or remake:
            lowers = []
            uppers = []
            for k, pair in enumerate(pair_list):
                lower, upper = get_bounds(univariate_marginals[k])
                np.savetxt(lower_files[k], lower.reshape(1,4))
                np.savetxt(upper_files[k], upper.reshape(1,4))
                lowers.append(lower)
                uppers.append(upper)
            lowers, uppers = map(np.array, (lowers, uppers))
        else:
            lowers = np.array(map(np.loadtxt, lower_files)).reshape(-1, 2, 2)
            uppers = np.array(map(np.loadtxt, upper_files)).reshape(-1, 2, 2)
        lower_tables.append(lowers)
        upper_tables.append(uppers)
        # finished with this file. move to next one

    #aggregate_counts = np.concatenate([x[None,...] for x in aggregate_counts], axis=0)
    aggregate_counts = np.array(aggregate_counts)
    lower_tables, upper_tables = map(np.array, (lower_tables, upper_tables))
    assert lower_tables.shape[0] == upper_tables.shape[0]
    
    return aggregate_counts, lower_tables, upper_tables

def aggregate_main(pair_list, aggregate_counts, lower_tables, upper_tables):
    # aggregate marginals
    aggregate_summed = aggregate_counts.sum(axis=0)
    aggregate_bivariates = aggregate_summed / aggregate_summed.sum(axis=-1).sum(axis=-1)[...,None,None]
    aggregate_univariates = np.hstack((aggregate_bivariates.sum(axis=-1), aggregate_bivariates.sum(axis=-2)))
    for k, pair in enumerate(pair_list):
        np.savetxt('agg_bimarg_%i-%i'%pair, aggregate_bivariates[k].reshape(1,4))
        np.savetxt('agg_unimarg_%i-%i'%pair, aggregate_univariates[k].reshape(1,4))

    # aggregate bounds
    for k, pair in enumerate(pair_list):
        lower, upper = get_bounds(aggregate_univariates[k])
        np.savetxt('agg_lower_%i-%i'%pair, lower.reshape(1,4))
        np.savetxt('agg_upper_%i-%i'%pair, upper.reshape(1,4))

    # average marginals
    print aggregate_counts.shape
    average_bivariates = aggregate_counts / aggregate_counts.sum(axis=-1).sum(axis=-1)[...,None,None]
    average_bivariates = average_bivariates.mean(axis=0)
    print average_bivariates.shape
    average_univariates = np.hstack((average_bivariates.sum(axis=-1), average_bivariates.sum(axis=-2)))
    print average_univariates.shape
    for k, pair in enumerate(pair_list):
        np.savetxt('avg_bimarg_%i-%i'%pair, average_bivariates[k].reshape(1,4))
        np.savetxt('avg_unimarg_%i-%i'%pair, average_univariates[k].reshape(1,4))

    average_lower = lower_tables.mean(axis=0)
    average_upper = upper_tables.mean(axis=0)
    for k, pair in enumerate(pair_list):
        np.savetxt('avg_lower_%i-%i'%pair, average_lower[k].reshape(1,4))
        np.savetxt('avg_upper_%i-%i'%pair, average_upper[k].reshape(1,4))

MUT_PAIRS = ['30-88', '54-82', '73-90', '46-82', '24-74', '35-36', '69-84',
             '24-46', '24-82', '13-33', '10-93', '12-19', '33-66', '10-46',
             '32-82', '24-64', '37-63', '33-60', '41-93', '30-35', '35-88',
             '32-46', '20-62', '63-93']
def main():
    pair_list = [tuple(map(int, p.split('-'))) for p in MUT_PAIRS]
    print pair_list
    #pair_list = [(30, 88), (54, 82), (73, 90)]
    args = sample_main(pair_list, remake=True)
    aggregate_main(pair_list, *args)

if __name__ == "__main__":
    main()
