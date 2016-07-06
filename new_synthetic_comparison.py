import os
import sys
import csv
import sqlite3
import numpy as np
import pandas as pd
from itertools import combinations
from utils.wildtype import WILDTYPE_PRO
from utils.mut_freqs import DB_LOCATION

MUT_PAIRS = ['30-88', '54-82', '73-90', '46-82', '24-74', '35-36', '69-84',
             '24-46', '24-82', '13-33', '10-93', '12-19', '33-66', '10-46',
             '32-82', '24-64', '37-63', '33-60', '41-93', '30-35', '35-88',
             '32-46', '20-62', '63-93']
MUT_POSITIONS = list(set([s.split('-')[0] for s in MUT_PAIRS] + [s.split('-')[1] for s in MUT_PAIRS]))
MUT_POSITIONS = sorted(map(int, MUT_POSITIONS))

seq2ind = dict(zip('.MW', range(-1, 2)))


"""
1. Get Counts for each sample. Whole Sample. All reads. Save as file.
2. Calculate marginals from counts for each sample.
3. Aggregate counts and get marginals for total.
4. Get bounds for each sample and aggregate.
"""

def main():
    import glob
    counts = []
    uppers = []
    lowers = []
    bis = []
    read_files = glob.glob('/u2/scripps/samples/*/*pro_reads')
    for f in read_files:
        try:
            bi_file = os.sep.join(f.split(os.sep)[:-1] + ['bi'])
            bi = np.loadtxt(bi_file)
            bis.append(bi)
            upper_file = os.sep.join(f.split(os.sep)[:-1] + ['upper'])
            uppers.append(np.loadtxt(upper_file))
            lower_file = os.sep.join(f.split(os.sep)[:-1] + ['lower'])
            lowers.append(np.loadtxt(lower_file))

            N_reads = sum(1 for line in open(f))
            count_file = os.sep.join(f.split(os.sep)[:-1] + ['bicounts'])
            np.savetxt(count_file, bi * N_reads)
            counts.append(bi * N_reads)
        except IOError:
            pass

    agg_counts = np.dstack(counts).sum(axis=-1)
    np.savetxt('agg_counts', agg_counts)

    agg_counts = agg_counts.reshape(-1, 2, 2)
    agg_bi = agg_counts / agg_counts.sum(axis=2).sum(axis=1)[:,None,None]

    N = len(MUT_POSITIONS)
    bi_frame = np.zeros((N, N, 2, 2))
    bi_frame[np.triu_indices(N, 1)] = agg_bi
    selected_bi = bi_frame[zip(*[(i,j) for i,j in combinations(range(N), 2) if j==i+1])]
    agg_uni = selected_bi.sum(axis=-1)
    last_uni = selected_bi[-1].sum(axis=0)
    print selected_bi.shape, agg_uni.shape, last_uni.shape
    agg_uni = np.concatenate((agg_uni, last_uni[None,:]))
    assert np.allclose(agg_uni.sum(axis=1), np.ones(N))

    paired_uni = np.array([(agg_uni[i,0], agg_uni[j,0]) for i,j in combinations(range(N), 2)])

    np.savetxt('agg_uni', agg_uni)
    np.savetxt('agg_bi', agg_bi.reshape(-1, 4))
    avg_bis = np.dstack(bis).mean(axis=-1)
    np.savetxt('avg_bis', avg_bis)

    lower_bound = paired_uni.sum(axis=1) - 1
    lower_bound[lower_bound < 0] = 0.
    estimate_lower_table = np.vstack((lower_bound, paired_uni[:,0] - lower_bound,
        paired_uni[:,1] - lower_bound, 1 + lower_bound - paired_uni.sum(axis=1))).T

    upper_bound = paired_uni.min(axis=1)
    estimate_upper_table = np.vstack((upper_bound, paired_uni[:,0] - upper_bound,
        paired_uni[:,1] - upper_bound, 1 + upper_bound - paired_uni.sum(axis=1))).T

    assert estimate_lower_table.shape == estimate_upper_table.shape

    np.savetxt('agg_lower_table', estimate_lower_table)
    np.savetxt('agg_upper_table', estimate_upper_table)

    avg_lower = np.dstack(lowers).mean(axis=2)
    avg_upper = np.dstack(uppers).mean(axis=2)

    np.savetxt('avg_lower', avg_lower)
    np.savetxt('avg_upper', avg_upper)

if __name__ == "__main__":
    main()
