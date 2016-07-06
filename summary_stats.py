#!/usr/bin/env python
"""
Created April 15, 2013

@author: Bill
"""

import os, sys
import csv
import collections
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
import matplotlib.gridspec as gridspec
from itertools import izip_longest
from collections import defaultdict
from read_mutation_files import get_mutations_aa
from utils.wildtype import WILDTYPE_PRO, WILDTYPE_GAG, AA_TO_IND, IND_TO_AA

POLYMORPHICS = [35, 37, 41, 63, 77, 93]
GAG_CS = [132, 363, 377, 432, 448]
FULL_CS = [0] + GAG_CS + [499]
GAG_LABELS = ['MA', 'CA', 'p2', 'NC', 'p1', 'p6']
GAG_LABEL_POS = [(FULL_CS[i+1] + FULL_CS[i]) / 2.
                 for i in range(len(GAG_CS) + 1)]
PRO_SIG_POS = [10, 20, 24, 30, 32, 36, 46, 48, 53, 54, 62, 63, 71, 73, 74, 82,
               88, 90, 93]
GAG_SIG_POS = [12, 62, 75, 76, 79, 81, 112, 128, 132, 200, 219, 360, 362, 363,
               368, 369, 370, 371, 373, 374, 375, 376, 381, 389, 390, 401, 409,
               428, 430, 431, 436, 437, 449, 451, 452, 453, 468, 474, 484, 487,
               497]

def get_series(directory, gag=False):
    samples = {}
    with open('/u2/scripps/samples/patient_groups-failures.txt', 'r') as f:
        current_patient = ''
        for line in f:
            line = line.strip().split()
            if line[0] == 'PATIENT-SAMPLE-ID': continue

            patient, sample = line[0].split('-')
            if patient != current_patient:
                samples[patient] = [line[0]]
                current_patient = patient
            else:
                samples[patient].append(line[0])

    series = []
    parent = 'Sample_%s'
    suffix = '-refined_pro_counts_aa'
    if gag: suffix = '-refined_gag_counts_aa'
    for key in samples:
        #if len(samples[key]) < 2: continue

        sample_filenames = []
        for sample in sorted(samples[key]):
            for filename in os.listdir(os.path.join(directory, parent%sample)):
                if filename.endswith(suffix):
                    full = os.path.join(directory, parent%sample, filename)
                    sample_filenames.append(full)
        #if len(sample_filenames) < 2:
        #    continue
        series.append(sample_filenames)

    return series

def per_position_stats(freqs, outfile, mono_threshold=.90, gag=False):
    L = 99
    WILDTYPE = WILDTYPE_PRO
    if gag:
        L = 499
        WILDTYPE = WILDTYPE_GAG

    N_patients = len(freqs)
    flips = [0 for i in range(L)]
    reversions = [0 for i in range(L)]
    monomorphic_patients = [0 for i in range(L)]
    not_wt = [0 for i in range(L)]
    monomorphic_samples = [0 for i in range(L)]
    sample_count = 0
    poly_count = 0
    poly_arr = []
    for i in range(N_patients):
        patient_samples = np.array(freqs[i])
        sample_count += len(patient_samples)
        for j in range(L):
            monomorphic = 1
            last_wildtype = 0
            last_AA = ''
            for k in range(patient_samples.shape[0]):
                sample = patient_samples[k, j]
                ind2, ind1 = sample.argsort()[-2:]

                # Monomorphic or not
                if sample[ind1] < mono_threshold:
                    monomorphic = 0
                else:
                    monomorphic_samples[j] += 1

                if k == 0:
                    last_AA = IND_TO_AA[ind1]

                is_wildtype = (IND_TO_AA[ind1] == WILDTYPE[j])
                if not is_wildtype:
                    not_wt[j] += 1
                if k == 0:
                    last_wildtype = is_wildtype
                else:
                    #if last_wildtype != is_wildtype:
                    if last_AA != IND_TO_AA[ind1]:
                        if sample[ind1] < .6:
                            poly_count += 1
                            poly_arr.append(j+1)
                        flips[j] += 1
                        if is_wildtype:#not last_wildtype:
                            reversions[j] += 1
                last_wildtype = IND_TO_AA[ind1]
            if monomorphic:
                monomorphic_patients[j] += 1
    monomorphic_patients = [m/float(N_patients) for m in monomorphic_patients]
    print sample_count
    print poly_count
    print len(set(poly_arr))
    d = defaultdict(int)
    for i in poly_arr:
        d[i] += 1
    print sorted(d.items(), key=lambda x: x[1])

    with open(outfile, 'wb') as csvfile:
        output_writer = csv.writer(csvfile, dialect='excel')
        to_write = zip(range(1, L+1), flips, reversions, not_wt,
                             [sample_count]*L, monomorphic_samples,
                             monomorphic_patients)
        output_writer.writerow(['Position', 'Flips', 'Reversions', 'Not WT',
                                'Sample Count', 'Monomorphic Samples',
                                'Monomorphic Patients'])
        output_writer.writerows(to_write)


def per_patient_stats(freqs, outfile, mono_threshold=.98, gag=False):

    PRO_POLY_MUTS = [10, 12, 13, 14, 15, 35, 36, 37, 41, 62, 63, 64, 77,
                    93]
    GAG_POLY_MUTS = [362, 370, 371, 373]
    GAG_SIG_MUTS = list(set([128, 132, 360, 362, 363, 368, 373, 374, 375, 376,
                             381, 428, 430, 431, 436, 437, 449, 451, 452, 453,
                             12, 62, 75, 76, 79, 81, 112, 200, 219, 369, 370,
                             371, 389, 390, 401, 409, 468, 474, 487, 497]))
    PRO_SIG_MUTS = [10, 11, 13, 20, 22, 23, 24, 30, 32, 33, 34, 35, 36, 43,
                    45, 46, 47, 48, 50, 53, 54, 55, 58, 60, 62, 63, 66, 71,
                    72, 73, 74, 75, 76, 77, 79, 82, 83, 84, 85, 88, 89, 90,
                    92, 93, 95]
    PRO_POLY_ONLY_MUTS = list(set(PRO_POLY_MUTS) - set(PRO_SIG_MUTS))
    PRO_SIG_ONLY_MUTS = list(set(PRO_SIG_MUTS) - set(PRO_POLY_MUTS))
    PRO_SIG_POLY_MUTS = list(set(PRO_POLY_MUTS) & set(PRO_SIG_MUTS))
    PRO_NON_ASSOC = list(set(range(99)) - set(PRO_POLY_MUTS) - set(PRO_SIG_MUTS))
    GAG_POLY_ONLY_MUTS = list(set(GAG_POLY_MUTS) - set(GAG_SIG_MUTS))
    GAG_SIG_ONLY_MUTS = list(set(GAG_SIG_MUTS) - set(GAG_POLY_MUTS))
    GAG_SIG_POLY_MUTS = list(set(GAG_POLY_MUTS) & set(GAG_SIG_MUTS))
    GAG_NON_ASSOC = list(set(range(499)) - set(GAG_POLY_MUTS) - set(GAG_SIG_MUTS))

    print len(GAG_SIG_ONLY_MUTS)

    L = 99
    WILDTYPE = WILDTYPE_PRO
    if gag:
        L = 499
        WILDTYPE = WILDTYPE_GAG

    N_patients = len(freqs)
    flips = []
    reversions = []
    monomorphic_positions = []
    N_samples = []
    N_muts = []

    poly_stats = []
    sig_stats = []
    both_stats = []

    for i in range(N_patients):
        N_flips = 0
        N_reverts = 0
        N_pos_monomorphic = 0
        patient_samples = np.array(freqs[i]) #shape = (N_samples, 99, 22)
        N_muts_patient = [0 for el in range(L)]
        N_poly, N_sig, N_both = [0, 0], [0,0], [0,0]
        for j in range(L):
            monomorphic = 1
            last_wildtype = 0
            last_AA = ''
            for k in range(patient_samples.shape[0]):
                sample = patient_samples[k, j]
                ind2, ind1 = sample.argsort()[-2:]

                # Monomorphic or not
                if sample[ind1] < mono_threshold:
                    monomorphic = 0

                is_wildtype = (IND_TO_AA[ind1] == WILDTYPE[j])
                if not is_wildtype:
                    N_muts_patient[j] = 1
                if k == 0:
                    last_AA = IND_TO_AA[ind1]
                # Set first sample as wildtype or not
                if k == 0:
                    last_wildtype = is_wildtype
                else:
                    if last_AA != IND_TO_AA[ind1]:#last_wildtype != is_wildtype:
                        N_flips += 1

                        if not gag:
                            if j+1 in PRO_POLY_ONLY_MUTS:
                                N_poly[0] += 1
                            elif j+1 in PRO_SIG_ONLY_MUTS:
                                N_sig[0] += 1
                            elif j+1 in PRO_SIG_POLY_MUTS:
                                N_both[0] += 1
                        else:
                            if j+1 in GAG_POLY_ONLY_MUTS:
                                N_poly[0] += 1
                            elif j+1 in GAG_SIG_ONLY_MUTS:
                                N_sig[0] += 1
                            elif j+1 in GAG_SIG_POLY_MUTS:
                                N_both[0] += 1

                        if is_wildtype:#not last_wildtype:
                            N_reverts += 1

                            if not gag:
                                if j+1 in PRO_POLY_ONLY_MUTS:
                                    N_poly[1] += 1
                                elif j+1 in PRO_SIG_ONLY_MUTS:
                                    N_sig[1] += 1
                                elif j+1 in PRO_SIG_POLY_MUTS:
                                    N_both[1] += 1
                            else:
                                if j+1 in GAG_POLY_ONLY_MUTS:
                                    N_poly[1] += 1
                                elif j+1 in GAG_SIG_ONLY_MUTS:
                                    N_sig[1] += 1
                                elif j+1 in GAG_SIG_POLY_MUTS:
                                    N_both[1] += 1

            if monomorphic:
                N_pos_monomorphic += 1
            #if j+1 in POLYMORPHICS and not monomorphic:
            #    N_pos_monomorphic += 1

        flips.append(N_flips)
        reversions.append(N_reverts)
        monomorphic_positions.append(N_pos_monomorphic/float(L))
        N_samples.append(patient_samples.shape[0])
        N_muts.append(sum(N_muts_patient))

        poly_stats.append(N_poly)
        sig_stats.append(N_sig)
        both_stats.append(N_both)

    #print (np.array(N_muts)/np.array(N_samples)).mean()
    #print sum(N_muts)
    #print sum(N_samples)

    poly_stats, sig_stats, both_stats = map(np.array, [poly_stats, sig_stats, both_stats])
    stats_func = lambda x: (x[:, 0] / x[:, 1]).mean()
    np.seterr(all='ignore')

    print "poly", poly_stats.mean(axis=0), stats_func(poly_stats)
    print "PI", sig_stats.mean(axis=0), stats_func(poly_stats)
    print "both", both_stats.mean(axis=0), stats_func(both_stats)
    print 'poly', 'pi', 'both'
    print poly_stats.sum(0) + sig_stats.sum(0) + both_stats.sum(0)
    print poly_stats.sum(0), sig_stats.sum(0), both_stats.sum(0)
    print sum(flips), sum(reversions)

    with open(outfile, 'wb') as csvfile:
        output_writer = csv.writer(csvfile, dialect='excel')
        to_write = zip(N_samples, flips, reversions, N_muts,
                       monomorphic_positions)
        output_writer.writerow(['Number of Samples', 'Flips', 'Reversions',
                                'Muts in Patient', 'Monomorphic Positions'])
        output_writer.writerows(to_write)

def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def plot_position(position_file, n_dist, gag=False):
    if gag:
        GAG_CS_MUTS = [128, 132, 360, 362, 363, 368, 373, 374, 375, 376,
                       381, 428, 430, 431, 436, 437, 449, 451, 452, 453]
        GAG_NONCS_MUTS = [12, 62, 75, 76, 79, 81, 112, 200, 219, 369, 370,
                          371, 389, 390, 401, 409, 468, 474, 487, 497]
        GAG_OTHER = list(set(range(1, 500)) - set(GAG_CS_MUTS) - set(GAG_NONCS_MUTS))
        GAG_PI_ASSOC = [125, 102, 252, 28, 93, 479, 380, 340]
        GAG_POLYMORP = [76, 418, 403, 12, 456, 79, 389]
        GAG_OTHER_INT = [15, 82, 248, 478]
    data = []
    with open(position_file, 'r') as f:
        f.next()
        for line in f:
            line = line.strip().split(',')
            x = map(int, line[:-1]) + [float(line[-1])]
            data.append(x)
    data = np.array(data)
    L = data.shape[0]

    norm = 0
    mut_norm = 0
    for key in n_dist:
        norm += (key-1) * n_dist[key]
        mut_norm += key * n_dist[key]
    norm = float(norm)
    mut_norm = float(mut_norm)
    print data[:,3].mean()

    revs = data[:,2]
    revs /= np.array(map(float, data[:,1]))
    revs[np.isnan(revs)] = 0

    flips = flatten([[i]*data[i,1] for i in range(L) if data[i, 1] > 0])
    flips = np.array([f for f in flips])
    #reversions = flatten([[i]*data[i,2] for i in range(L) if data[i,2] > 0])
    reversions = flatten([[i]*revs[i] for i in range(L) if revs[i] > 0])
    reversions = np.array([r for r in reversions])
    muts = flatten([[i]*data[i, 3] for i in range(L) if data[i, 3] > 0])
    muts = np.array([m for m in muts])

    fig = plt.figure(figsize=(28, 10))
    #if gag:
    #    fig.suptitle('Statistics across GAG per position', size=20)
    #else:
    #    fig.suptitle('Statistics across PRO per position', size=20)
    gs = gridspec.GridSpec(3, 2)
    gs.update(left=.06, right=.94, top=.98, bottom=.05, wspace=.09)
    ax1 = plt.subplot(gs[0, :])
    ax12 = ax1.twinx()
    ax13 = ax1.twiny()
    ax2 = plt.subplot(gs[1, :])
    ax22 = ax2.twinx()
    ax23 = ax2.twiny()
    ax3 = plt.subplot(gs[2, :])
    ax32 = ax3.twinx()
    ax33 = ax3.twiny()
    #ax3 = plt.subplot(gs[2, 0])
    #ax4 = plt.subplot(gs[2, 1])

    if gag:
        #for max, ax in zip((72, 9., 1.09), (ax1, ax2, ax3)):
        for max, ax in zip((84, 16.5, 1.11), (ax1, ax2, ax3)):
            for cs in GAG_CS:
                ax.axvline(x=(cs-.5), linewidth=2, color='k', ls='--')
            #for pos in GAG_PI_ASSOC:
            #    ax.axvline(x=(pos-1), linewidth=1, color='green')
            #for pos in GAG_POLYMORP:
            #    ax.axvline(x=(pos-1), linewidth=1, color='m')
            #for pos in GAG_OTHER_INT:
            #    ax.axvline(x=(pos-1), linewidth=1, color='Orange')
            #max = np.array(trans).max()
            for i in range(len(GAG_LABELS)):
                ax.text(GAG_LABEL_POS[i]-1, max, GAG_LABELS[i], color='k',
                        horizontalalignment='center', size='large')


    sig = np.array(GAG_SIG_POS)-1
    sig_func = lambda x: np.array([x[i] if i in sig else 0
                                   for i in range(data.shape[0])])
    params = dict(bins=L, range=(0,L), align='left')
    #mutn, _, _ = ax1.hist(muts, color='g', label='Mutations', **params)
    #trans, _, _ = ax2.hist(flips, color='b', label='Transitions', **params)
    #ax3.hist(reversions, color='c', label='Reversions', **params)
    ax1.bar(range(499), data[:,3], color='.5', label="Mutations", align='center', width=1)
    ax2.bar(range(499), data[:,1], color='.5', label="Transitions", align='center', width=1)
    ax3.bar(range(499), revs, color='.5', label="Reversion Rate", align='center', width=1)
    trans = range(499)

    pi = np.array(GAG_PI_ASSOC)-1
    poly = np.array(GAG_POLYMORP)-1
    #ax12.bar(pi, data[:,3][pi], color='r', align='center', width=1, label="Mutations")
    #ax22.bar(pi, data[:,1][pi], color='r', align='center', width=1, label="Transitions")
    #ax32.bar(pi, revs[pi], color='r', align='center', width=1, label="Reversion Rate")
    #ax12.bar(poly, data[:,3][poly], color='b', align='center', width=1)
    #ax22.bar(poly, data[:,1][poly], color='b', align='center', width=1)
    #ax32.bar(poly, revs[poly], color='b', align='center', width=1)
    #ax1.legend(loc=2)
    #ax2.legend(loc=2)
    #ax3.legend(loc=2)

    tick = lambda x: ticker.FuncFormatter(lambda y, pos: ('%.1f%%')%(y*100/x))
    for ax in (ax1, ax2, ax3):
        ax.set_xlim(0, len(trans))
        ax.set_xticks(np.arange(4, len(trans), 5))
        ax.set_xticklabels(np.arange(5, len(trans), 5))
        if gag:
            ax.set_xlim(0, len(trans))
            ax.set_xticks(np.arange(19, len(trans), 20))
            ax.set_xticklabels(np.arange(20, len(trans), 20))
        ax.grid()
    ax1.set_xlabel('Sequence Position')
    ax2.set_xlabel('Sequence Position')
    ax3.set_xlabel('Sequence Position')

    ax1.set_ylabel('Number of Samples in\nwhich a Mutation Fixes')
    ax12.set_yticks(ax1.get_yticks())
    ax12.yaxis.set_major_formatter(tick(mut_norm))
    ax12.set_ylabel('Percentage of %i Total Samples'%mut_norm)
    ax2.set_ylabel('Number of Transitions')
    ax22.set_yticks(ax2.get_yticks())
    ax22.yaxis.set_major_formatter(tick(norm))
    ax22.set_ylabel('Percentage of %i Possible Transitions'%norm)
    ax3.set_ylabel('Reversion Rate')
    ax32.set_ylabel('Reversion Rate')
    ax3.set_yticks(np.arange(0, 1.2, .2))
    ax32.set_yticks(np.arange(0, 1.2, .2))
    ax3.set_ylim(0, 1.2)
    ax32.set_ylim(0, 1.2)
    ax3.yaxis.set_major_formatter(tick(1))
    ax32.yaxis.set_major_formatter(tick(1))
    #ax3.set_ylabel('Number of Reversions')
    #ax3.set_yticks(ax2.get_yticks())
    #ax32.set_yticks(ax2.get_yticks())
    #ax32.yaxis.set_major_formatter(tick(norm))
    #ax32.set_ylabel('Percentage of %i Possible Reversions'%norm)

    for i, ax in enumerate((ax13, ax23, ax33)):
        ax.set_xlim(0, len(trans))
        #if gag:
        #    ax.set_xticks(np.array(GAG_OTHER_INT)-1)
        #if not gag:
        #    ax.set_xticks(np.array(PRO_SIG_POS)-1)
        #else:
        #    ax.set_xticks(np.array(GAG_SIG_POS)-1)
        #ax.tick_params(axis='x', colors='red', width=4, length=10)
        ax.xaxis.tick_bottom()
        ax.set_xticklabels([])
    plt.show()
    #fig.savefig('poster_plot.png', dpi=300)

    import sys
    sys.exit()
    if gag:
        f = plt.figure()
        ax = f.add_subplot(111)
        mut_list = []
        for pos in range(499):
            tup = (pos, data[pos-1,3], data[pos-1,2], data[pos-1,1], data[pos-1,2]/data[pos-1,1])
            if tup[-1] > 1 or np.isnan(tup[-1]):
                continue
            #if tup[3] < 2: continue
            mut_list.append(tup)
        mut_list = sorted(mut_list, key=lambda x: x[1])[::-1]
        for el in mut_list:
            print "%i\t%.0f\t%.0f\t%.0f\t%.2f"%el

        mut_list = np.array(mut_list)
        vector = np.array([0, 1, 0, 0, 1], dtype=bool)
        mut_list = mut_list[:, vector]
        ax.scatter(mut_list[:,1], mut_list[:,0])
        ax.grid()
        ax.set_xticks(np.arange(0., 1.05, .1))
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(0, 100)
        plt.show()
    else:
        f = plt.figure()
        ax = f.add_subplot(111)
        mut_list = []
        for pos in range(99):
            tup = (pos, data[pos-1,3], data[pos-1,2], data[pos-1,1], data[pos-1,2]/data[pos-1,1])
            if tup[-1] > 1 or np.isnan(tup[-1]):
                continue
            #if tup[3] < 2: continue
            mut_list.append(tup)
        mut_list = sorted(mut_list, key=lambda x: x[1])[::-1]
        for el in mut_list:
            print "%i\t%.0f\t%.0f\t%.0f\t%.2f"%el

        mut_list = np.array(mut_list)
        vector = np.array([0, 1, 0, 0, 1], dtype=bool)
        mut_list = mut_list[:, vector]
        ax.scatter(mut_list[:,1], mut_list[:,0])
        ax.grid()
        ax.set_xticks(np.arange(0., 1.05, .1))
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(0, 100)
        plt.show()

    #params = dict(color='b', range=(0, 16), histtype='stepfilled',
    #              align='mid')
    #ax3.hist(np.array(trans), bins=16, **params)
    #ax32 = ax3.twinx()
    #ax3.set_ylim(0, L)
    #ax32.set_yticks(ax3.get_yticks())
    #ax32.yaxis.set_major_formatter(tick(len(mutn)))
    #ax32.hist(np.array(trans), range=(0, 16), bins=16, normed=100, cumulative=True, histtype='step')
    #ax4.hist(mutn, color='g', bins=85, range=(0, 85))
    #ax3.hist(np.array(n1-n2)/norm, facecolor='b', **params)
    #ax4.hist(np.array(n2)/norm, facecolor='g', **params)

    #for ax in (ax3, ax4):
    #    ax.yaxis.set_major_formatter(tick(len(n2)))
    #    ax.set_xlabel('Number of Transitions')
    #    ax.set_ylabel('Percentage of All %i Positions'%len(n2))
    #    ax.set_ylim(0, len(n2)*.85)
    #    ax.set_xticks(np.arange(0, .2, .01))
    #    #ax.set_xticks(np.arange(0, 15,1))
    #    labels = ax.get_xticklabels()
    #    for label in labels: label.set_rotation(30)
    #    ax.set_yticks(np.arange(0, len(n2)*.85, len(n2)*0.05))
    #    ax.grid()

    #ax3.set_title('Histogram of All Transitions')
    #ax4.set_title('Histogram of Mutations')


def plot_patient(patient_file, n_dist, gag=False):
    print patient_file
    data = []
    with open(patient_file, 'r') as f:
        f.next()
        for line in f:
            line = line.strip().split(',')
            x = map(int, line[:-1]) + [float(line[-1])]
            data.append(x)
    data = np.array(data)
    L = data.shape[0]

    N_pos = 99
    if gag: N_pos = 499

    print 'Avg Trans per Pat', data[:,1].mean()
    print 'Avg Rev per Pat', data[:, 2].mean()
    print 'Avg Rev per Pat exclude 0',  np.array([data[i, 2] for i in range(L) if data[i, 2] > 0]).mean()
    print 'Total trans', data[:, 1].sum()
    print 'Total rev', data[:, 2].sum()

    print data[:, 3].mean()
    print data[:, 3].sum()
    muts = flatten([[i]*data[i, 3] for i in range(L) if data[i, 3] > 0])
    muts = np.array([m for m in muts])
    flips = flatten([[i]*data[i,1] for i in range(L) if data[i, 1] > 0])
    flips = np.array([f for f in flips])
    reversions = flatten([[i]*data[i,2] for i in range(L) if data[i,2] > 0])
    reversions = np.array([r for r in reversions])
    series3p = np.array([i for i in range(L) if data[i, 0] > 2])
    series5p = np.array([i for i in range(L) if data[i, 0] > 3])

    norm = 0
    for key in n_dist:
        norm += (key-1) * n_dist[key]
    norm = float(norm)

    fig = plt.figure()
    if gag:
        fig.suptitle('Statistics across GAG per patient', size=20)
    else:
        fig.suptitle('Statistics across PRO per patient', size=20)
    gs = gridspec.GridSpec(3, 3)
    gs.update(left=.04, right=.96, top=.96, bottom=.05, wspace=.10, hspace=.15)
    ax1 = plt.subplot(gs[0, :])
    ax12 = ax1.twinx()
    ax13 = ax1.twiny()
    ax2 = plt.subplot(gs[1, :])
    ax22 = ax2.twinx()
    ax23 = ax2.twiny()
    ax3 = plt.subplot(gs[2, :])
    ax32 = ax3.twinx()
    ax33 = ax3.twiny()
    #ax2 = fig.add_subplot(gs[1, 0])
    #ax3 = fig.add_subplot(gs[1, 1])
    #ax4 = fig.add_subplot(gs[1, 2])

    params = dict(bins=L, range=(0, L), histtype='stepfilled', align='left')

    ax1.hist(muts, color='g', label='Mutations', **params)
    ax2.hist(flips, color='b', label='Transitions', **params)
    ax3.hist(reversions, color='c', label='Reversions', **params)

    tick = lambda x: ticker.FuncFormatter(lambda y, pos: ('%.1f%%')%(y*100/x))

    ax1.set_ylabel('Number of Positions seeing a Mutation')
    ax12.set_ylabel('Percentage of %i Positions'%N_pos)
    ax12.set_yticks(ax1.get_yticks())
    ax12.yaxis.set_major_formatter(tick(N_pos))

    ax2.set_ylabel('Number of Transitions')
    ax22.set_ylabel('Percentage of %i Possible Transitions'%norm)
    ax22.set_yticks(ax2.get_yticks())
    ax22.yaxis.set_major_formatter(tick(norm))

    ax3.set_ylim(ax2.get_ylim())
    ax3.set_yticks(ax2.get_yticks())
    ax3.set_ylabel('Number of Reversions')
    ax32.set_ylabel('Percentage of %i Possible Reversions'%norm)
    ax32.set_yticks(ax2.get_yticks())
    ax32.yaxis.set_major_formatter(tick(norm))

    for ax in (ax1, ax2, ax3):
        ax.legend(loc=2)
        ax.set_xticks(np.arange(0, L, 5))
        ax.set_xlim(0, L)
        ax.set_xlabel('Patient')
        ax.grid()

    for ax in (ax13, ax23, ax33):
        ax.set_xlim(0, L)
        ax.set_xticks(series3p)
        ax.tick_params(axis='x', colors='red', width=2, length=6)
        ax.xaxis.tick_bottom()
        ax.set_xticklabels([])
        fakeax = ax.twiny()
        fakeax.set_xlim(0, L)
        fakeax.set_xticks(series5p)
        fakeax.tick_params(axis='x', colors='red', width=4, length=10)
        fakeax.set_xticklabels([])
        fakeax.xaxis.tick_bottom()

    #params = dict(bins=39, range=(0, .39), histtype='stepfilled')
    #ax2.hist(np.array(n1)/norm, facecolor='r', **params)
    #ax3.hist(np.array(n1-n2)/norm, facecolor='b', **params)
    #ax4.hist(np.array(n2)/norm, facecolor='g', **params)

    #for ax in (ax2, ax3, ax4):
    #    ax.yaxis.set_major_formatter(tick(len(n2)))
    #    ax.set_xlabel('Frequency')
    #    ax.set_ylabel('Percentage of All %i Patients'%len(n2))
    #    ax.set_ylim(0, len(n2)*.85)
    #    ax.set_xticks(np.arange(0, .42, .02))
    #    labels = ax.get_xticklabels()
    #    for label in labels: label.set_rotation(30)
    #    ax.set_yticks(np.arange(0, len(n2)*.85, len(n2)*0.05))
    #    ax.grid()

    #ax1.set_title('Transitions across PRO per patient', size=20)
    #ax2.set_title('Histogram of All Transitions')
    #ax3.set_title('Histogram of Mutations')
    #ax4.set_title('Histogram of Reversions')
    #if gag:
    #    ax1.set_title('Transitions across GAG per patient', size=20)
    plt.show()




def plot_monomorphic(infiles):
    pro_pos_data, pro_pat_data = [], []
    gag_pos_data, gag_pat_data = [], []
    for inf, l in zip(infiles, (pro_pos_data, gag_pos_data,
                                pro_pat_data, gag_pat_data)):
        with open(inf, 'r') as f:
            f.next()
            for line in f:
                line = line.strip().split(',')
                l.append(map(float, line[1::2]))
    pro_pos_data, pro_pat_data = map(np.array, (pro_pos_data, pro_pat_data))
    gag_pos_data, gag_pat_data = map(np.array, (gag_pos_data, gag_pat_data))

    pro_pos_mono = pro_pos_data[:,1]
    pro_pat_mono = pro_pat_data[:,1]
    gag_pos_mono = gag_pos_data[:,1]
    gag_pat_mono = gag_pat_data[:,1]
    no_trans_mono = lambda arr: [arr[i, 1] for i in range(arr.shape[0])
                                 if arr[i,0] == 0]
    pro_pos_tmono = no_trans_mono(pro_pos_data)
    pro_pat_tmono = no_trans_mono(pro_pat_data)
    gag_pos_tmono = no_trans_mono(gag_pos_data)
    gag_pat_tmono = no_trans_mono(gag_pat_data)

    f, axarr = plt.subplots(2, 2)
    f.suptitle('Fraction of Monomorphic Positions/Patients over GAG/PRO',
                size=22)
    params = dict(range=(.3, 1.), bins=70, histtype='stepfilled', alpha=0.5)
    params2 = dict(range=(.3, 1.), bins=70, histtype='step', lw=3, color='k')

    _, _, p1 = axarr[0, 0].hist(pro_pos_mono, **params)
    axarr[0, 1].hist(gag_pos_mono, **params)
    axarr[1, 0].hist(pro_pat_mono, **params)
    axarr[1, 1].hist(gag_pat_mono, **params)
    _, _, p2 = axarr[0, 0].hist(pro_pos_tmono, **params2)
    axarr[0, 1].hist(gag_pos_tmono, **params2)
    axarr[1, 0].hist(pro_pat_tmono, **params2)
    axarr[1, 1].hist(gag_pat_tmono, **params2)
    f.legend((p1, p2), ('All Patients','No Transition Patients'), 'upper right')

    axarr[0, 0].set_title('PRO Positions')
    axarr[0, 1].set_title('GAG Positions')
    axarr[1, 0].set_title('PRO Patients')
    axarr[1, 1].set_title('GAG Patients')
    axarr[0, 0].set_xlabel('Fraction of Patients that are Monomorphic')
    axarr[0, 1].set_xlabel('Fraction of Patients that are Monomorphic')
    axarr[1, 0].set_xlabel('Fraction of Positions that are Monomorphic')
    axarr[1, 1].set_xlabel('Fraction of Positions that are Monomorphic')
    axarr[0, 0].set_ylabel('Fraction of %i Positions'% len(pro_pos_mono))
    axarr[0, 1].set_ylabel('Fraction of %i Positions'% len(gag_pos_mono))
    axarr[1, 0].set_ylabel('Fraction of %i Patients'% len(pro_pat_mono))
    axarr[1, 1].set_ylabel('Fraction of %i Patients'% len(gag_pat_mono))

    tick = lambda x: ticker.FuncFormatter(lambda y, pos: '%.2f'%(y/x))
    for ax, l in zip(axarr.flat, (pro_pos_data, gag_pos_data,
                                  pro_pat_data, gag_pat_data)):
        L = len(l)
        ax.set_xlim(.3, 1)
        ax.set_ylim(0, L*.4)
        ax.set_yticks(np.arange(0, L*.4, L*0.05))
        ax.yaxis.set_major_formatter(tick(L))
        ax.set_xticks(np.arange(.3, 1., 0.05))
        ax.grid()

    f.subplots_adjust(left=0.04, right=0.96, top=0.95, bottom=0.05)
    plt.show()

def plot_reversion_frac(infiles):
    pro_pos_data, pro_pat_data = [], []
    gag_pos_data, gag_pat_data = [], []
    for inf, l in zip(infiles, (pro_pos_data, gag_pos_data,
                                pro_pat_data, gag_pat_data)):
        with open(inf, 'r') as f:
            f.next()
            for line in f:
                line = line.strip().split(',')
                l.append(map(float, line[1:3]))
    pro_pos_data, pro_pat_data = map(np.array, (pro_pos_data, pro_pat_data))
    gag_pos_data, gag_pat_data = map(np.array, (gag_pos_data, gag_pat_data))

    rev_func = lambda arr: [arr[i,1]/arr[i,0] for i in range(arr.shape[0])
                            if arr[i,0] > 0]
    pro_pos_rev = rev_func(pro_pos_data)
    pro_pat_rev = rev_func(pro_pat_data)
    gag_pos_rev = rev_func(gag_pos_data)
    gag_pat_rev = rev_func(gag_pat_data)

    f, axarr = plt.subplots(2, 2)
    f.suptitle('Number of Reversions per Transition for each Position/Patient '
               'over PRO/GAG', size=22)
    params = dict(bins=100, range=(0, 1.), color='g',
                  histtype='stepfilled', alpha=0.7)
    axarr[0, 0].hist(pro_pos_rev, **params)
    axarr[0, 1].hist(gag_pos_rev, **params)
    axarr[1, 0].hist(pro_pat_rev, **params)
    axarr[1, 1].hist(gag_pat_rev, **params)

    axarr[0, 0].set_title('PRO Positions')
    axarr[0, 1].set_title('GAG Positions')
    axarr[1, 0].set_title('PRO Patients')
    axarr[1, 1].set_title('GAG Patients')
    axarr[0, 0].set_ylabel('Fraction of %i Positions seeing transitions'%
                            len(pro_pos_rev))
    axarr[0, 1].set_ylabel('Fraction of %i Positions seeing transitions'%
                            len(gag_pos_rev))
    axarr[1, 0].set_ylabel('Fraction of %i Patients seeing transitions'%
                            len(pro_pat_rev))
    axarr[1, 1].set_ylabel('Fraction of %i Patients seeing transitions'%
                            len(gag_pat_rev))

    tick = lambda x: ticker.FuncFormatter(lambda y, pos: '%.2f'%(y/x))
    for ax, l in zip(axarr.flat,(pro_pos_rev, gag_pos_rev,
                                 pro_pat_rev, gag_pat_rev)):
        L = len(l)
        ax.set_ylim(0, L)
        ax.set_yticks(np.arange(0, L, L*0.05))
        ax.yaxis.set_major_formatter(tick(L))

        ax.set_xticks(np.arange(0, 1., 0.05))
        ax.set_xlabel('Fraction of (Reversions) / (Total Transitions)')
        ax.grid()

    f.subplots_adjust(left=.04, right=.96, bottom=.05, top=.90, wspace=.10)
    plt.show()

def auxillary_stats():
    directory = '/u2/scripps/samples'
    filenames = [os.path.join(root, f) for root, _, files in os.walk(directory)
                 for f in files if f.endswith('refined_pro_counts_aa')]

    thres = 0.9
    muts = []
    monomorphic = []
    n_muts_pro = []
    for filename in filenames:
        sample_mut = []
        sample_monomorphic = []
        sample_n_muts = 0
        arr, counts = get_mutations_aa(filename, 0.01)
        if np.any(counts < 1e-8):
            continue
        for i in range(arr.shape[0]):
            arr[i] /= counts[i]
            max_ind = arr[i].argsort()[-1]
            max_freq = arr[i][max_ind]
            is_not_wt = (WILDTYPE_PRO[i] != IND_TO_AA[max_ind])
            mono = (max_freq > thres)
            if is_not_wt: sample_n_muts += 1
            sample_mut.append(is_not_wt)
            sample_monomorphic.append(mono)

        muts.append(sample_mut)
        monomorphic.append(sample_monomorphic)
        n_muts.append(sample_n_muts)

    muts = np.array(muts)
    monomorphic = np.array(monomorphic)
    print muts.shape
    print monomorphic.shape
    print muts.sum(axis=0).shape
    print muts.sum(axis=0)
    print monomorphic.sum(axis=0)/163.

def plot_statics(infiles):
    pro_pos_data, gag_pos_data = [], []
    for inf, l in zip(infiles, (pro_pos_data, gag_pos_data)):
        with open(inf, 'r') as f:
            f.next()
            for line in f:
                line = line.strip().split(',')
                x = map(int, line[:-1]) + [float(line[-1])]
                l.append(x)
    pro_data, gag_data = map(np.array, (pro_pos_data, gag_pos_data))
    L_pro, L_gag = pro_data.shape[0], gag_data.shape[0]
    pro_samp, gag_samp = map(float, (pro_data[0, 4], gag_data[0, 4]))
    print pro_samp, gag_samp

    pro_muts = pro_data[:, 3]
    gag_muts = gag_data[:, 3]

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    fig.subplots_adjust(left=.04, right=.96, top=.94, bottom=.04, wspace=.10, hspace=.20)

    params1 = dict(bins=100, range=(0, 1), histtype='bar', align='mid')
    params2 = dict(bins=20, range=(0, 1), histtype='bar', align='mid')
    ax1.hist(pro_muts/pro_samp, **params2)
    ax2.hist(gag_muts/gag_samp, **params2)
    ax1.hist(pro_muts/pro_samp, color='c', **params1)
    ax2.hist(gag_muts/gag_samp, color='c', **params1)

    tick = lambda x: ticker.FuncFormatter(lambda y, pos: ('%.0f%%')%(y*100/x))
    for ax, samp, ylim in zip((ax1, ax2), (pro_samp, gag_samp), (99, 499)):
        ax.set_xlim(0, 1)
        ax.set_xticks(np.arange(0, 1, .05))
        ax.xaxis.set_major_formatter(tick(1))
        ax.set_ylim(0, ylim)
        ax.set_yticks(np.arange(0, ylim, ylim/10.))
        ax.yaxis.set_major_formatter(tick(ylim))
        ax.grid()
        ax.set_xlabel('Mutation Rate over %i samples'%samp)
        ax.set_ylabel('Percentage of %i Positions'%ylim)
    ax1.set_title('PRO')
    ax2.set_title('GAG')
    plt.show()

def plot_sample_muts():
    PRO_POLY_MUTS = [10, 12, 13, 14, 15, 35, 36, 37, 41, 62, 63, 64, 77,
                    93]#, 18, 19, 57, 67]
    GAG_POLY_MUTS = [362, 370, 371, 373]
    GAG_SIG_MUTS = list(set([128, 132, 360, 362, 363, 368, 373, 374, 375, 376,
                             381, 428, 430, 431, 436, 437, 449, 451, 452, 453,
                             12, 62, 75, 76, 79, 81, 112, 200, 219, 369, 370,
                             371, 389, 390, 401, 409, 468, 474, 487, 497]))
    PRO_SIG_MUTS = [10, 11, 13, 20, 22, 23, 24, 30, 32, 33, 34, 35, 36, 43,
                    45, 46, 47, 48, 50, 53, 54, 55, 58, 60, 62, 63, 66, 71,
                    72, 73, 74, 75, 76, 77, 79, 82, 83, 84, 85, 88, 89, 90,
                    92, 93, 95]
    PRO_POLY_ONLY_MUTS = list(set(PRO_POLY_MUTS) - set(PRO_SIG_MUTS))
    PRO_SIG_ONLY_MUTS = list(set(PRO_SIG_MUTS) - set(PRO_POLY_MUTS))
    PRO_SIG_POLY_MUTS = list(set(PRO_POLY_MUTS) & set(PRO_SIG_MUTS))
    PRO_NON_ASSOC = list(set(range(1,100)) - set(PRO_POLY_MUTS) - set(PRO_SIG_MUTS))
    GAG_POLY_ONLY_MUTS = list(set(GAG_POLY_MUTS) - set(GAG_SIG_MUTS))
    GAG_SIG_ONLY_MUTS = list(set(GAG_SIG_MUTS) - set(GAG_POLY_MUTS))
    GAG_SIG_POLY_MUTS = list(set(GAG_POLY_MUTS) & set(GAG_SIG_MUTS))
    GAG_NON_ASSOC = list(set(range(1, 500)) - set(GAG_SIG_POLY_MUTS) - set(GAG_SIG_ONLY_MUTS))
    print GAG_POLY_ONLY_MUTS
    print GAG_SIG_ONLY_MUTS
    print GAG_SIG_POLY_MUTS
    N_muts_all_pro = defaultdict(int)
    N_muts_poly_pro = defaultdict(int)
    N_muts_sig_pro = defaultdict(int)
    N_muts_both_pro = defaultdict(int)
    N_muts_all_gag = defaultdict(int)
    N_muts_poly_gag = defaultdict(int)
    N_muts_sig_gag = defaultdict(int)
    N_muts_both_gag = defaultdict(int)
    N_muts_non_gag = defaultdict(int)
    N_muts_non_pro = defaultdict(int)

    directory = '/u2/scripps/samples/'
    pro_samp, gag_samp  = 0, 0
    pro_muts, gag_muts = [], []
    pro_n_muts, gag_n_muts = [], []
    pro_n_muts_41 = []
    pro_n_muts_45 = []
    sig_list_41 = [10, 11, 13, 20, 34, 35, 36, 43, 45, 55, 58, 60, 63, 71, 74,
                   75, 77, 79, 83, 85, 89, 91, 93, 95, 23, 24, 30, 32, 33, 46,
                   47, 48, 50, 53, 54, 73, 76, 82, 84, 88, 90]
    sig_list_45 = sig_list_41 + [22, 62, 66, 72, 92]#, 12, 14, 15, 19, 37, 41, 57, 64, 67, 69, 70]))
    sig_list_45.remove(91)
    for syst in ('pro', 'gag'):
        gag = 0
        WILDTYPE = WILDTYPE_PRO
        if syst == 'gag':
            gag = 1
            WILDTYPE = WILDTYPE_GAG
        filenames = [os.path.join(r, f) for r, _, files in os.walk(directory)
                     for f in files if f.endswith('refined_%s_counts_aa'%syst)]
        #filenames = get_series(directory, gag)
        #filenames = [sample for series in filenames for sample in series]
        muts = []
        n_muts = []
        n_muts_41 = []
        n_muts_45 = []
        samp = 0
        n_dist = defaultdict(int)
        for series in [0]:#filenames:
            #broke = 0
            #if len(series) < 2: continue
            to_add_muts = []
            to_add_n_muts = []
            to_add_n_muts_41 = []
            to_add_n_muts_45 = []
            for z, filename in enumerate(filenames):#enumerate(series):#filenames:
                #if z > 0:
                #    broke = 1
                #    continue
                sample_muts = []
                sample_n_muts = 0
                sample_n_muts_41 = 0
                sample_n_muts_45 = 0
                arr, counts = get_mutations_aa(filename, 0.01, gag=gag)
                if np.any(counts < 1e-8):
                    broke = 1
                    continue
                for i in range(arr.shape[0]):
                    arr[i] /= counts[i]
                    max_ind = arr[i].argsort()[-1]
                    max_freq = arr[i][max_ind]
                    is_not_wt = (WILDTYPE[i] != IND_TO_AA[max_ind])
                    if is_not_wt:
                        sample_n_muts += 1
                        if (i+1) in sig_list_41:
                            sample_n_muts_41 += 1
                        if (i+1) in sig_list_45:
                            sample_n_muts_45 += 1
                        if not gag:
                            N_muts_all_pro[i+1] += 1
                            if (i+1) in PRO_POLY_ONLY_MUTS:
                                N_muts_poly_pro[i+1] += 1
                            if (i+1) in PRO_SIG_ONLY_MUTS:
                                N_muts_sig_pro[i+1] += 1
                            if (i+1) in PRO_SIG_POLY_MUTS:
                                N_muts_both_pro[i+1] += 1
                            if (i+1) in PRO_NON_ASSOC:
                                N_muts_non_pro[i+1] += 1
                        else:
                            N_muts_all_gag[i+1] += 1
                            if (i+1) in GAG_POLY_ONLY_MUTS:
                                N_muts_poly_gag[i+1] += 1
                            if (i+1) in GAG_SIG_ONLY_MUTS:
                                N_muts_sig_gag[i+1] += 1
                            if (i+1) in GAG_SIG_POLY_MUTS:
                                N_muts_both_gag[i+1] += 1
                            if (i+1) in GAG_NON_ASSOC:
                                N_muts_non_gag[i+1] += 1
                    sample_muts.append(is_not_wt)
                to_add_n_muts.append(sample_n_muts)
                to_add_n_muts_41.append(sample_n_muts_41)
                to_add_n_muts_45.append(sample_n_muts_45)
                to_add_muts.append(sample_muts)
            #if broke:
            #    continue
            for i in range(len(to_add_muts)):
                n_muts.append(to_add_n_muts[i])
                n_muts_41.append(to_add_n_muts_41[i])
                n_muts_45.append(to_add_n_muts_45[i])
                muts.append(to_add_muts[i])
            #n_dist[len(series)] += 1
            #samp += len(series)
            samp += 1
        if not gag:
            pro_muts = np.array(muts).sum(axis=0)
            pro_n_muts = n_muts
            pro_n_muts_41 = n_muts_41
            pro_n_muts_45 = n_muts_45
            pro_samp = float(samp)
        else:
            gag_muts = np.array(muts).sum(axis=0)
            gag_n_muts = n_muts
            gag_samp = float(samp)
        #print n_dist.items()
        #print sum([v*k for k,v in n_dist.items()])
        #print len(pro_n_muts)

    print sum([x > 0 for x in N_muts_all_pro])
    print sum([x > 0 for x in N_muts_all_gag])
    print sum([x > 0 for x in N_muts_sig_pro])
    print sum([x > 0 for x in N_muts_sig_gag])
    #import sys
    #sys.exit()

    display_func = lambda x: '%i\t%i\t%.2f'%(len(x), sum(x.values()), np.mean(x.values()))
    print 'N_pos\tN_muts\tAvg Muts'
    print display_func(N_muts_all_pro)
    print display_func(N_muts_sig_pro)
    print display_func(N_muts_both_pro)
    #print N_muts_both_pro.items()
    print display_func(N_muts_poly_pro)
    print display_func(N_muts_non_pro)

    print display_func(N_muts_all_gag)
    print display_func(N_muts_sig_gag)
    print display_func(N_muts_both_gag)
    print display_func(N_muts_non_gag)

    pro_avg_muts_samp = np.array(pro_n_muts).mean()
    gag_avg_muts_samp = np.array(gag_n_muts).mean()
    pro_avg_muts_pos = np.array(pro_muts).mean()
    gag_avg_muts_pos = np.array(gag_muts).mean()
    print len(pro_n_muts), len(gag_n_muts)
    print 'TOTAL'
    print sum(pro_n_muts), sum(pro_muts)
    print sum(gag_n_muts), sum(gag_muts)

    print 'PRO/GAG AVG/SAMP'
    print pro_avg_muts_samp, gag_avg_muts_samp
    print 'PRO/GAG AVG/POS'
    print pro_avg_muts_pos, gag_avg_muts_pos

    pro_d = defaultdict(int)
    gag_d = defaultdict(int)
    pro_d_41 = defaultdict(int)
    pro_d_45 = defaultdict(int)
    for i in pro_n_muts:
        pro_d[i] += 1
    for i in gag_n_muts:
        gag_d[i] += 1
    for i in pro_n_muts_41:
        pro_d_41[i] += 1
    for i in pro_n_muts_45:
        pro_d_45[i] += 1
    pro_x = np.array(sorted(pro_d.items()))[:, 0]
    gag_x = np.array(sorted(gag_d.items()))[:, 0]
    pro_x_41 = np.array(sorted(pro_d_41.items()))[:, 0]
    pro_x_45 = np.array(sorted(pro_d_45.items()))[:, 0]
    p_n_muts = np.array(sorted(pro_d.items()))[:,1]
    g_n_muts = np.array(sorted(gag_d.items()))[:,1]
    p_n_muts_41 = np.array(sorted(pro_d_41.items()))[:, 1]
    p_n_muts_45 = np.array(sorted(pro_d_45.items()))[:, 1]
    pro_samp = len(pro_n_muts)

    fig = plt.figure()
    rc('font', size=16)
    fig.suptitle('Figure 1. Number of Dominant Mutations across 41-Drug Associated Protease Positions', size=26)
    ax1 = fig.add_subplot(111)
    #ax2 = fig.add_subplot(122)
    #ax3 = fig.add_subplot(223)
    #ax4 = fig.add_subplot(224)
    fig.subplots_adjust(left=.06, right=.96, top=.94, bottom=.06, wspace=.10, hspace=.20)

    omar_x = range(19)
    omar_naive = [.07, .18, .25, .225, .15, .07, .02, .01, .005] + [0]*10
    omar_treated = [.02, .04, .07, .06, .07, .06, .08, .09, .12, .11, .10, .07,
                    .05, .03, .02, .01] + [0]*3
    omar_1 = [.025, 0.1, 0.155, 0.15, 0.13, 0.12, 0.095, 0.075, 0.065, 0.04, 0.02, 0.01, 0.01] + [0]*6
    omar_naive = np.array(omar_naive)*pro_samp
    omar_treated = np.array(omar_treated)*pro_samp
    omar_1 = np.array(omar_1)*pro_samp
    omar_both = (omar_naive + omar_treated)/2.
    pro_only_1_samples = np.array([12, 12, 13, 11, 6, 11, 11, 7,5, 3, 1, 0])#*163/92.
    pro_2_plus_samples = np.array(p_n_muts_41)-pro_only_1_samples


    #ax1.plot(pro_x, p_n_muts, 'k-s', label='Scripps all pos', lw=2)
    ax1.plot(pro_x_41, p_n_muts_41, 'k-s', label='Scripps All Samples', lw=2)
    #ax1.plot(pro_x_41, pro_only_1_samples, 'k--s', label='Scripps 1 PI Samples', lw=2)
    #ax1.plot(pro_only_1_samples, 'm--s')
    #ax1.plot(pro_2_plus_samples, 'm-s')
    ax1.plot(omar_x, omar_naive, 'b-o', label='Stanford Naive', lw=2)
    ax1.plot(omar_x, omar_1, 'r-o', label='Stanford Treated 1 PI', lw=2)
    #ax1.plot(pro_x_45, p_n_muts_45, 'r--s', label='Scripps 45 pos', lw=2)
    ax1.plot(omar_x, omar_treated, 'g-o', label='Stanford Treated 2+ PI', lw=2)
    #ax1.plot(omar_x, omar_both, 'c--o', label='Omar Both', lw=2)
    #ax2.plot(gag_x, g_n_muts, 'r-D', label='Gag', lw=2)
    #ax2.plot(pro_x*5, p_n_muts*.2, 'k-s', label='Pro projected over Gag', lw=2)
    ax1.legend(loc=1)
    #ax2.legend(loc=5)

    ax1.fill_between(omar_x, omar_naive, 0, color='blue', alpha=0.1)
    ax1.fill_between(omar_x, omar_treated, 0, color='green', alpha=0.1)
    ax1.fill_between(omar_x, omar_1, 0, color='red', alpha=0.1)

    tick = lambda x: ticker.FuncFormatter(lambda y, pos: ('%.1f%%')%(y*100/x))
    ax1.set_xlim(0, max(pro_d.keys()))
    #ax2.set_xlim(0, 75)#max(gag_d.keys()))
    for ax, ylim in zip((ax1,), (pro_samp,)):#zip((ax1, ax2), (pro_samp, gag_samp)):
        ax.set_ylim(0, ylim/4.)
        ax.set_yticks(np.arange(0, ylim/4., ylim/40.))
        ax.yaxis.set_major_formatter(tick(ylim))
        ax.grid()
        ax.set_xlabel('Number of Dominant Mutations')
        ax.set_ylabel('Percentage of %i Samples'%ylim)
    #ax1.set_title('PRO')
    #ax2.set_title('GAG')

    #more10 = [m for m in pro_n_muts if m >= 10]
    ##print sum(more10)
    ##print len(more10)
    ##print len([m for m in gag_n_muts if m >= 50])

    #col_labels = ['PRO', 'GAG']
    #row_labels = ['# of Muts     0 /  0 <  2 / 10',
    #              '# of Muts     2 / 10 <  4 / 20',
    #              '# of Muts     4 / 20 <  6 / 30',
    #              '# of Muts     6 / 30 <  8 / 40',
    #              '# of Muts     8 / 40 < 10 / 50',
    #              '# of Muts    10 / 50 < 15 / 75']
    #row_labels = ['# of Muts     0 <  2,     0 < 10',
    #              '# of Muts     2 <  4,   10 < 20',
    #              '# of Muts     4 <  6,   20 < 30',
    #              '# of Muts     6 <  8,   30 < 40',
    #              '# of Muts     8 < 10,  40 < 50',
    #              '# of Muts   10 < 15,  50 < 75']
    #vals = [(len([m for m in pro_n_muts if i <= m < i+2]),
    #         len([m for m in gag_n_muts if i*5 <= m < (i+2)*5]))
    #         for i in range(0, 12, 2)]
    #vals[-1] = (len([m for m in pro_n_muts if 10 <= m]),
    #            len([m for m in gag_n_muts if 50 <=m]))
    #table = ax2.table(cellText=vals, colWidths = [0.2]*2, rowLabels=row_labels,
    #                  colLabels=col_labels, loc='upper right', cellLoc='center')
    #table.update(dict(alpha=1.))
    #table.set_fontsize(20)

    #params = dict(histtype='stepfilled', alpha=0.8)#, align='left')
    #muts_5 = flatten([[i]*pro_n_muts[i]*5 for i in range(len(pro_n_muts))
    #                  if pro_n_muts[i] > 0])
    #muts_5 = np.array([m for m in muts_5])
    #tmp_muts = flatten([[i]*pro_n_muts[i] for i in range(len(pro_n_muts))
    #                    if pro_n_muts[i] > 0])
    #p_muts = np.array([m for m in tmp_muts])
    #tmp_muts = flatten([[i]*gag_n_muts[i] for i in range(len(gag_n_muts))
    #                    if gag_n_muts[i] > 0])
    #g_muts = np.array([m for m in tmp_muts])

    #ax3.hist(muts_5, bins=pro_samp, range=(0, pro_samp), color='r', **params)
    #ax3.hist(p_muts, bins=pro_samp, range=(0, pro_samp), **params)

    #ax4.hist(muts_5, bins=pro_samp, range=(0, pro_samp), color='r', **params)
    #ax4.hist(g_muts, bins=gag_samp, range=(0, gag_samp), color='g', **params)

    #tick = lambda x: ticker.FuncFormatter(lambda y, pos: ('%.0f%%')%(y*100/x))
    #for ax, samp, ylim in zip((ax3, ax4), (pro_samp, gag_samp), (99, 499)):
    #    ax.set_xlim(0, 1)
    #    ax.set_xticks(np.arange(0, 1, .05))
    #    ax.xaxis.set_major_formatter(tick(1))
    #    ax.set_ylim(0, ylim)
    #    ax.set_yticks(np.arange(0, ylim, ylim/10.))
    #    ax.yaxis.set_major_formatter(tick(ylim))
    #    ax.grid()
    #    ax.set_xlabel('Mutation Rate over %i samples'%samp)
    #    ax.set_ylabel('Percentage of %i Positions'%ylim)
    #ax3.set_title('PRO')
    #ax4.set_title('GAG')

    plt.show()


def run_main(directory, gag=False):
    position_outfile = '/u2/scripps/samples/statistics_per_position.csv'
    patient_outfile = '/u2/scripps/samples/statistics_per_patient.csv'
    if gag:
        position_outfile = position_outfile[:-4] + '_gag.csv'
        patient_outfile = patient_outfile[:-4] + '_gag.csv'
    series = get_series(directory, gag)
    samples = []
    n_dist = defaultdict(int)

    for s in series:
        if len(s) < 2: continue
        series_freqs = []
        allgood = 1
        for sample in s:
            arr, counts = get_mutations_aa(sample, 0.05, gag=gag)
            if np.any(counts < 1e-8):
                allgood = 0
                continue
            fs = arr / counts.reshape(counts.shape[0], 1)
            series_freqs.append(fs)
        if allgood:
            n_dist[len(s)] += 1
            samples.append(series_freqs)
    print n_dist.items()
    #per_position_stats(samples, position_outfile, gag=gag)
    #per_patient_stats(samples, patient_outfile, gag=gag)
    plot_position(position_outfile, n_dist, gag=gag)
    #plot_patient(patient_outfile, n_dist, gag=gag)

def run_secondary():
    pro_pos = '/u2/scripps/samples/statistics_per_position.csv'
    pro_pat = '/u2/scripps/samples/statistics_per_patient.csv'
    gag_pos = '/u2/scripps/samples/statistics_per_position_gag.csv'
    gag_pat = '/u2/scripps/samples/statistics_per_patient_gag.csv'

    #plot_monomorphic([pro_pos, gag_pos, pro_pat, gag_pat])
    #plot_reversion_frac([pro_pos, gag_pos, pro_pat, gag_pat])
    #plot_statics((pro_pos, gag_pos))

def run_auxillary():
    auxillary_stats()

def run_polymorphic():
    pro_poly_muts = [12, 14, 15, 19, 37, 41, 57, 64, 69]
    pro_nonpoly_muts = [10, 11, 13, 20, 22, 23, 24, 30, 32, 33, 34, 35, 36,
                        43, 45, 46, 47, 48, 50, 53, 54, 55, 58, 60, 62, 63,
                        66, 71, 72, 73, 74, 75, 76, 77, 79, 82, 83, 84, 85,
                        88, 89, 90, 92, 93, 95]
    gag_nonpoly_muts = [128, 132, 360, 362, 363, 368, 373, 374, 375, 376,
                        381, 428, 430, 431, 436, 437, 449, 451, 452, 453,
                        12, 62, 75, 76, 79, 81, 112, 200, 219, 369, 370,
                        371, 389, 390, 401, 409, 468, 474, 487, 497]

    PRO_POLY_MUTS = [10, 12, 13, 14, 15, 35, 36, 37, 41, 62, 63, 64, 77,
                    93]#, 18, 19, 57, 67]
    GAG_POLY_MUTS = [362, 370, 371, 373]
    GAG_SIG_MUTS = list(set([128, 132, 360, 362, 363, 368, 373, 374, 375, 376,
                             381, 428, 430, 431, 436, 437, 449, 451, 452, 453,
                             12, 62, 75, 76, 79, 81, 112, 200, 219, 369, 370,
                             371, 389, 390, 401, 409, 468, 474, 487, 497]))
    PRO_SIG_MUTS = [10, 11, 13, 20, 22, 23, 24, 30, 32, 33, 34, 35, 36, 43,
                    45, 46, 47, 48, 50, 53, 54, 55, 58, 60, 62, 63, 66, 71,
                    72, 73, 74, 75, 76, 77, 79, 82, 83, 84, 85, 88, 89, 90,
                    92, 93, 95]
    PRO_POLY_ONLY_MUTS = list(set(PRO_POLY_MUTS) - set(PRO_SIG_MUTS))
    PRO_SIG_ONLY_MUTS = list(set(PRO_SIG_MUTS) - set(PRO_POLY_MUTS))
    PRO_SIG_POLY_MUTS = list(set(PRO_POLY_MUTS) & set(PRO_SIG_MUTS))
    PRO_NON_ASSOC = list(set(range(1,100)) - set(PRO_POLY_MUTS) - set(PRO_SIG_MUTS))
    GAG_POLY_ONLY_MUTS = list(set(GAG_POLY_MUTS) - set(GAG_SIG_MUTS))
    GAG_SIG_ONLY_MUTS = list(set(GAG_SIG_MUTS) - set(GAG_POLY_MUTS))
    GAG_SIG_POLY_MUTS = list(set(GAG_POLY_MUTS) & set(GAG_SIG_MUTS))
    GAG_NON_ASSOC = list(set(range(1, 500)) - set(GAG_SIG_POLY_MUTS) - set(GAG_SIG_ONLY_MUTS))

    pro_file = '/u2/scripps/samples/statistics_per_position.csv'
    gag_file = '/u2/scripps/samples/statistics_per_position_gag.csv'
    pro_data = []
    gag_data = []
    pro_poly = []
    pro_nonpoly = []
    pro_all = []
    gag_poly = []
    gag_nonpoly = []
    gag_all = []
    with open(pro_file, 'r') as f:
        f.next()
        for line in f:
            line = line.strip().split(',')
            pro_data.append(map(int, line[1:4]))
    with open(gag_file, 'r') as f:
        f.next()
        for line in f:
            line = line.strip().split(',')
            gag_data.append(map(int, line[1:4]))
    pro_data = np.array(pro_data)
    gag_data = np.array(gag_data)

    pro_poly_only, pro_sig_only, pro_both = [], [] ,[]
    gag_poly_only, gag_sig_only, gag_both = [], [], []
    pro_non_assoc, gag_non_assoc = [], []

    for pos in range(99):
        f, r, m = pro_data[pos]
        #if not r: continue
        frac = float(r)/f
        tup = (pos+1, r, f, frac, m)
        pro_all.append(tup)
        if pos+1 in pro_poly_muts:
            pro_poly.append(tup)
        if pos+1 in pro_nonpoly_muts:
            pro_nonpoly.append(tup)
        if pos+1 in PRO_POLY_ONLY_MUTS:
            pro_poly_only.append(tup)
        if pos+1 in PRO_SIG_ONLY_MUTS:
            pro_sig_only.append(tup)
        if pos+1 in PRO_SIG_POLY_MUTS:
            pro_both.append(tup)
        if pos+1 in PRO_NON_ASSOC and f > 0:
            pro_non_assoc.append(tup)

    pro_poly_only, pro_sig_only, pro_both = map(np.array,
                [pro_poly_only, pro_sig_only, pro_both])
    pro_poly_only[np.isnan(pro_poly_only)] = 0
    pro_sig_only[np.isnan(pro_sig_only)] = 0
    pro_both[np.isnan(pro_both)] = 0
    pro_non_assoc = np.array(pro_non_assoc)
    pro_non_assoc[np.isnan(pro_non_assoc)] = 0

    np.seterr(all='ignore')
    stat_func = lambda x: x[:, 3].mean()
    print "PRO"
    print 'poly', stat_func(pro_poly_only), len(np.nonzero(pro_poly_only[:,2])[0])
    print sum(pro_poly_only[:, 1]), sum(pro_poly_only[:, 2])
    print 'sig', stat_func(pro_sig_only), len(np.nonzero(pro_sig_only[:,2])[0])
    print 'both', stat_func(pro_both), len(np.nonzero(pro_both[:,2])[0])
    print 'non', stat_func(pro_non_assoc), len(np.nonzero(pro_non_assoc[:,2])[0])
    print sum(pro_non_assoc[:, 1]), sum(pro_non_assoc[:, 2])



    pro_poly = np.array(pro_poly)
    pro_nonpoly = np.array(pro_nonpoly)
    pro_all = np.array(pro_all)
    pro_poly[np.isnan(pro_poly)] = 0
    pro_nonpoly[np.isnan(pro_nonpoly)] = 0
    pro_all[np.isnan(pro_all)] = 0

    #print "pro"
    #print pro_poly[:,3].mean()
    #print pro_nonpoly[:,3].mean()
    print pro_all[:,3].mean()

    gag_non_assoc = []
    for pos in range(499):
        f, r, m = gag_data[pos]
        #if not r: continue
        frac = float(r)/f
        tup = (pos+1, r, f, frac, m)
        gag_all.append(tup)
        if pos+1 in gag_nonpoly_muts:
            gag_nonpoly.append(tup)
        if f > 0 and pos+1 not in gag_nonpoly_muts:
            gag_poly.append(tup)
        if pos+1 in GAG_POLY_ONLY_MUTS:
            gag_poly_only.append(tup)
        elif pos+1 in GAG_SIG_ONLY_MUTS:
            gag_sig_only.append(tup)
        elif pos+1 in GAG_SIG_POLY_MUTS:
            gag_both.append(tup)
        elif pos+1 in GAG_NON_ASSOC and f > 0:
            gag_non_assoc.append(tup)

    gag_poly_only, gag_sig_only, gag_both = map(np.array,
                [gag_poly_only, gag_sig_only, gag_both])
    gag_poly_only[np.isnan(gag_poly_only)] = 0
    gag_sig_only[np.isnan(gag_sig_only)] = 0
    gag_both[np.isnan(gag_both)] = 0
    gag_non_assoc = np.array(gag_non_assoc)

    np.seterr(all='ignore')
    stat_func = lambda x: x[:, 3].mean()
    print "GAG"
    #print 'poly', stat_func(gag_poly_only)
    print 'sig', stat_func(gag_sig_only), len(np.nonzero(gag_sig_only[:,2])[0])
    print 'both', stat_func(gag_both), np.nonzero(gag_both), len(gag_both)
    print 'non', stat_func(gag_non_assoc)
    print sum(gag_non_assoc[:,1]), sum(gag_non_assoc[:, 2]), len(gag_non_assoc)
    print sorted(GAG_SIG_ONLY_MUTS)

    gag_poly = np.array(gag_poly)
    gag_nonpoly = np.array(gag_nonpoly)
    gag_all = np.array(gag_all)
    gag_poly[np.isnan(gag_poly)] = 0
    gag_nonpoly[np.isnan(gag_nonpoly)] = 0
    gag_all[np.isnan(gag_all)] = 0

    #print "gag"
    #print gag_poly[:, 3].mean()
    #print gag_nonpoly[:, 3].mean()
    print 'all', gag_all[:, 3].mean()

    count = 0
    for position in range(499):
        if gag_all[position, 3] > .67 and gag_all[position, 2] > 1:
            count += 1
            print gag_all[position]
    print 'pos that meet criteria', count

    sys.exit()

    print gag_all[:, 1].sum(), gag_all[:, 2].sum()
    print sorted(gag_all.tolist(), key=lambda k: k[2])[-5:]

    print 'SUMMARY STATS'
    print pro_nonpoly[:, 1].sum(), pro_nonpoly[:, 2].sum(), pro_poly[:, 1].sum(), pro_poly[:, 2].sum()
    print gag_nonpoly[:, 1].sum(), gag_nonpoly[:, 2].sum()

    print "MUTS"
    print pro_nonpoly[:, -1].sum(), pro_nonpoly[:, -1].mean()
    print gag_nonpoly[:, -1].sum(), gag_nonpoly[:, -1].mean()


if __name__ == "__main__":
    directory = '/u2/scripps/samples/'
    #run_main(directory, gag=0)
    run_main(directory, gag=1)

    #run_secondary()
    #run_auxillary()
    #plot_sample_muts()
    #run_polymorphic()
