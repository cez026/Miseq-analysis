import os, sys
import csv
import numpy as np
from itertools import izip_longest
from collections import defaultdict
from read_mutation_files import get_mutations_aa
from utils.wildtype import WILDTYPE_PRO, WILDTYPE_GAG, AA_TO_IND, IND_TO_AA

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
        sample_filenames = []
        for sample in sorted(samples[key]):
            for filename in os.listdir(os.path.join(directory, parent%sample)):
                if filename.endswith(suffix):
                    full = os.path.join(directory, parent%sample, filename)
                    sample_filenames.append(full)
        series.append(sample_filenames)

    return series

def transitions_per_position(series, outfile, gag=False):
    L = 99
    wildtype = WILDTYPE_PRO
    if gag:
        L = 499
        wildtype = WILDTYPE_GAG

    freqs = []
    dist_of_samples = defaultdict(int)
    for s in series:
        if len(s) < 2: continue
        series_freqs = []
        all_good = 1
        for sample in s:
            counts, total = get_mutations_aa(sample, 0.05, gag=gag)
            if np.any(total < 1e-8):
                all_good = 0
                continue
            counts /= total.reshape(total.shape[0], 1)
            series_freqs.append(counts)
        if all_good:
            dist_of_samples[len(s)] += 1
            freqs.append(series_freqs)

    print 'Distribution of Samples', dist_of_samples.items()
    print 'Number of Patients', sum(dist_of_samples.values())
    print 'Total Samples', sum([k*v for k,v in dist_of_samples.items()])

    def print_freqs(freqs, i, j, xi, xj):
        freqs_i = freqs[:, i-1, AA_TO_IND[xi]]
        freqs_j = freqs[:, j-1, AA_TO_IND[xj]]
        if np.all(freqs_i < 1e-8) and np.all(freqs_j < 1e-8):
            pass
        else:
            print freqs_i
            print freqs_j
            print

    for s in freqs:
        series_freqs = np.array(s)
        #print_freqs(series_freqs, 30, 88, 'N', 'D')
        #print_freqs(series_freqs, 54, 82, 'V', 'A')
        #print_freqs(series_freqs, 30, 82, 'N', 'A')
        print_freqs(series_freqs, 71, 90, 'V', 'M')
        #print_freqs(series_freqs, 10, 54, 'I', 'V')

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


if __name__ == "__main__":
    gag = 0
    series = get_series('/u2/scripps/samples', gag)
    transitions_per_position(series, None, gag)
