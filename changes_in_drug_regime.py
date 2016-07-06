import os, sys, time, re
import numpy as np
import sqlite3
from collections import defaultdict
from utils.mut_freqs import read_database_freqs, DB_LOCATION
from utils.wildtype import WILDTYPE_GAG

AA_TO_IND = dict(zip('ACDEFGHIKLMNPQRSTVWY*X', range(22)))
IND_TO_AA = dict(zip(range(22), 'ACDEFGHIKLMNPQRSTVWY*X'))

def get_time_series(threshold=0.01):
    conn = sqlite3.connect(DB_LOCATION, detect_types=sqlite3.PARSE_DECLTYPES)
    cursor = conn.cursor()

    cursor.execute('SELECT id, patientID, date, timepoint, vload from samples')

    patient_series = defaultdict(list)
    for (sample, patient, d, t, vload,) in cursor.fetchall():
        patient_series[patient].append( (t, sample, d, vload) )

    cursor.execute('SELECT id, classM, regimen1, regimen2 from patients')
    regimens = defaultdict(list)
    for (patient, success, reg1, reg2,) in cursor.fetchall():
        regimens[patient].append(success)
        f = lambda x: x.rstrip('/r')
        reg1s = tuple(map(f, reg1.split(':', 1)))
        reg2s = tuple(map(f, reg2.split(':', 1)))
        regimens[patient].append(reg1s)
        regimens[patient].append(reg2s)
        #for reg in (reg1, reg2):
        #    [regimens[patient].append(r) for r in re.findall(r"[\w']+", reg)
        #                                 if len(r) > 2]

    freqs = {}
    cursor.execute('SELECT id from samples')
    for (sampleID,) in cursor.fetchall():
        sample_counts = np.zeros((500, 22))
        cursor.execute('SELECT position, count, wt, residue FROM sequence\
                        WHERE sampleID = ? and protein = "gag"', (sampleID,))

        for position, count, wt_flag, residue in cursor.fetchall():
            aa_ind = AA_TO_IND[residue]
            sample_counts[position-1, aa_ind] = count

        if not np.any(sample_counts > 0): continue

        position_counts = sample_counts.sum(axis=1)
        if threshold is not None:
            for i in range(len(sample_counts)):
                sub = 0
                for j in range(len(sample_counts[i])):
                    freq = sample_counts[i, j] / position_counts[i]
                    if freq < threshold:
                        sub += sample_counts[i, j]
                        sample_counts[i, j] = 0.
                position_counts[i] -= sub
                position_counts[i] = float(position_counts[i])
                if position_counts[i] < 1:
                    position_counts[i] = 1.

        sample_freqs = sample_counts /\
                        position_counts.reshape(position_counts.shape[0], 1)

        freqs[sampleID] = sample_freqs

    full_dict = defaultdict(list)

    for patient in patient_series.keys():
        regimen_info = regimens[patient]
        success, regimen_data = regimen_info[0], regimen_info[1:]

        patient_samples = patient_series[patient]
        patient_freqs = []
        for (t, sample, d, vload) in sorted(patient_samples):
            patient_freqs.append(freqs[sample])
        patient_freqs = np.array(patient_freqs).reshape(len(patient_freqs), 500, 22)

        full_dict[patient] += regimen_data
        full_dict[patient].append(patient_freqs)

    # sort full_dict by regimen
    # regimen is the second one that is switched to
    indinavir = {}
    saquinavir = {}
    nelfinavir = {}
    for pat, vals in full_dict.iteritems():
        if 'IDV' in vals[1]:
            indinavir[pat] = vals
        elif 'SQV' in vals[1]:
            saquinavir[pat] = vals
        elif 'NFV' in vals[1]:
            nelfinavir[pat] = vals

    return indinavir, saquinavir, nelfinavir

def find_resistance_mutations():
    idv, sqv, nfv = get_time_series()
    def find_diff(d):
        muts = defaultdict(lambda: [0, []])
        for patient in d:
            freqs = d[patient][-1]

            first, last = freqs[0,:,:], freqs[-1,:,:]
            for i in range(500):
                diff = last[i] - first[i]
                for xi in range(20):
                    if IND_TO_AA[xi] == WILDTYPE_GAG[i]: continue
                    if abs(diff[xi]) > 0.3:
                        s = str(i+1)+'-'+IND_TO_AA[xi]
                        muts[s][0] += 1
                        muts[s][1] += [patient, '%.2f'%last[i][xi],
                                                '%.2f'%first[i][xi]]
        return muts

    idv_muts = find_diff(idv)
    sqv_muts = find_diff(sqv)
    nfv_muts = find_diff(nfv)

    k = lambda x: int(x[0].split('-', 1)[0])
    print 'IDV\n'
    for mut, score in sorted(idv_muts.items(), key=k):
        if score[0] > 2:
            print mut, score
    print 'SQV\n'
    for mut, score in sorted(sqv_muts.items(), key=k):
        if score[0] > 2:
            print mut, score
    print 'NFV\n'
    for mut, score in sorted(nfv_muts.items(), key=k):
        if score[0] > 2:
            print mut, score


if __name__ == "__main__":
    find_resistance_mutations()
