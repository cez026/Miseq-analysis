import sqlite3
import warnings
import numpy as np
from wildtype import AA_TO_IND, IND_TO_AA, WILDTYPE_GAG, WILDTYPE_PRO

warnings.simplefilter('ignore')

DB_LOCATION = '/u2/scripps/samples/hiv.db'
PRO_START = 1463
GAG_END = 1498

def read_database_coverage_new(gag=True):
    conn = sqlite3.connect(DB_LOCATION, detect_types=sqlite3.PARSE_DECLTYPES)
    cursor = conn.cursor()

    pro_coverage = {}
    gag_coverage = {}
    for protein, d in zip(('gag', 'pr'), (gag_coverage, pro_coverage)):
        cursor.execute('SELECT id FROM samples')
        for (sampleID, ) in cursor.fetchall():
            counts = np.zeros(99) if protein == 'pr' else np.zeros(500)
            cursor.execute('SELECT position, sum(count) FROM sequence\
                            WHERE sampleID = ? and protein = ?\
                            GROUP BY position', (sampleID, protein,))
            for position, count in cursor.fetchall():
                counts[position - 1] += count

            d[sampleID] = counts
    return gag_coverage if gag else pro_coverage

def read_database_coverage(gag=True):
    conn = sqlite3.connect(DB_LOCATION, detect_types=sqlite3.PARSE_DECLTYPES)
    cursor = conn.cursor()

    pro_coverage = {}
    gag_coverage = {}
    for protein, d in zip(('gag', 'pr'), (gag_coverage, pro_coverage)):
        cursor.execute('SELECT id FROM samples')
        for (sampleID,) in cursor.fetchall():
            counts = [0,]

            cursor.execute('SELECT position, count FROM sequence\
                            WHERE sampleID = ? and protein = "%s"'%protein, (sampleID,))
            cur_pos = 1
            for position, count in cursor.fetchall():
                if position == cur_pos:
                    counts[position-1] += count
                else:
                    cur_pos += 1
                    counts.append(count)

            d[sampleID] = counts

    #coverage = {}
    #for sample in pro_coverage.keys():
    #    counts = gag_coverage[sample][:PRO_START/3]
    #    counts += pro_coverage[sample]
    #    coverage[sample] = counts

    return gag_coverage if gag else pro_coverage

def read_database_freqs(gag=True, threshold=None, exclude_wt=False):
    protein = 'gag' if gag else 'pr'
    wildtype = WILDTYPE_GAG if gag else WILDTYPE_PRO

    conn = sqlite3.connect(DB_LOCATION, detect_types=sqlite3.PARSE_DECLTYPES)
    cursor = conn.cursor()
    cursor.execute('SELECT id FROM samples')

    freqs = {}
    for (sampleID,) in cursor.fetchall():
        sample_counts = np.zeros((500, 22)) if gag else np.zeros((99, 22))

        cursor.execute('SELECT position, count, wt, residue FROM sequence\
                        WHERE sampleID = ? and protein = "%s"'%protein, (sampleID,))
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

        if exclude_wt:
            s1, s2 = sample_freqs.shape
            for i in range(sample_freqs.shape[0]):
                sample_freqs[i, AA_TO_IND[wildtype[i]]] = 0.
            #sample_freqs = sample_freqs[sample_freqs >= 0].reshape(s1, s2-1)

        freqs[sampleID] = sample_freqs

    return freqs

def read_database_freqs_new(gag=True, threshold=None, exclude_wt=False):
    protein = 'gag' if gag else 'pr'
    wildtype = WILDTYPE_GAG if gag else WILDTYPE_PRO

    conn = sqlite3.connect(DB_LOCATION, detect_types=sqlite3.PARSE_DECLTYPES)
    cursor = conn.cursor()
    cursor.execute('SELECT id FROM samples')

    freqs = {}
    for (sampleID,) in cursor.fetchall():
        sample_freqs = np.zeros((500, 22)) if gag else np.zeros((99, 22))

        cursor.execute('SELECT position, frequency, count, wt, residue FROM sequence\
                        WHERE sampleID = ? and protein = "%s"'%protein, (sampleID,))
        for position, frequency, count, wt_flag, residue in cursor.fetchall():
            aa_ind = AA_TO_IND[residue]
            sample_freqs[position-1, aa_ind] = frequency

        if threshold is not None:
            for i in range(len(sample_freqs)):
                for j in range(len(sample_freqs[i])):
                    if sample_freqs[i,j] < threshold:
                        sample_freqs[i, j] = 0.

        if exclude_wt:
            s1, s2 = sample_freqs.shape
            for i in range(sample_freqs.shape[0]):
                sample_freqs[i, AA_TO_IND[wildtype[i]]] = -1
            sample_freqs = sample_freqs[sample_freqs >= 0].reshape(s1, s2-1)

        freqs[sampleID] = sample_freqs

    return freqs

if __name__ == "__main__":
    freqs = read_database_freqs(gag=0)
    print len(freqs.keys())
    print freqs.values()[0]
