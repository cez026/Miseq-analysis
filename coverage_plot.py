import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from utils.mut_freqs import read_database_coverage_new

def avg_coverage(coverage_dict):
    avg = np.zeros((len(coverage_dict), len(coverage_dict.values()[0])))
    for i, sample_coverage in enumerate(coverage_dict.itervalues()):
        avg[i] = sample_coverage
    avg = np.array(avg)
    print avg.mean()

def gag_coverage():
    coverage = read_database_coverage_new(gag=True)
    avg_coverage(coverage)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for sample in coverage.itervalues():
        ax.plot(sample)

    locs = [(132, 'MA'), (363, 'CA'), (377, 'p2'), (432, 'NC'), (448, 'p1'), (500, 'p6')]
    x_start = 0
    colors = ('#AAAAAA', '#EEEEEE')
    for i, (junc, name) in enumerate(locs):
        color = colors[i%2]
        width = junc - x_start
        rect = Rectangle((x_start, 1700000), width, 100000, edgecolor='None', color=color)
        ax.add_patch(rect)
        ax.text(x_start + width/2., 1750000, name, ha='center', va='center')
        x_start = junc
    #ax.set_yscale('log')

    plt.show()

def pro_coverage():
    coverage = read_database_coverage(gag=False)
    avg_coverage(coverage)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for sample in coverage:
        ax.plot(coverage[sample])
    #ax.set_yscale('log')

    plt.show()

def both_coverage():
    gag = read_database_coverage_new(gag=True)
    pr = read_database_coverage_new(gag=False)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    data = []
    full_data = []
    for g, p in zip(gag.itervalues(), pr.itervalues()):
        full_data.append(np.hstack((g,p)))
        data.append(np.hstack((g,p)).mean())
        ax.plot(np.hstack((g, p)), alpha=0.5)
    data = np.array(data)
    full_data = np.array(full_data)
    print full_data[np.where(data < 300)][full_data[np.where(data < 300)].nonzero()].mean()
    print sorted(data)[:10]
    print sorted(data)[-10:]
    print 'mean:', data.mean()
    print 'std:', data.std()
    print 'min:', data.min()
    print '25%:', np.percentile(data, 25)
    print '50%:', np.percentile(data, 50)
    print '75%:', np.percentile(data, 75)

    locs = [(132, 'MA'), (363, 'CA'), (377, 'p2'), (432, 'NC'), 
            (448, 'p1'), (500, 'p6'), (599, 'PR')]
    x_start = 0
    colors = ('#AAAAAA', '#EEEEEE')
    for i, (junc, name) in enumerate(locs):
        color = colors[i%2]
        width = junc - x_start
        rect = Rectangle((x_start, 1700000), width, 100000, edgecolor='None', color=color)
        ax.add_patch(rect)
        ax.text(x_start + width/2., 1750000, name, ha='center', va='center')
        ax.axvline(x_start, ls='--', color='k', lw=0.5)
        x_start = junc

    ax.set_xlabel('Sequence Position')
    ax.set_ylabel('Coverage Depth')
    mpl.rc('font', **{'size': 22})

    plt.show()

if __name__ == "__main__":
    both_coverage()
    #pro_coverage()
