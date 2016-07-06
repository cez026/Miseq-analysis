import re
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

def main4(pairs):
    for pair in pairs:
        assert len(re.findall(r'^[0-9]{1,2}-[0-9]{1,2}$', pair)) > 0

    data = [] #ultimate N_samples x N_pairs X data_size
    for d in glob.glob('/u2/scripps/samples/*'):
        if d.split(os.sep)[-1][0] not in map(str, range(10)): continue

        basepath = d
        if os.path.exists(os.path.join(basepath, 'skip')): continue

        sample_data = []
        for pair in pairs:
            bi = np.loadtxt(os.path.join(basepath, 'pair_bimarg_%s'%pair))
            lower = np.loadtxt(os.path.join(basepath, 'pair_lower_%s'%pair))
            upper = np.loadtxt(os.path.join(basepath, 'pair_upper_%s'%pair))
            assert 0 <= lower[0] <= bi[0] <= upper[0] <= 1 \
                    or np.allclose(lower[0], bi[0]) or np.allclose(bi[0], upper[0])

            sample_data.append([bi[0], bi[0] - lower[0], upper[0] - bi[0]])
        data.append(sample_data)
    data = np.array(data)
    data = data.reshape(len(pairs), -1, 3)
    N_pairs, N_samples, _ = data.shape
    data = data[(np.repeat(range(N_pairs), N_samples).reshape(N_pairs, -1), data[:,:,0].argsort())]

    x = np.arange(N_samples)
    y = data[:,:,0]
    yerr = data[:,:,1:]

    for k in range(0, N_pairs, 2):
        fig, axarr = plt.subplots(2, 1, sharex=True, sharey=True)
        #fig.suptitle('Location of true double mutant probability within estimated bounds', size=20)


        params = {'barsabove': True, 'ecolor': 'k', 'elinewidth': 1.8, 'markeredgewidth': 1., 'alpha':0.5,
                'markersize': 6, 'markeredgecolor': 'r', 'markerfacecolor': 'r'}
        for i, (ax, pair) in enumerate(zip(axarr, pairs)):
            ax.plot([0, N_samples], [0, 0], 'k-', alpha=0.6)
            ax.plot([0, N_samples], [1, 1], 'k-', alpha=0.6)
            ax.errorbar(x, y[i+k], yerr=yerr[i+k].T, fmt='s', **params)

            ax.set_xlim(-0.5 + N_samples/2, N_samples - 0.5)
            ax.set_ylim(-0.1, 1.1)

            ax.set_title(pair, size=20)

        [ax.axis('off') for i, ax in enumerate(axarr.flat) if i >= N_pairs]
        axarr[0].set_xlabel('Sample', size=20)
        [ax.set_ylabel('Double mutant probability', size=20) for ax in axarr]

        fig.subplots_adjust(left=0.04, right=0.98, wspace=0.02, hspace=0.12)

    plt.show()

def main3():
    pairs = ['30-88', '54-82', '73-90', '46-82']
    for pair in pairs:
        assert len(re.findall(r'^[0-9]{1,2}-[0-9]{1,2}$', pair)) > 0

    data = []
    for d in glob.glob('/u2/scripps/samples/*'):
        if d.split(os.sep)[-1][0] not in map(str, range(10)): continue

        basepath = d
        if os.path.exists(os.path.join(basepath, 'skip')): continue

        sample_data = []
        for pair in pairs:
            bi = np.loadtxt(os.path.join(basepath, 'pair_bimarg_%s'%pair))
            lower = np.loadtxt(os.path.join(basepath, 'pair_lower_%s'%pair))
            upper = np.loadtxt(os.path.join(basepath, 'pair_upper_%s'%pair))
            assert 0 <= lower[0] <= bi[0] <= upper[0] <= 1 \
                    or np.allclose(lower[0], bi[0]) or np.allclose(bi[0], upper[0])

            sample_data.append([bi[0], bi[0] - lower[0], upper[0] - bi[0]])
        data.append(sample_data)
    data = np.array(data)
    data = data.reshape(len(pairs), -1, 3)
    N_pairs, N_samples, _ = data.shape
    data = data[(np.repeat(range(N_pairs), N_samples).reshape(N_pairs, -1), data[:,:,0].argsort())]

    x = np.arange(N_samples)
    y = data[:,:,0]
    yerr = data[:,:,1:]
    print x.shape, y.shape, yerr.shape

    fig, ax = plt.subplots()

    params = {'barsabove': True, 'ecolor': 'k', 'elinewidth': 1., 'markeredgewidth': 1., 'alpha':0.4,
            'markersize': 4, 'markeredgecolor': 'r', 'markerfacecolor': 'r'}
    colors = ['r', 'b', 'g', 'c']
    for i, pair in enumerate(pairs):
        params['markerfacecolor'] = colors[i] 
        params['markeredgecolor'] = colors[i]
        ax.plot([0, N_samples], [0, 0], 'k-', alpha=0.6)
        ax.plot([0, N_samples], [1, 1], 'k-', alpha=0.6)
        ax.errorbar(x, y[i], yerr=yerr[i].T, fmt='s', label=pair, **params)

        ax.set_xlim(-0.5 + N_samples/2, N_samples - 0.5)
        ax.set_ylim(-0.1, 1.1)

    ax.set_xlabel('Sample', size=20)
    ax.set_ylabel('Double mutant probability', size=20)

    ax.legend(loc='best', fancybox=True)

    fig.subplots_adjust(left=0.04, right=0.98, wspace=0.02, hspace=0.12)
    plt.show()

def main2(pairs):
    for pair in pairs:
        assert len(re.findall(r'^[0-9]{1,2}-[0-9]{1,2}$', pair)) > 0

    data = [] #ultimate N_samples x N_pairs X data_size
    for d in glob.glob('/u2/scripps/samples/*'):
        if d.split(os.sep)[-1][0] not in map(str, range(10)): continue

        basepath = d
        if os.path.exists(os.path.join(basepath, 'skip')): continue

        sample_data = []
        for pair in pairs:
            bi = np.loadtxt(os.path.join(basepath, 'pair_bimarg_%s'%pair))
            lower = np.loadtxt(os.path.join(basepath, 'pair_lower_%s'%pair))
            upper = np.loadtxt(os.path.join(basepath, 'pair_upper_%s'%pair))
            assert 0 <= lower[0] <= bi[0] <= upper[0] <= 1 \
                    or np.allclose(lower[0], bi[0]) or np.allclose(bi[0], upper[0])

            sample_data.append([bi[0], bi[0] - lower[0], upper[0] - bi[0]])
        data.append(sample_data)
    data = np.array(data)
    data = data.reshape(len(pairs), -1, 3)
    N_pairs, N_samples, _ = data.shape
    data = data[(np.repeat(range(N_pairs), N_samples).reshape(N_pairs, -1), data[:,:,0].argsort())]

    x = np.arange(N_samples)
    y = data[:,:,0]
    yerr = data[:,:,1:]

    for k in range(0, N_pairs, 4):
        fig, axarr = plt.subplots(4, 1, sharex=True, sharey=True)
        fig.suptitle('Location of true double mutant probability within estimated bounds', size=20)

        #axarr = axarr[::-1, :]
        axarr = axarr[:,None]

        params = {'barsabove': True, 'ecolor': 'k', 'elinewidth': 1., 'markeredgewidth': 1., 'alpha':0.4,
                'markersize': 4, 'markeredgecolor': 'r', 'markerfacecolor': 'r'}
        print pairs[k:k+4]
        for i, (ax, pair) in enumerate(zip(axarr.flat, pairs[k:k+4])):
            ax.plot([0, N_samples], [0, 0], 'k-', alpha=0.6)
            ax.plot([0, N_samples], [1, 1], 'k-', alpha=0.6)
            ax.errorbar(x, y[i+k], yerr=yerr[i+k].T, fmt='s', **params)

            ax.set_xlim(-0.5 + N_samples/2, N_samples - 0.5)
            ax.set_ylim(-0.1, 1.1)

            ax.set_title(pair, size=20)

        [ax.axis('off') for i, ax in enumerate(axarr.flat) if i >= N_pairs]
        [ax.set_xlabel('Sample', size=20) for ax in axarr[0,:].flat]
        [ax.set_ylabel('Double mutant probability', size=20) for ax in axarr[:,0].flat]
        plt.setp([ax.get_xticklabels() for ax in axarr[1:,:].flat], visible=False)
        plt.setp([ax.get_yticklabels() for ax in axarr[:,1:].flat], visible=False)

        fig.subplots_adjust(left=0.04, right=0.98, wspace=0.02, hspace=0.12)

    plt.show()


def main5(pairs):
    for pair in pairs:
        assert len(re.findall(r'^[0-9]{1,2}-[0-9]{1,2}$', pair)) > 0

    data = [] #ultimate N_samples x N_pairs X data_size
    for d in glob.glob('/u2/scripps/samples/*'):
        if d.split(os.sep)[-1][0] not in map(str, range(10)): continue

        basepath = d
        if os.path.exists(os.path.join(basepath, 'skip')): continue

        sample_data = []
        for pair in pairs:
            bi = np.loadtxt(os.path.join(basepath, 'pair_bimarg_%s'%pair))
            lower = np.loadtxt(os.path.join(basepath, 'pair_lower_%s'%pair))
            upper = np.loadtxt(os.path.join(basepath, 'pair_upper_%s'%pair))
            assert 0 <= lower[0] <= bi[0] <= upper[0] <= 1 \
                    or np.allclose(lower[0], bi[0]) or np.allclose(bi[0], upper[0])

            sample_data.append([bi[0], bi[0] - lower[0], upper[0] - bi[0]])
        data.append(sample_data)
    data = np.array(data)
    data = data.reshape(len(pairs), -1, 3)
    N_pairs, N_samples, _ = data.shape
    data = data[(np.repeat(range(N_pairs), N_samples).reshape(N_pairs, -1), data[:,:,0].argsort())]

    x = np.arange(N_samples)
    y = data[:,:,0]
    yerr = data[:,:,1:]

    for k in range(0, N_pairs, 4):
        fig, axarr = plt.subplots(4, 1, sharex=True, sharey=True)
        #fig.suptitle('Location of true double mutant probability within estimated bounds', size=20)

        #axarr = axarr[::-1, :]

        params = {'barsabove': True, 'ecolor': 'k', 'elinewidth': 1., 'markeredgewidth': 1., 'alpha':0.4,
                'markersize': 4, 'markeredgecolor': 'r', 'markerfacecolor': 'r'}
        print pairs[k:k+4]
        for i, (ax, pair) in enumerate(zip(axarr, pairs[k:k+4])):
            ax.plot([0, N_samples], [0, 0], 'k-', alpha=0.6)
            ax.plot([0, N_samples], [1, 1], 'k-', alpha=0.6)
            ax.errorbar(x, y[i+k], yerr=yerr[i+k].T, fmt='s', **params)

            ax.set_xlim(-0.5 + N_samples/2, N_samples - 0.5)
            ax.set_ylim(-0.1, 1.1)

            ax.set_title(pair, size=20)

        [ax.axis('off') for i, ax in enumerate(axarr) if i >= N_pairs]
        fig.text(0.5, 0.04, 'Sample', size=20, ha='center', va='center')
        fig.text(0.06, 0.5, 'Double mutant probability', size=20, ha='center', va='center', rotation='vertical')
        #[ax.set_xlabel('Sample', size=20) for ax in axarr[0,:].flat]
        #[ax.set_ylabel('Double mutant probability', size=20) for ax in axarr[:,0].flat]
        plt.setp([ax.get_xticklabels() for ax in axarr[:-1]], visible=False)
        #plt.setp([ax.get_yticklabels() for ax in axarr[:]], visible=False)
        plt.setp([ax.get_xticklabels() for ax in axarr], fontsize=16)
        plt.setp([ax.get_yticklabels() for ax in axarr], fontsize=16)

        fig.subplots_adjust(left=0.08, right=0.98, wspace=0.02, hspace=0.14)

    plt.show()


def main(pair):
    assert len(re.findall(r'^[0-9]{1,2}-[0-9]{1,2}$', pair)) > 0
    
    data = []
    for bi_file in glob.glob('/u2/scripps/samples/*/pair_bimarg_%s'%pair):
        basepath = os.path.dirname(bi_file)
        if os.path.exists(os.path.join(basepath, 'skip')): continue
        lower_file = os.path.join(basepath, 'pair_lower_%s'%pair)
        upper_file = os.path.join(basepath, 'pair_upper_%s'%pair)

        bi = np.loadtxt(bi_file)
        lower = np.loadtxt(lower_file)
        upper = np.loadtxt(upper_file)

        sample_data = [bi[0], bi[0] - lower[0], upper[0] - bi[0]]
        assert 0 <= lower[0] <= bi[0] <= upper[0] <= 1

        data.append(sample_data)

    data = np.array(data)
    data = data[data[:,0].argsort()]

    x = np.arange(data.shape[0])
    y = data[:,0]
    yerr = data[:,1:].T

    plt.plot([0, data.shape[0]], [0, 0], 'k-')
    plt.plot([0, data.shape[0]], [1, 1], 'k-', alpha=0.4)
    params = {'barsabove': True, 'ecolor': 'r', 'elinewidth': 1.5,
              'markeredgewidth': 1.5}
    plt.errorbar(x, y, yerr=yerr, fmt='o', alpha=0.6, ** params)

    plt.xlim(-0.5, data.shape[0] - 0.5)
    plt.ylim(-0.1, 1.1)

    plt.xlabel('Sample', size=20)
    plt.ylabel('Double mutant probability', size=20)
    plt.title('Location of true %s double mutant probability within estimated bounds'%pair, size=20)
    plt.show()

#if __name__ == "__main__":
#    if len(sys.argv) != 2:
#        print "Usage:  python %s [pair_name]"%argv[0]
#        sys.exit(1)
#
#    main(sys.argv[1])
MUT_PAIRS = ['30-88', '54-82', '73-90', '46-82', '24-74', '35-36', '69-84',
             '24-46', '24-82', '13-33', '10-93', '12-19', '33-66', '10-46',
             '32-82', '24-64', '37-63', '33-60', '41-93', '30-35', '35-88',
             '32-46', '20-62', '63-93']
if __name__ == "__main__":
    main5(MUT_PAIRS)
    #main3()
    #main4(['30-88','46-82'])
