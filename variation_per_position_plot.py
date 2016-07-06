import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from utils.mut_freqs import read_database_freqs

#cdict = matplotlib.cm.get_cmap('spectral_r')._segmentdata
#cdict = matplotlib.cm.get_cmap('gray_r')._segmentdata
#cdict = dict((k, list(v)) for k,v in cdict.iteritems())
#print cdict['red']
#cdict['red'][0] = (0, 1, 1)
#cdict['green'][0] = (0, 1, 1)
#cdict['blue'][0] = (0, 1, 1)
#my_cmap = LinearSegmentedColormap('bill', cdict)

ALPHABET = 'ACDEFGHIKLMNPQRSTVWY*X'
AA_TO_IND = dict(zip('ACDEFGHIKLMNPQRSTVWY*X', range(22)))
GAG_RESISTANCE_MUTS = ['E12K', 'G62R', 'L75R', 'R76K', 'Y79F', 'T81A', 'K112E',
                  'V128I', 'V128T', 'V128A', 'Y132F', 'M200I', 'H219Q', 'H219Q',
                  'A360V', 'V362I', 'L363M', 'L363F', 'L363C', 'L363N', 'L363Y',
                  'S368C', 'S368N', 'Q369H', 'V370A', 'V370M', 'S373P', 'S373Q',
                  'S373T', 'A374P', 'A374S', 'T375N', 'T375S', 'I376V', 'G381S',
                  'I389T', 'V390A', 'V390D', 'I401T', 'I401V', 'R409K', 'E428G',
                  'Q430R', 'A431V', 'K436E', 'K436R', 'I437T', 'I437V', 'L449F',
                  'L449P', 'L449V', 'S451T', 'S451G', 'S451R', 'R452S', 'R452K',
                  'P453A', 'P453L', 'P453T', 'E468K', 'Q474L', 'A487S', 'P497L']

GAG_RESISTANCE_MUTS_SIMP = [m[1:] for m in GAG_RESISTANCE_MUTS]

GAG_RESISTANCE = [int(m[1:-1]) for m in GAG_RESISTANCE_MUTS]

def get_freqs(gag=False, remove_rm=False):
    freqs = read_database_freqs(gag=gag, threshold=0.01, exclude_wt=True)

    freqs = np.array([freqs[sample] for sample in freqs])
    if remove_rm:
        rms = [(int(m[1:-1]), AA_TO_IND[m[-1]] - int(m[0] < m[-1])) for m in GAG_RESISTANCE_MUTS]
        N, L, Q = freqs.shape
        for i in range(L):
            for j in range(Q):
                if (i+1, j) not in rms:
                    freqs[:, i, j] = 0.
    #freqs = freqs.sum(axis=2)
    #freqs.sort(axis=0)
    #freqs = freqs[::-1]

    return freqs

def plot_cs_variants(save=False):
    matplotlib.rc('font', **{'size': 16})
    locs = [(132, 'MA'), (363, 'CA'), (377, 'p2'),
            (432, 'NC'), (448, 'p1'), (500, 'p6')]

    freqs = get_freqs(gag=True)
    N_samples, L = freqs.shape

    fig, axarr = plt.subplots(2, 3)

    quarter = N_samples / 4.
    cs_locs = [(128, 138, 'MA/CA'), (359, 369, 'CA/p2'), (373, 383, 'p2/NC'),
               (428, 438, 'NC/p1'), (444, 454, 'p1/p6')]
    for i, (cs, ax) in enumerate(zip(cs_locs, axarr.flatten())):
        s, e, name = cs
        ax.pcolor(freqs[:,s-1:e-1], cmap=my_cmap)
        ax.set_title(name)
        ax.set_xticks(np.arange(0.5, e-s, 1))
        ax.set_xticklabels(range(s, e, 1), ha='center')
        ax.set_yticks(np.arange(0, N_samples + quarter, quarter))
        ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax.set_ylim(0, 0.8*N_samples)
        ax.axvline(5, ls='--', color='k', lw=0.5)
        plt.tick_params(axis='x', which='both', bottom='off', top='off')
    axarr[0,0].set_ylabel('Samples with variants')
    axarr[1,1].set_xlabel('Sequence Position')
    axarr[1,2].axis('off')
    plt.show()

def plot_gag_variants_no_rm(save=False):
    from matplotlib.gridspec import GridSpec

    freqs = get_freqs(gag=True)
    N_samples, L = freqs.shape

    #freqs_rm = get_freqs(gag=True, remove_rm=True)
    freqs_rm = freqs.copy()
    for i in range(L):
        if i in GAG_RESISTANCE: continue
        freqs_rm[:, i-1] = 0.
    N_samples, L = freqs.shape

    fig = plt.figure(figsize=(24, 12))
    fig.suptitle('Variant Frequencies in Gag', size=20)
    gs = GridSpec(3, 2, height_ratios=[10,1,10], width_ratios=[100, 1])
    ax1 = plt.subplot(gs[0, :-1])
    ax3 = plt.subplot(gs[1, :-1])
    ax2 = plt.subplot(gs[2, :-1])
    matplotlib.rc('font', size=24)
    matplotlib.rc('axes', linewidth=5)
    [ax.set_xlabel('Sequence Position') for ax in (ax1, ax2)]
    ax1.set_ylabel('Samples with variants')
    ax2.set_ylabel('Samples with variants\nat PI-associated positions')

    plt1 = ax1.pcolor(freqs, cmap=my_cmap)
    plt2 = ax2.pcolor(freqs_rm, cmap=my_cmap)

    plt.setp(ax1.get_xticklabels(), visible=False)
    quarter = N_samples / 4.
    [ax.set_yticks(np.arange(0, N_samples + quarter, quarter)) for ax in (ax1, ax2)]
    [ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%']) for ax in (ax1, ax2)]
    [ax.set_ylim(0, 0.8*N_samples) for ax in (ax1, ax2)]

    locs = [(132, 'MA'), (363, 'CA'), (377, 'p2'), (432, 'NC'), (448, 'p1'), (500, 'p6')]
    x_start = 0
    colors = ('#AAAAAA', '#EEEEEE')
    for i, (junc, name) in enumerate(locs):
        color = colors[i%2]
        width = junc - x_start
        rect = Rectangle((x_start, 0), width, 1, edgecolor='None', color=color)
        ax3.add_patch(rect)
        ax3.text(x_start + width/2., 1/2., name, ha='center', va='center')
        ax1.axvline(x_start, color='k', ls='--', lw=1.0)
        ax2.axvline(x_start, color='k', ls='--', lw=1.0)
        x_start = junc
    ax3.set_xlim(0, L)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    [ax.tick_params(top=False, left=False, right=False, bottom=False) for ax in (ax1, ax2, ax3)]

    # add region brackets
    #ma_l, ma_h = 105, 135
    #nc_l, nc_h = 375, 400
    #p6_l, p6_h = 455, 490
    #[ax.axvspan(ma_l, ma_h, ec='k', color='y', alpha=0.2, lw=0) for ax in (ax1, ax2, ax3)]
    #[ax.axvspan(nc_l, nc_h, ec='k', color='y', alpha=0.2, lw=0) for ax in (ax1, ax2, ax3)]
    #[ax.axvspan(p6_l, p6_h, ec='k', color='y', alpha=0.2, lw=0) for ax in (ax1, ax2, ax3)]
    #ax1.text((ma_l+5), N_samples*0.7, 'A', ha='center', va='center', fontsize=20)
    #ax1.text((nc_l+5), N_samples*0.7, 'B', ha='center', va='center', fontsize=20)
    #ax1.text((p6_l+5), N_samples*0.7, 'C', ha='center', va='center', fontsize=20)

    ax2.invert_yaxis()

    #cax = plt.subplot(gs[0, -1])
    #cbar = fig.colorbar(plt1, cax=cax, ticks=np.arange(0, 1.1, 0.1))
    #cbar.outline.set_linewidth(0.5)
    #cbar.set_label("Variant frequency")
    #cbar.ax.set_yticklabels(['%i%%'%i for i in range(0, 110, 10)])

    fig.subplots_adjust(hspace=0., wspace=0.02, left=0.08, right=0.95, top=0.95, bottom=0.06)
    if not save:
        plt.show()
    else:
        pdf = PdfPages('/home/wflynn/Documents/publications/gag-pol_deep_sequencing/number_of_variants_vs_position_gag_colorized.pdf')
        pdf.savefig(fig)
        pdf.close()


def plot_gag_variants(remove_rm=False, save=False):
    import pylab
    freqs = get_freqs(gag=True)
    if remove_rm:
        for rm in GAG_RESISTANCE:
            freqs[:,rm-1] = 0.
    N_samples, L = freqs.shape

    fig =  pylab.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    matplotlib.rc('font', size=14)
    matplotlib.rc('axes', linewidth=5)
    ax.set_xlabel('Sequence Position')
    ax.set_ylabel('Samples with variants')

    pylab.pcolor(freqs, cmap='pink_r')
    pylab.yticks([0, 0.25, 0.5, 0.75, 1.0], ['0%', '25%', '50%', '75%', '100%'], transform=ax.get_xaxis_transform())

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    [t.set_visible(False) for t in ax.get_xticklines()]
    [t.set_visible(False) for t in ax.get_yticklines()]

    #add gag bar on top
    locs = [(132, 'MA'), (363, 'CA'), (377, 'p2'), (432, 'NC'), (448, 'p1'), (500, 'p6')]
    height = 0.1
    y_start = 0.9
    x_start = 0
    colors = ('.8', '.9')
    for i, (junc, name) in enumerate(locs):
        color = colors[i%2]
        width = junc - x_start
        rect = Rectangle((x_start, y_start), width, height, edgecolor='None', color=color, transform=ax.get_xaxis_transform())
        ax.add_patch(rect)
        ax.text(x_start + width/2., y_start + height/2., name, ha='center', va='center', transform=ax.get_xaxis_transform())
        x_start = junc

    if save:
        pylab.savefig('/home/wflynn/Documents/publications/gag-pol_deep_sequencing/number_of_variants_vs_position_gag.png', dpi=300)
    else:
        pylab.show()

def plot_gag_insets(save=False):

    def inset(freqs, axes, lo, hi, no_y_labels=False):
        N_samples, L = freqs.shape
        freqs_rm = freqs.copy()
        for i in range(L):
            if i in GAG_RESISTANCE: continue
            freqs_rm[:, i-1] = 0.

        #plt1 = axes[0].pcolor(freqs[:, lo:hi], cmap=my_cmap)#, vmin=0.01)
        #plt2 = axes[1].pcolor(freqs_rm[:, lo:hi], cmap=my_cmap)#, vmin=0.01)
        #my_cmap = matplotlib.cm.get_cmap('jet_r')
        #my_cmap.set_under('w')
        plt1 = axes[0].pcolor(freqs[:, lo:hi], cmap=my_cmap)#, vmin=0.01)
        plt2 = axes[1].pcolor(freqs_rm[:, lo:hi], cmap=my_cmap)#, vmin=0.01)

        [ax.set_xlabel('Sequence Position') for ax in axes[:-1]]
        if not no_y_labels:
            axes[0].set_ylabel('Samples with variants')
            axes[1].set_ylabel('Samples with variants at PI-associated positions')

        plt.setp([ax.get_xticklabels() for ax in (axes[0], axes[-1])], visible=False)
        axes[1].set_xticks(np.arange(0, hi-lo+1, 5))
        axes[1].set_xticklabels(np.arange(lo, hi+1, 5))

        quarter = N_samples / 4.
        [ax.set_yticks(np.arange(0, N_samples + quarter, quarter)) for ax in axes[:-1]]
        [ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%']) for ax in axes[:-1]]
        [ax.set_ylim(0, 0.8*N_samples) for ax in axes[:-1]]
        if no_y_labels:
            plt.setp([ax.get_yticklabels() for ax in axes], visible=False)

        locs = [(132, 'MA'), (363, 'CA'), (377, 'p2'), (432, 'NC'), (448, 'p1'), (500, 'p6')]
        x_start = lo
        colors = ('#AAAAAA', '#EEEEEE')
        for i, (junc, name) in enumerate(locs):
            if junc < lo: continue
            color = colors[i%2]
            width = min(junc, hi) - x_start
            rect = Rectangle((x_start, 0), width, 1, edgecolor='None', color=color)
            axes[-1].add_patch(rect)
            axes[-1].text(x_start + width/2., 1/2., name, ha='center', va='center')
            x_start = junc
        axes[-1].set_xlim(lo, hi)
        plt.setp(axes[-1].get_xticklabels(), visible=False)
        plt.setp(axes[-1].get_yticklabels(), visible=False)

        [ax.tick_params(top=False, left=False, right=False, bottom=False) for ax in axes]
        #axes[0].grid(); axes[1].grid()
        return (plt1, plt2)

    from matplotlib.gridspec import GridSpec

    freqs = get_freqs(gag=True)
    N_samples, L = freqs.shape

    fig = plt.figure()#figsize=(12,8))
    fig.suptitle('Regions of dense variation in Gag proteins', size=20)
    gs = GridSpec(3, 4, height_ratios=[10,1,10], width_ratios=[30,30,30,1])
    # Left
    ma_l, ma_h = 105, 135
    ax1 = plt.subplot(gs[0,0])
    ax3 = plt.subplot(gs[1, 0])
    ax2 = plt.subplot(gs[2:,0])
    ax1.text(1, N_samples*0.75, 'A', fontsize=20)
    # Center
    nc_l, nc_h = 375, 400
    ax4 = plt.subplot(gs[0,1])
    ax6 = plt.subplot(gs[1, 1])
    ax5 = plt.subplot(gs[2:,1])
    ax4.text(1, N_samples*0.75, 'B', fontsize=20)
    # Center
    p6_l, p6_h = 455, 490
    ax7 = plt.subplot(gs[0,2])
    ax9 = plt.subplot(gs[1, 2])
    ax8 = plt.subplot(gs[2:,2])
    ax7.text(1, N_samples*0.75, 'C', fontsize=20)
    # Colorbar
    cax = plt.subplot(gs[0, 3])

    plts1 = inset(freqs, [ax1, ax2, ax3], ma_l, ma_h)
    plts2 = inset(freqs, [ax4, ax5, ax6], nc_l, nc_h, no_y_labels=1)
    plts3 = inset(freqs, [ax7, ax8, ax9], p6_l, p6_h, no_y_labels=1)
    [ax.invert_yaxis() for ax in (ax2, ax5, ax8)]
    matplotlib.rc('font', size=20)

    cbar = fig.colorbar(plts3[0], cax=cax, ticks=np.arange(0, 1.1, 0.1))
    cbar.set_label("Variant frequency")
    cbar.ax.set_yticklabels(['%i%%'%i for i in range(0, 110, 10)])

    fig.subplots_adjust(hspace=0., wspace=0.05, left=0.05, right=0.95, top=0.95, bottom=0.05)
    if not save:
        plt.show()
    else:
        fig.savefig('/home/wflynn/Documents/publications/gag-pol_deep_sequencing/regions_of_dense_variation_gag_colorized.pdf', dpi=300, orientation='landscape', format='pdf')

def plot_pro_variants(save=False):
    import pylab
    freqs = get_freqs(gag=False)
    N_samples, L = freqs.shape

    fig =  pylab.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    matplotlib.rc('font', size=14)
    matplotlib.rc('axes', linewidth=5)
    ax.set_xlabel('Sequence Position')
    ax.set_ylabel('Samples with variants')

    im=pylab.pcolor(freqs, cmap=my_cmap)#'pink_r')
    pylab.xticks(np.arange(0, 100, 10))
    pylab.yticks([0, 0.25, 0.5, 0.75, 1.0], ['0%', '25%', '50%', '75%', '100%'], transform=ax.get_xaxis_transform())

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    [t.set_visible(False) for t in ax.get_xticklines()]
    [t.set_visible(False) for t in ax.get_yticklines()]

    cax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
    fig.subplots_adjust(hspace=0., wspace=0.05, left=0.05, right=0.90, top=0.95, bottom=0.05)
    cbar = pylab.colorbar(im, cax=cax, ticks=np.arange(0, 1.1, 0.1))
    cbar.set_label("Variant frequency")
    cbar.ax.set_yticklabels(['%i%%'%i for i in range(0, 110, 10)])

    if save:
        pylab.savefig('/home/wflynn/Documents/publications/gag-pol_deep_sequencing/number_of_variants_vs_position_pro.png', dpi=300)
    else:
        pylab.show()

def publication_plot_old(save=False):
    from matplotlib.gridspec import GridSpec

    freqs = get_freqs(gag=True)
    N_samples, L, _ = freqs.shape

    #freqs_rm = get_freqs(gag=True, remove_rm=True)
    freqs_rm = freqs.copy()
    #for i in range(L):
    #    if i in GAG_RESISTANCE: continue
    #    freqs_rm[:, i-1] = 0.
    for i in range(L):
        for xi, ch in enumerate(ALPHABET[:20]):
            if str(i+1) + ch in GAG_RESISTANCE_MUTS_SIMP: continue
            freqs_rm[:, i, xi] = 0.
    freqs = freqs.sum(axis=2)
    freqs.sort(axis=0)
    freqs = freqs[::-1]
    freqs_rm = freqs_rm.sum(axis=2)
    freqs_rm.sort(axis=0)
    freqs_rm = freqs_rm[::-1]


    #lanl_freqs = np.genfromtxt('lanl/gag_freqs_nowt_experienced.txt')
    lanl_freqs = np.genfromtxt('lanl/gag_freqs_nowt_naive.txt')
    lanl_freqs = lanl_freqs.reshape(-1,500,20)
    lanl_rm = lanl_freqs.copy()
    for i in range(L):
        for xi, ch in enumerate(ALPHABET[:20]):
            if str(i+1) + ch in GAG_RESISTANCE_MUTS_SIMP: continue
            lanl_rm[:, i, xi] = 0.
    lanl_freqs = lanl_freqs.sum(axis=2)
    lanl_rm = lanl_rm.sum(axis=2)
    print lanl_freqs.shape, lanl_rm.shape
    print freqs.shape


    #fig = plt.figure(figsize=(24, 12))
    fig = plt.figure(figsize=(7.5, 4))
    fig.suptitle('Variant Frequencies in Gag', size=20)
    gs = GridSpec(3, 2, height_ratios=[10,1,10], width_ratios=[100, 1])
    ax1 = plt.subplot(gs[0, :-1])
    ax3 = plt.subplot(gs[1, :-1])
    ax2 = plt.subplot(gs[2, :-1])
    matplotlib.rc('font', size=24)
    matplotlib.rc('axes', linewidth=5)
    [ax.set_xlabel('Sequence Position') for ax in (ax1, ax2)]
    ax1.set_ylabel('Samples with variants')
    ax2.set_ylabel('LANL PI-naive\nsequences with variants')


    red_cmap = ListedColormap(['red'])
    black_cmap = ListedColormap(['black'])
    red_cmap.set_bad((0,0,0,0))
    freqs_ma = np.ma.masked_where(freqs < 0.01, freqs)
    freqs_rm_ma = np.ma.masked_where(freqs_rm < 0.01, freqs_rm)
    plt1 = ax1.pcolor(freqs_ma, cmap=black_cmap, edgecolors='none')
    plt2 = ax1.pcolor(freqs_rm_ma, cmap=red_cmap, alpha=0.75, edgecolors='none')
    plt3 = ax2.bar(np.arange(L), lanl_freqs.sum(axis=0), width=1, fc='0.5', ec='0.5')
    plt4 = ax2.bar(np.arange(L), lanl_rm.sum(axis=0), width=1, fc='r', ec='none', alpha=0.75)

    plt.setp(ax1.get_xticklabels(), visible=False)
    quart1 = N_samples / 4.
    quart2 = lanl_freqs.shape[0] / 4.
    ax1.set_yticks(np.arange(0, N_samples+quart1, quart1))
    ax2.set_yticks(np.arange(0, N_samples+quart1, quart1))
    ax2.set_yticks(np.arange(0, lanl_freqs.shape[0]+quart2, quart2))
    quarter = N_samples / 4.
    #[ax.set_yticks(np.arange(0, N_samples + quarter, quarter)) for ax in (ax1, ax2)]
    [ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%']) for ax in (ax1, ax2)]
    ax1.set_ylim(0, 1.5*N_samples)
    #ax2.set_ylim(0, 0.8*N_samples)
    ax2.set_ylim(0, 0.9*lanl_freqs.shape[0])

    locs = [(132, 'MA'), (363, 'CA'), (377, 'p2'), (432, 'NC'), (448, 'p1'), (500, 'p6')]
    x_start = 0
    colors = ('#AAAAAA', '#EEEEEE')
    for i, (junc, name) in enumerate(locs):
        color = colors[i%2]
        width = junc - x_start
        rect = Rectangle((x_start, 0), width, 1, edgecolor='None', color=color)
        ax3.add_patch(rect)
        ax3.text(x_start + width/2., 1/2., name, ha='center', va='center')
        ax1.axvline(x_start, color='k', ls='--', lw=1.0)
        ax2.axvline(x_start, color='k', ls='--', lw=1.0)
        x_start = junc
    ax3.set_xlim(0, L)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    [ax.tick_params(top=False, left=False, right=False, bottom=False) for ax in (ax1, ax2, ax3)]

    ax2.invert_yaxis()

    fig.subplots_adjust(hspace=0., wspace=0.02, left=0.08, right=0.95, top=0.95, bottom=0.06)
    if not save:
        plt.show()
    else:
        fig.savefig('PUB_VARIANTS_PLOT.tiff', dpi=600)
        #pdf = PdfPages('/home/wflynn/Documents/publications/gag-pol_deep_sequencing/number_of_variants_vs_position_gag_colorized.pdf')
        #pdf.savefig(fig)
        #pdf.close()

def si_plot(save=False):
    from matplotlib.gridspec import GridSpec

    freqs = get_freqs(gag=False)
    N_samples, L, _ = freqs.shape
    print N_samples, L

    # placeholder for PR resistance mutations
    freqs_rm = freqs.copy()
    for i in range(L):
        for xi, ch in enumerate(ALPHABET[:20]):
            pass

    freqs = np.ma.masked_where(freqs < 0.01, freqs)
    freqs = np.any(freqs >= 0.01, axis=2)
    freqs = freqs.sum(axis=0)
    #rm_freqs = np.ma.masked_where(freqs < 0.01, freqs)
    #rm_freqs = np.any(freqs >= 0.01, axis=2)
    #rm_freqs = freqs.sum(axis=0)

    treated_freqs = np.loadtxt('stanford/pr_freqs_nowt_treated.txt')
    treated_freqs = treated_freqs.reshape(-1, 99, 20)
    N_treated = treated_freqs.shape[0]
    print N_treated

    untreated_freqs = np.loadtxt('stanford/pr_freqs_nowt_untreated.txt')
    untreated_freqs = untreated_freqs.reshape(-1, 99, 20)
    N_untreated = untreated_freqs.shape[0]
    print N_untreated

    #place holder for pr resistance muts
    treated_rm = treated_freqs.copy()
    untreated_rm = untreated_freqs.copy()
    for i in range(L):
        for xi, ch in enumerate(ALPHABET[:20]):
            pass

    treated_freqs = np.ma.masked_where(treated_freqs < 0.01, treated_freqs)
    treated_freqs = np.any(treated_freqs >= 0.01, axis=2)
    treated_freqs = treated_freqs.sum(axis=0)
    #treated_rm = np.ma.masked_where(treated_rm < 0.01, treated_rm)
    #treated_rm = np.any(treated_rm >= 0.01, axis=2)
    #treated_rm = treated_rm.sum(axis=0)
    untreated_freqs = np.ma.masked_where(untreated_freqs < 0.01, untreated_freqs)
    untreated_freqs = np.any(untreated_freqs >= 0.01, axis=2)
    untreated_freqs = untreated_freqs.sum(axis=0)
    #untreated_rm = np.ma.masked_where(untreated_rm < 0.01, untreated_rm)
    #untreated_rm = np.any(untreated_rm >= 0.01, axis=2)
    #untreated_rm = untreated_rm.sum(axis=0)

    fig = plt.figure(figsize=(12, 6))
    matplotlib.rc('font', size=12)

    #fig.suptitle('Variant Frequencies in PR')
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Sequence Position')
    ax1.set_ylabel('Samples with variants')

    x = np.arange(L)
    f1 = freqs/float(N_samples)
    f2 = treated_freqs/float(N_treated)
    f3 = untreated_freqs/float(N_untreated)

    d12 = f1 - f2
    interest_threshold = 0.05
    interest_inds = np.where(d12 > interest_threshold)[0]
    background_inds = list(set(x) - set(interest_inds))

    f11 = f1.copy()
    f12 = f1.copy()
    f21 = f2.copy()
    f22 = f2.copy()
    f31 = f3.copy()
    f32 = f3.copy()

    f11[background_inds] = 0.
    f21[background_inds] = 0.
    f31[background_inds] = 0.
    f12[interest_inds] = 0.
    f22[interest_inds] = 0.
    f32[interest_inds] = 0.

    alpha = 0.7
    c1 = '#666666'
    c2 = '#88CCEE'
    c3 = '#CC6677'
    width1 = 0.8
    width2 = width1/2
    params = {'ec': 'none'}

    plt1 = ax1.bar(x, f11, 1, fc=c1, label='Deep Sequencing', **params)
    plt2 = ax1.bar(x+0.1, f21, width2, fc=c2, label='Stanford Treated', **params)
    plt3 = ax1.bar(x+0.5, f31, width2, fc=c3, label='Stanford Untreated', **params)
    ax1.bar(x, f12, 1, fc=c1, alpha=alpha/2, **params)
    ax1.bar(x+0.1, f22, width2, fc=c2, alpha=alpha, **params)
    ax1.bar(x+0.5, f32, width2, fc=c3, alpha=alpha, **params)

    ax1.set_xticks(range(0,100,20)+[98])
    ax1.set_xticklabels(['1'] + map(str, range(20, 100, 20)) + ['99'])
    ax1.set_yticks([0, .25, .5, .75, 1])
    ax1.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])

    ax1.tick_params(top=False, left=False, right=False, bottom=False)

    ax1.legend()#(plt1, plt2, plt3), ('Deep Sequencing', 'Stanford Treated', 'Stanford Untreated'))

    plt.show()

def publication_plot(save=False):
    from matplotlib.gridspec import GridSpec

    freqs = get_freqs(gag=True)
    N_samples, L, _ = freqs.shape

    freqs_rm = freqs.copy()
    for i in range(L):
        for xi, ch in enumerate(ALPHABET[:20]):
            if str(i+1) + ch in GAG_RESISTANCE_MUTS_SIMP: continue
            freqs_rm[:, i, xi] = 0.

    freqs = np.ma.masked_where(freqs < 0.01, freqs)
    freqs = np.any(freqs >= 0.01, axis=2)
    freqs = freqs.sum(axis=0)
    freqs_rm = np.ma.masked_where(freqs_rm < 0.01, freqs_rm)
    freqs_rm = np.any(freqs_rm >= 0.01, axis=2)
    freqs_rm = freqs_rm.sum(axis=0)

    #lanl_freqs = np.loadtxt('lanl/gag_freqs_nowt_experienced.txt')
    lanl_freqs = np.loadtxt('lanl/gag_freqs_nowt_naive.txt')
    lanl_freqs = lanl_freqs.reshape(-1,500,20)
    N_lanl_samples = lanl_freqs.shape[0]
    lanl_rm = lanl_freqs.copy()
    for i in range(L):
        for xi, ch in enumerate(ALPHABET[:20]):
            if str(i+1) + ch in GAG_RESISTANCE_MUTS_SIMP: continue
            lanl_rm[:, i, xi] = 0.
    lanl_freqs = np.ma.masked_where(lanl_freqs < 0.01, lanl_freqs)
    lanl_freqs = np.any(lanl_freqs >= 0.01, axis=2)
    lanl_freqs = lanl_freqs.sum(axis=0)
    lanl_rm = np.ma.masked_where(lanl_rm < 0.01, lanl_rm)
    lanl_rm = np.any(lanl_rm >= 0.01, axis=2)
    lanl_rm = lanl_rm.sum(axis=0)

    fig = plt.figure(figsize=(12, 6))
    matplotlib.rc('font', size=14)
    #fig.suptitle('Variant Frequencies in Gag')
    gs = GridSpec(3, 2, height_ratios=[10,1,10], width_ratios=[100, 1])
    ax1 = plt.subplot(gs[0, :-1])
    ax3 = plt.subplot(gs[1, :-1])
    ax2 = plt.subplot(gs[2, :-1])
    #matplotlib.rc('axes', linewidth=5)
    [ax.set_xlabel('Sequence Position') for ax in (ax1, ax2)]
    ax1.set_ylabel('Samples with variants')
    ax2.set_ylabel('LANL PI-naive\nsequences with variants')

    x = np.arange(L)
    f1 = freqs/float(N_samples)
    f2 = lanl_freqs/float(N_lanl_samples)
    d = f1 - f2
    interest_threshold = 0.1
    interest_inds = np.where(d > interest_threshold)[0]
    background_inds = list(set(x) - set(interest_inds))
    if 1: # plot with gray background
        f11 = freqs.copy()
        f12 = freqs.copy()
        f13 = freqs_rm.copy()
        f14 = freqs_rm.copy()
        f11[background_inds] = 0.
        f12[interest_inds] = 0.
        f13[background_inds] = 0.
        f14[interest_inds] = 0.

        f21 = lanl_freqs
        f22 = lanl_freqs.copy()
        f23 = lanl_rm.copy()
        f24 = lanl_rm.copy()
        f21[background_inds] = 0.
        f22[interest_inds] = 0.
        f23[background_inds] = 0.
        f24[interest_inds] = 0.

        alpha = 0.2
        plt11 = ax1.bar(x, f11, fc='k')
        plt12 = ax1.bar(x, f12, fc='k', ec='none', alpha=alpha)
        plt13 = ax1.bar(x, f13, fc='r', ec='none')
        plt14 = ax1.bar(x, f14, fc='r', ec='none', alpha=alpha)
        plt21 = ax2.bar(x, f21, fc='k')
        plt22 = ax2.bar(x, f22, fc='k', ec='none', alpha=alpha)
        plt23 = ax2.bar(x, f23, fc='r', ec='none')
        plt24 = ax2.bar(x, f24, fc='r', ec='none', alpha=alpha)
    else: # plot with nothing in between areas of interest
        f11 = freqs.copy()
        f12 = freqs_rm.copy()
        f11[background_inds] = 0.
        f12[background_inds] = 0.

        f21 = lanl_freqs.copy()
        f22 = lanl_rm.copy()
        f21[background_inds] = 0.
        f22[background_inds] = 0.

        plt11 = ax1.bar(x, f11, fc='k')
        plt12 = ax1.bar(x, f12, fc='r', ec='none')
        plt21 = ax2.bar(x, f21, fc='k')
        plt22 = ax2.bar(x, f22, fc='r', ec='none')

    plt.setp(ax1.get_xticklabels(), visible=False)
    quart1 = N_samples / 4.
    quart2 = N_lanl_samples / 4.
    ax1.set_yticks(np.arange(0, N_samples+quart1, quart1))
    ax2.set_yticks(np.arange(0, N_samples+quart1, quart1))
    ax2.set_yticks(np.arange(0, N_lanl_samples+quart2, quart2))
    quarter = N_samples / 4.
    [ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%']) for ax in (ax1, ax2)]
    ax1.set_ylim(0, 0.8*N_samples)
    ax2.set_ylim(0, 0.8*N_lanl_samples)

    locs = [(132, 'MA'), (363, 'CA'), (377, 'p2'), (432, 'NC'), (448, 'p1'), (500, 'p6')]
    x_start = 0
    colors = ('#AAAAAA', '#E1E1E1')
    for i, (junc, name) in enumerate(locs):
        color = colors[i%2]
        width = junc - x_start
        rect = Rectangle((x_start, 0), width, 1, edgecolor='None', color=color)
        ax3.add_patch(rect)
        va = 'center'
        yy = .5
        if name.startswith('p'):
            yy = .425
        ax3.text(x_start + width/2., yy, name, ha='center', va=va, fontsize=12)
        ax1.axvline(x_start, color='k', ls='--', lw=1.0)
        ax2.axvline(x_start, color='k', ls='--', lw=1.0)
        x_start = junc
    ax3.set_xlim(0, L)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    [ax.tick_params(top=False, left=False, right=False, bottom=False) for ax in (ax1, ax2, ax3)]

    ax2.invert_yaxis()

    fig.subplots_adjust(hspace=0., wspace=0.02, left=0.10, right=0.97, top=0.95, bottom=0.12)
    fig.savefig('figure_2.svg', dpi=600)
    plt.show()

if __name__ == "__main__":
    #plot_gag_variants(remove_rm=0, save=0)
    #plot_gag_variants_no_rm(save=0)
    #plot_gag_insets(save=0)
    #plot_pro_variants()
    #plot_cs_variants()
    publication_plot(save=True)
    #si_plot()
