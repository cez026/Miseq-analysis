import os
import sys
import re
import numpy as np
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path


def draw_network(G, pos, ax, sg=None):
    for n in G:
        c = patches.Circle(pos[n], radius=0.02, alpha=0.5)
        ax.add_patch(c)
        G.node[n]['patch'] = c
        x, y = pos[n]

    seen = {}
    for (u, v, d) in G.edges(data=True):
        n1 = G.node[u]['patch']
        n2 = G.node[v]['patch']
        rad = 0.25
        if (u,v) in seen:
            rad = seen.get((u,v))
            rad = (rad + np.sign(rad)*0.1) * -1
        alpha = 0.5
        color = 'k'

        e = patches.FancyArrowPatch(n1.center, n2.center, patchA=n1, patchB=n2,
                                    arrowstyle='-',
                                    connectionstyle='arc3,rad=%s'%rad,
                                    mutation_scale=10.0,
                                    lw=2,
                                    alpha=alpha,
                                    color=color)
        seen[(u,v)] = rad
        ax.add_patch(e)
    return e

def circular_main():
    fig, ax = plt.subplots()

    regiondict = dict(zip(range(1,600), ['MA']*(133-1) + ['CA']*(364-133) + ['p2']*(378-364) + ['NC']*(433-378) + ['p1']*(449-433) + ['p6']*(501-449) + ['PR']*(600-501)))

    N_lines = 50
    muts = get_muts('gag-gag') + get_muts('gag-pr')
    muts = [mut for mut in muts if regiondict[mut[1]] != regiondict[mut[0]]]
    muts.sort(key=lambda x: x[-1], reverse=True)
    min_mi = muts[N_lines-1][2]

    mutdict = {}
    for i in range(1, 600):
        mutdict[i] = {"name": str(i), "size": 1, "imports": [], "region": regiondict[i]}

    for i, j, mi in muts[:N_lines*2]:
        mutdict[i]["size"] += 1
        mutdict[i]["imports"].append(str(j))
        mutdict[j]["size"] += 1
        mutdict[j]["imports"].append(str(i))

    with open('correlations.json', 'w') as fout:
        fout.write('[\n')
        lines = []
        for i, j, mi in muts[:N_lines*2]:
            lines.append('{"source":"%d", "target":"%d", "mi":"%.3f"}'%(i, j, mi))
        fout.write(',\n'.join(lines))
        fout.write('\n]\n')

    with open('corr_data.json', 'w') as fout:
        fout.write('[\n')
        for i in range(1, 600):
            fout.write('{"name":"%s","size":%s,"imports":[%s]},\n'%(mutdict[i]["name"], str(mutdict[i]["size"]), '' if not mutdict[i]["imports"] else '"'+'","'.join(map(str,mutdict[i]["imports"]))+'"'))
        fout.write(']\n')

    with open('cyto_data.csv', 'w') as fout:
        fout.write(','.join(['i', 'ri', 'j', 'rj', 'mi', 's'])+'\n')
        for i, j, mi in muts[:N_lines*2]:
            fout.write(','.join(map(str, [i, regiondict[i], j, regiondict[j], mi, 1]))+'\n')


    mut_pos = []
    for (i,j,mi) in muts[:N_lines]:
        mut_pos.append((i,j))

    G = nx.cycle_graph(599)
    G.add_edges_from(mut_pos)
    pos = nx.circular_layout(G)
    print G.number_of_edges()

    draw_network(G, pos, ax)

    ax.autoscale()
    plt.show()

def get_muts(region='gag-gag'):
    muts = []
    infile = '/home/wflynn/code/mi/mw/%s_MI_MW.txt'%region
    r1, r2 = region.split('-')
    max_mi = 0
    with open(infile, 'r') as f:
        f.next() # skip header
        for line in f:
            #line = re.split('- |\t |, ', line.rstrip())
            line = line.rstrip().split('\t')
            i, j = map(int, line[0].split('-'))
            if r1 == 'pr': i += 500
            if r2 == 'pr': j += 500
            mi = float(line[1])
            if mi > max_mi: max_mi = mi
            muts.append((i, j, mi))
    print max_mi
    return muts

def main():
    matplotlib.rc('font', size=12)
    fig = plt.figure(figsize=(16,9))
    gs = GridSpec(2, 1, height_ratios=[20, 1])#, 20])
    gs.update(hspace=0., wspace=0.)

    ax1 = plt.subplot(gs[0])
    label_ax = plt.subplot(gs[1])

    [ax.set_xlim(0, 599) for ax in (ax1, label_ax)]
    ax1.set_ylim(0, 350)

    # New way
    regiondict = dict(zip(range(1,600), ['MA']*(133-1) + ['CA']*(364-133) + ['p2']*(378-364) + ['NC']*(433-378) + ['p1']*(449-433) + ['p6']*(501-449) + ['PR']*(600-501)))

    N_lines = 50
    muts = get_muts('gag-gag') + get_muts('gag-pr')
    muts = [mut for mut in muts if regiondict[mut[1]] != regiondict[mut[0]]]
    muts.sort(key=lambda x: x[-1], reverse=True)
    min_mi = muts[N_lines-1][2]

    counter = 0
    for mut in muts[:N_lines]:
        r1, r2 = regiondict[mut[0]], regiondict[mut[1]]
        c = 'r' if r2 == 'PR' else 'b'
        ax1.add_patch(make_circ(*mut, ec=c))
        counter += 1

    print counter

    r = range(1)
    proxy1 = plt.Line2D(r, r, color='b', markerfacecolor='none', lw=3)
    proxy2 = plt.Line2D(r, r, color='r', markerfacecolor='none', lw=3)
    ax1.legend((proxy1, proxy2), ('Gag-Gag', 'Gag-PR'))

    # Add x-axis boxes
    locs = [(132, 'MA'), (363, 'CA'), (377, 'p2'), (432, 'NC'), (448, 'p1'), (500, 'p6'),
            (599, 'PR')]
    x_start = 0
    colors = ('#AAAAAA', '#EEEEEE')
    for i, (junc, name) in enumerate(locs):
        color = colors[i%2]
        width = junc - x_start
        rect = patches.Rectangle((x_start, 0), width, 1, color=color)
        label_ax.add_patch(rect)
        label_ax.text(x_start + width/2., 1/2., name, ha='center', va='center')
        x_start = junc
    label_ax.set_xlim(0, 599)
    label_ax.set_xticks([1]+range(50, 650, 50))
    label_ax.set_xticklabels([1]+range(50, 500, 50)+[1]+range(50, 150, 50))

    [plt.setp(ax.get_xticklabels(), visible=False) for ax in (ax1, )]
    [plt.setp(ax.get_yticklabels(), visible=False) for ax in (ax1, label_ax)]
    [ax.tick_params(top=False, left=False, right=False, bottom=False) for ax in (ax1, label_ax)]
    ax1.tick_params(bottom=True)
    ax1.set_xticks(np.arange(0, 599, 10))

    label_ax.set_xlabel('Sequence position')

    plt.show()

if __name__ == "__main__":
    circular_main()
