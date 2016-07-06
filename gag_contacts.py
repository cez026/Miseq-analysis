import csv
from pdb_parser import PDBAtomAtomDistanceReader

POS2GAG = dict(zip(range(1, 501), ['MA']*500))
for i in range(133, 501):
    if i > 448: POS2GAG[i] = 'p6'
    elif i > 432: POS2GAG[i] = 'p1'
    elif i > 377: POS2GAG[i] = 'NC'
    elif i > 363: POS2GAG[i] = 'p2'
    else: POS2GAG[i] = 'CA'

def main():
    gag_gag_file = 'gag-gag_correlations_mw.csv'
    pr_pr_file = 'pr-pr_correlations_mw.csv'

    #pdb_MA = "2H3F_ma_1chain.pdb"
    pdb_MA = "2H3F_ma.pdb"
    pdb_MA3 = "1HIW_ma_trimer.pdb"
    pdb_CA = "3MGE_ca.pdb"
    pdb_CA2 = "2M8L_ca_dimer.pdb"
    pdb_CA5 = "3P05_ca_pentamer.pdb"
    pdb_CA6 = "3MGE_ca_hexamer.pdb"
    #pdb_NC = "2EXF_nc_1chain.pdb"#"1ESK_nc_1chain.pdb"
    pdb_NC = "2EXF_nc.pdb"#"1ESK_nc_1chain.pdb"
    pdb_PR2 = "1ODW_pr.pdb"

    regions = ("MA", "CA", "NC", "PR")

    print "Reading structures..."
    pdbwrapper = lambda x: PDBAtomAtomDistanceReader(*x)

    # 1 only if monomer structures
    if 0:
        structures = dict(zip(regions,
                              [(pdbwrapper((pdb_MA, 'MA')),),
                               (pdbwrapper((pdb_CA, 'CA')),),
                               (pdbwrapper((pdb_NC, 'NC')),),
                               (pdbwrapper((pdb_PR2, 'PR')),)]))
    else:
        structures = dict(zip(regions,
                              [map(pdbwrapper, zip((pdb_MA, pdb_MA3), ('MA', 'MA3'))),
                               map(pdbwrapper, zip((pdb_CA, pdb_CA2, pdb_CA5,
                                   pdb_CA6), ('CA', 'CA2', 'CA5', 'CA6'))),
                               map(pdbwrapper, zip((pdb_NC,), ('NC',))),
                               map(pdbwrapper, zip((pdb_PR2,), ('PR2',)))]))

    # PDBAtomAtomDistanceReader.get_single_distance() output
    #return {'protein': self.region, 'p1': p1, 'p2': p2, 'label': self.label,
    #        'r1': min_dist['r1'], 'r2': min_dist['r2'],
    #        'min_dist': min_dist,
    #        'avg_dist': avg_dist,
    #        'max_dist': max_dist}
    pair_distances = [('Protein', 'Pos 1', 'Pos 2', 'Res 1', 'Res 2', 'MI',
                       'Smallest Rij', 'Structure', 'Model', 'Chains',
                       'Atoms')]
    structure_report = []
    for region, pairsfile in zip(('Gag','PR'), (gag_gag_file, pr_pr_file)):
        print "Intra %s distances for strongest %s correlations"%(region,
                                                                  region)

        # Read correlated pair files and get intra-protein pairs
        correlated_pairs = []
        with open(pairsfile, 'r') as f:
            f.next()
            f.next()
            for line in f:
                if region == "PR":
                    p1, p2, _, _, _, mi, _ = line.split(',', 6)
                    r1 = r2 = "PR"
                    p1, p2 = map(int, (p1, p2))
                else:
                    p1, p2, _, _, mi, _ = line.split(',', 5)
                    p1, p2 = map(int, (p1, p2))
                    r1 = POS2GAG[p1]
                    r2 = POS2GAG[p2]
                true_p1, true_p2 = p1, p2

                # exclude inter-protein pairs
                if r1 != r2: continue

                # change protein numbering for CA and NC
                if r1 == "CA":
                    p1 -= 132
                    p2 -= 132
                elif r2 == "NC":
                    p1 -= 363
                    p2 -= 363

                correlated_pairs.append((true_p1, true_p2, p1, p2, r1, mi))

        for tp1, tp2, p1, p2, r, mi in correlated_pairs:
            if r in ('p1', 'p2', 'p6'): continue
            #struct = structures[r]
            print "Current pair: %s %i %i"%(r, p1, p2)
            structlist = structures[r]
            minimum = {'min_dist': {'dist': 100}}
            for struct in structlist:
                dist = struct.get_single_distance(p1, p2, exclude_backbone=True,
                                                  CA_only=False)
                if dist['min_dist']['dist'] < minimum['min_dist']['dist']:
                    minimum = dist

                s_result = (dist['protein'], dist['label'], tp1, tp2,
                            dist['p1'], dist['p2'], dist['r1'], dist['r2'], mi,
                            dist['min_dist']['dist'], dist['avg_dist'],
                            dist['max_dist']['dist'])
                structure_report.append(s_result)

            result = (minimum['protein'], tp1, tp2,
                      minimum['r1'], minimum['r2'], mi,
                      minimum['min_dist']['dist'], minimum['label'],
                      minimum['min_dist']['model_id'],
                      minimum['min_dist']['chains'],
                      minimum['min_dist']['atoms'])

            pair_distances.append(result)

    pair_distances.sort(key=lambda x: x[5], reverse=1)
    with open('correlated_pair_distances.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(pair_distances)

    structure_report.sort()
    structure_report.insert(0, ('Protein', 'Structure',
                                'Gag Pos 1', 'Gag Pos 2','Pos 1', 'Pos 2',
                                'Res 1', 'Res 2', 'MI', 'Min Rij', 'Avg Rij',
                                'Max Rij'))
    with open('pdb_structure_report.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(structure_report)

def odds_test():
    def print_results(d, name):
        N_pairs = len(d)
        close_pairs = len([p for p in d if p[2] < 8])
        print name
        print "total\t<8\tprob"
        print N_pairs, close_pairs, float(close_pairs) / N_pairs

    pdb_PR2 = "1ODW_pr.pdb"
    pdb_CA5 = "3P05_ca_pentamer.pdb"
    pdb_MA = "2H3F_ma.pdb"

    pr_structure = PDBAtomAtomDistanceReader(pdb_PR2, 'PR2')
    ca5_structure = PDBAtomAtomDistanceReader(pdb_CA5, 'CA5')
    ma_structure = PDBAtomAtomDistanceReader(pdb_MA, 'MA')

    print_results(pr_structure.get_pair_distances(), 'PR dimer')
    print_results(ma_structure.get_pair_distances(), 'MA monomer')
    print_results(ca5_structure.get_pair_distances(), 'CA pentamer')

if __name__ == "__main__":
    #main()
    odds_test()

