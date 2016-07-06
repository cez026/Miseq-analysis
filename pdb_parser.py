import numpy as np
from itertools import combinations, product
from Bio.PDB import PDBParser

class PDBModel(object):
    def __init__(self, model):
        self.model = model
        self.model_id = self.model.get_id()
        self.parent = self.model.get_parent()
        print "Model %i of structure %s loaded."%(self.model_id, self.parent)

        self.residues = list(self.model.get_residues())

    def get_single_distance(self, p1, p2, exclude_backbone=False,
            CA_only=False):
        error_return = {'dist': 100, 'chains': 'XX', 'r1': 'X', 'r2': 'X',
                        'model_id': self.model_id, 'atoms': 'X-X'}

        res1, res2 = [], []
        for res in self.residues:
            if res.get_id()[1] == p1:
                res1.append(res)
            elif res.get_id()[1] == p2:
                res2.append(res)

        min_dist = 100
        min_chain = 'AA'
        min_atoms = ""
        if res1 == [] or res2 == []:
            print "Problem with %i %i in %s: Model %s"%(p1, p2, self.parent,
                                                        self.model_id)
            return error_return

        for r1, r2 in product(res1, res2):
            for atom1 in r1.child_list:
                if 'H' in atom1.get_name(): continue
                if exclude_backbone and atom1.get_name() in ('C', 'N'): continue
                if CA_only and atom1.get_name() != 'CA': continue

                coord1 = atom1.get_coord()
                for atom2 in r2.child_list:
                    if 'H' in atom2.get_name(): continue
                    if exclude_backbone and atom2.get_name() in ('C', 'N'): continue
                    if CA_only and atom2.get_name() != "CA": continue

                    coord2 = atom2.get_coord()

                    dist = np.sqrt(((coord1 - coord2)**2).sum())
                    if dist < min_dist:
                        min_dist = dist
                        min_chain = r1.get_parent().id + r2.get_parent().id
                        min_atoms = "%s-%s"%(atom1.get_name(), atom2.get_name())

        return {'dist': min_dist, 'chains': min_chain,
                'r1': r1.get_resname(), 'r2': r2.get_resname(),
                'model_id': self.model_id, 'atoms': min_atoms}

    def get_pair_distances(self):
        distances = {}

        for resi, resj in combinations(self.residues, 2):
            posi = resi.get_id()[1]
            posj = resj.get_id()[1]
            chaini = resi.get_parent().id
            chainj = resj.get_parent().id
            pair = (posi, posj)

            min_dist = 100
            for atomi in resi.child_list:
                coordi = atomi.get_coord()
                for atomj in resj.child_list:
                    coordj = atomj.get_coord()
                    dist = np.sqrt(((coordi - coordj)**2).sum())
                    if dist < min_dist:
                        min_dist = dist
            if min_dist == 100:
                print "Bad info for %s-%s: %i %i"%(self.parent, self.model_id, resi, resj)

            value = (min_dist, chaini, chainj, resi.get_resname(),
                     resj.get_resname())
            distances[pair] = value

        return distances


class PDBAtomAtomDistanceReader(object):
    def __init__(self, pdbname, label):
        self.pdbname = pdbname
        self.label = label
        self.region = label.translate(None, '0123456789')

        self.parser = PDBParser(QUIET=True)
        self.structure = self._load_structure()
        self.models = [PDBModel(model) for model in self.structure.get_list()]

    def _load_structure(self):
        return self.parser.get_structure(self.label, self.pdbname)

    def get_single_distance(self, p1, p2, exclude_backbone=False,
            CA_only=False):
        distances = []

        for model in self.models:
            dist_info = model.get_single_distance(p1, p2, exclude_backbone,
                                                  CA_only)
            distances.append(dist_info)

        distances.sort(key=lambda x: x['dist'])
        min_dist = distances[0]
        max_dist = distances[-1]
        avg_dist = np.array([d['dist'] for d in distances]).mean()

        return {'protein': self.region, 'p1': p1, 'p2': p2, 'label': self.label,
                'r1': min_dist['r1'], 'r2': min_dist['r2'],
                'min_dist': min_dist,
                'avg_dist': avg_dist,
                'max_dist': max_dist}

    def get_pair_distances(self):
        distances = []

        print "Calculating distances..."
        for model in self.models:
            dist_info = model.get_pair_distances()
            distances.append(dist_info)

        if len(self.models) == 1:
            return [k + v for k, v in distances[0].iteritems()]

        final_distances = {}

        residues = self.models[0].residues
        for r1, r2 in combinations(residues, 2):
            p1 = r1.get_id()[1]
            p2 = r2.get_id()[1]
            pair = (p1, p2)
            final_distances[pair] = min([d[pair] for d in distances])

        final_distances = [k + v for k, v in final_distances.iteritems()]

        return final_distances

def tests():
    s = PDBAtomAtomDistanceReader('1ODW_pr.pdb', 'PR')
    s.get_single_distance(30, 88)

if __name__ == "__main__":
    tests()

