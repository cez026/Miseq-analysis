#!/usr/bin/env python
'''
Created on Mar 7, 2013

@author: Bill
'''

import os
import sys
import csv
import pysam
import argparse
import numpy as np
from Bio.Seq import translate
from tempfile import NamedTemporaryFile
from itertools import permutations, product, combinations

NT_REF_PATH = os.sep.join([os.getcwd(), '..', 'samples', 'consensus.fa'])
with open(NT_REF_PATH, 'r') as ref_file:
    GAGPOL_NT_REFERENCE = ref_file.readlines()[1:]
    GAGPOL_NT_REFERENCE = ''.join(  map(lambda x: x.strip('\n'), 
                                    GAGPOL_NT_REFERENCE))

PRO_START = 1463
PRO_END = 1760
GAG_START = 0
GAG_END = 1498

PRO_NT_REFERENCE = GAGPOL_NT_REFERENCE[PRO_START:PRO_END + 1]
PRO_AA_REFERENCE = translate(PRO_NT_REFERENCE)
GAG_NT_REFERENCE = GAGPOL_NT_REFERENCE[GAG_START:GAG_END + 1]
GAG_AA_REFERENCE = translate(GAG_NT_REFERENCE)
#print GAG_AA_REFERENCE

class Analyzer:
    def __init__(self, sam_path, use_temp=True, sortn=False):
        print "Name: %s"%sam_path
        self.samdirname = os.path.splitext(sam_path)[0]
        self.patientid = sam_path.split(os.sep)[-1][:-4]
        print 'dirname: %s'%self.samdirname
        self.pro_reads = self.samdirname + '_pro_reads'
        self.pro_counts = self.samdirname + '_pro_counts'
        self.gag_counts = self.samdirname + '_gag_counts'
        self.pro_pair_counts = self.samdirname + '_pro_pair_counts'
        self.pro_plts = self.samdirname + '_pro_plts'
        self.gag_plts = self.samdirname + '_gag_plts'
        self.sig_pro_plts = self.samdirname + '_sig_pro_plts'
        self.sig_gag_plts = self.samdirname + '_sig_gag_plts'
        self.pileup_name = self.samdirname + '_pileup'

        bam_extension = '_sorted'
        self.use_temp = use_temp

        if use_temp:
            if sam_path.endswith('.sam'):
                print "make bam...",
                bam_content = pysam.view('-bS', sam_path)

                # write BAM to a temp file
                self.temp_file = NamedTemporaryFile(delete=False)
                self.temp_filename = self.temp_file.name
                self.temp_file.writelines(bam_content)
                self.temp_file.close()

                # sort BAM file
                print "sort...",
                pysam.sort(self.temp_file.name, 
                           self.temp_file.name+bam_extension)
                print "index...",
                pysam.index('%s%s.bam' % (self.temp_file.name, bam_extension))
                print "make sam!"

                self.samfile = pysam.Samfile(self.temp_file.name
                                             + bam_extension +'.bam', 'rb')
            else:
                self.use_temp = False
                if sortn:
                    sorted_path = sam_path + '_nsorted'
                    if not os.path.exists(sorted_path):
                        print "sorting by query name"
                        pysam.sort('-n', sam_path, sorted_path)
                    self.samfile = pysam.Samfile(sorted_path, 'rb')
                else:
                    self.samfile = pysam.Samfile(sam_path, 'rb')
        else:
            print 'storing bam files'
            if sam_path.endswith('.sam'):
                print "make bam...",
                bam_content = pysam.view('-bS', sam_path)

                # write BAM to a temp file
                self.bam_file_name = self.samdirname+'.bam'
                self.bam_file = open(self.bam_file_name, 'w+')
                self.bam_file.writelines(bam_content)
                self.bam_file.close()

                # sort BAM file
                print "sort...",
                pysam.sort(self.bam_file_name, self.bam_file_name+bam_extension)
                print "index...",
                pysam.index('%s%s.bam' % (self.bam_file_name, bam_extension))
                print "make sam!"

                self.samfile = pysam.Samfile(self.bam_file_name
                                             + bam_extension +'.bam', 'rb')
            else:
                if sortn:
                    sorted_path = sam_path + '_nsorted'
                    if not os.path.exists('%s.bam'%sorted_path):
                        print "sorting by query name..."
                        pysam.sort('-n', sam_path, sorted_path)
                    self.samfilen = pysam.Samfile('%s.bam'%sorted_path, 'rb')
                self.samfile = pysam.Samfile(sam_path, 'rb')

    def __del__(self):
        if self.use_temp:
            os.unlink(self.temp_filename)
            os.unlink(self.temp_filename + '_sorted.bam')
            os.unlink(self.temp_filename + '_sorted.bam.bai')

    
    def sam_stats(self):
        mapped = self.samfile.mapped
        unmapped = self.samfile.unmapped
        total = float(mapped + unmapped)

        print 'filename, mapped, unmapped, percent mapped'
        print '%s, %d, %d, %.2f%% map'%(self.samfile.filename, mapped, unmapped,
                                        100 * mapped/total)


    def sam_coverage(self):
        # This process doesn't work properly because samtools limits the max 
        # read depth to 8000 (or so) reads. The pysam committers said it's 
        # samtools, not pysam, that's the problem.
        pileup_iter = self.samfile.pileup('CONSENSUS_B_GAG_POL', GAG_START, PRO_END, maxdepth=1e6)
        return [p.n for p in pileup_iter]

    def trim_read(self, read, start, end, codon=True):
        """
        M    BAM_CMATCH         0
        I    BAM_CINS           1
        D    BAM_CDEL           2
        N    BAM_CREF_SKIP      3
        S    BAM_CSOFT_CLIP     4
        H    BAM_CHARD_CLIP     5
        P    BAM_CPAD           6
        =    BAM_CEQUAL         7
        X    BAM_CDIFF          8
        """
    
        if read.pos > end or read.aend < start:
            if codon: 
                return '', 0
            else:
                return ''
    
        aligned_seq = ''
        read_pos = 0
        for code, n in read.cigar:
            if code == 7:
                raise Exception(KeyError, "Exact match?")
            if code == 0:
                aligned_seq += read.seq[read_pos:read_pos + n]
            if code == 1:
                pass
            if code == 2:
                aligned_seq += 'N' * n
                read_pos -= n
            if code == 3:
                raise Exception(KeyError, "This shouldn't happen...")
            if code == 4:
                pass
            if code == 5:
                pass
            read_pos += n
    
        trimmed_seq = aligned_seq
        l_offset = start - read.pos
        r_offset = read.pos + len(aligned_seq) - end
        frame_offset = 0
        if l_offset > 0:
            trimmed_seq = trimmed_seq[l_offset:]
        if r_offset > 0:
            trimmed_seq = trimmed_seq[0:len(trimmed_seq) - r_offset]
        if not codon:
            return trimmed_seq

        if l_offset < 0:
            frame_offset = (start - read.pos) % 3
        return trimmed_seq, frame_offset


    def translate_read(self, read, start, end):
        trimmed_read, offset = self.trim_read(read, start, end, codon=True)
        prot_seq = translate(trimmed_read[offset:])
    
        prot_start = (read.pos + offset - start) / 3
        if prot_start < 0:
            prot_start = 0
    
        return prot_seq, prot_start

#------------------------------------------------------------------------------

    def _dNdS_sites(self, codon):
        syn = 0
        non = 0
        alphabet = 'ACGT'
        aa = translate(codon)
        if len(codon) < 3: return syn, non

        for i in range(3):
            for mut in alphabet:
                if mut == codon[i]: continue
                mut_codon = codon[:i] + mut + codon[i+1:]

                syn_flag = (aa == translate(mut_codon))
                syn += syn_flag
                non += (not syn_flag)
        
        syn /= 3.
        non /= 3.
        assert syn + non == 3
        return syn, non

    def dNdS(self, reference, start, end):
        ps, pn = 0, 0
        n_codon = (end - start) / 3
        for i in range(start, end, 3):
            ref_codon = reference[i-start:i-start+3]
            ref_aa = translate(ref_codon)
            s_i, n_i = self._dNdS_sites(ref_codon)
            if s_i == 0: continue
            
            m_i = 0
            inner_s, inner_n = 0, 0
            reads = self.samfile.fetch('CONSENSUS_B_GAG_POL', start, end)
            for read in reads:
                trimmed_read = self.trim_read(read, i, i+3, codon=False)
                if len(trimmed_read) < 3: continue

                m_i += 1
                cur_pos = read.pos - i
                if cur_pos < 0: 
                    cur_pos = 0

                sij, nij = 0, 0
                for j, nt in enumerate(trimmed_read):
                    if nt == ref_codon[j]: continue
                    mut_codon = ref_codon[:j] + nt + ref_codon[j+1:]
                    if translate(mut_codon) == ref_aa:
                        sij += 1
                    else:
                        nij += 1
            
                inner_s += sij / s_i
                inner_n += nij / n_i
            
            ps += inner_s / m_i
            pn += inner_n / m_i

        ps /= float(n_codon)
        pn /= float(n_codon)

        ds = -.75 * np.log(1 - 4*ps/3)
        dn = -.75 * np.log(1 - 4*pn/3)
        print ds/dn

#------------------------------------------------------------------------------

    def nucleotide_counts(self, reference, start, end):
        reads = self.samfile.fetch('CONSENSUS_B_GAG_POL', start, end)

        mutations = []
        for nt in reference:
            mutations.append(dict(zip(  ('ref','A','C','G','T','N'), 
                                        (nt, 0, 0, 0, 0, 0))))

        for read in reads:
            trimmed_read = self.trim_read(read, start, end, codon=False)
            if trimmed_read == '': continue
            cur_pos = read.pos - start
            if cur_pos < 0:
                cur_pos = 0 
            for nt in trimmed_read:
                if nt not in ('A', 'C', 'G', 'T', 'N'):
                    pass
                else:
                    mutations[cur_pos][nt] = mutations[cur_pos].get(nt, 0) + 1
                cur_pos += 1

        return mutations

    def protein_counts(self, reference, start, end):
        reads = self.samfile.fetch('CONSENSUS_B_GAG_POL', start, end)
        
        mutations = []
        for aa in reference:
            d = dict(zip(   ('ACDEFGHIKLMNPQRSTVWYX*'),
                            np.zeros(22)))
            d['ref'] = aa
            mutations.append(d)

        for read in reads:
            trans_read, trans_start = self.translate_read(read, start, end)
            if trans_read == '': continue
            cur_pos = trans_start

            for aa in trans_read:
                mutations[cur_pos][aa] = mutations[cur_pos].get(aa, 0) + 1
                cur_pos += 1

        return mutations

    def export_reads(self, reference, start, end):
        reads = self.samfile.fetch('CONSENSUS_B_GAG_POL', start, end)

        with open(self.pro_reads, 'w') as outf:
            for read in reads:
                trans_read, trans_start = self.translate_read(read, start, end)
                if trans_read == '': continue
                outf.write('%d,%s\n'%(trans_start, trans_read))

    def export_sequences2(self, reference, start, end):
        def write_se(read, start, end, outf):
            L = (end - start) / 3
            read1, start1 = self.translate_read(read, start, end)
            if read1 == '': return
            len1 = len(read1)
            seq = '.'*start1 + read1 + '.'*(L-len1-start1)
            outf.write("%s\n"%seq)

        def write_pe(read, mate, start, end, outf):
            L = (end - start)/3
            if read.qname != mate.qname:
                write_se(read, start, end, outf)
                write_se(mate, start, end, outf)

            read1, s1 = self.translate_read(read, start, end)
            read2, s2 = self.translate_read(mate, start, end)
            if read1 == '' and read2 == '': return 
            if s1 > s2:
                read1, read2 = read2, read1
                s1, s2 = s2, s1

            len1, len2 = len(read1), len(read2)
            if s2 >= s1 + len1:
                if s2 > L: s2 = L
                seq = '.'*s1 + read1 + '.'*(s2-len1-s1) + read2 + '.'*(L-len2-s2)
            else:
                seq = '.'*s1 + read1[:s2-s1] + read2
                seq += '.'*(L-len(seq))
            outf.write("%s\n"%seq)

        #L = (end - start) / 3

        #count = 0
        #found_mate = True
        mate1, mate2 = None, None
        with open(self.pro_reads, 'w') as outf:
            for read in self.samfilen:
                if not read.is_proper_pair or not read.is_paired:
                    write_se(read, start, end, outf)
                elif read.is_proper_pair and read.is_read1:
                    if mate1 is not None: write_se(mate1, start, end, outf)
                    mate1 = read
                elif read.is_proper_pair and read.is_read2:
                    mate2 = read
                    if mate1 and mate2:
                        write_pe(mate1, mate2, start, end, outf)
                    else:
                        write_se(mate2, start, end, outf)
                    mate1, mate2 = None, None

                ## get read1
                ## if previous read1 didn't get a pair, write it out
                #if read.is_proper_pair and read.is_read1:
                #    if not found_mate:
                #        found_mate = True
                #        if read1 == '': continue
                #        len1 = len(read1)
                #        seq = '.'*start1 + read1 + '.'*(L-len1-start1)
                #        outf.write('%s\n'%seq)
                #    else:
                #        read1, start1 = self.translate_read(read, start, end)
                #        found_mate = False
                #        continue
                ## get read2
                #elif read.is_proper_pair and read.is_read2:
                #    found_mate = True
                #    read2, start2 = self.translate_read(read, start, end)
                #    if read1 == '' and read2 == '': continue
                #    if start2 < start1:
                #        read1, read2 = read2, read1
                #        start1, start2 = start2, start1

                #    len1 = len(read1)
                #    len2 = len(read2)
                #    # read2 is separated from read 1
                #    if start2 >= start1 + len1:
                #        if start2 > L: start2 = L
                #        seq = '.'*start1 + read1 + '.'*(start2-len1-start1) +\
                #              read2 + '.'*(L-len2-start2)
                #    # read2 and read1 overlap
                #    else:
                #        seq = '.'*start1 + read1[:start2-start1] + read2
                #        seq += '.'*(L-len(seq))
                #    
                #    outf.write('%s\n'%seq)

                #elif not read.is_proper_pair:
                #    read1, start1 = self.translate_read(read, start, end)
                #    found_mate = True
                #    if read1 == '': continue
                #    len1 = len(read1)
                #    seq = '.'*start1 + read1 + '.'*(L-len1-start1)
                #    
                #outf.write('%s\n'%seq)

    def export_sequences(self, reference, start, end):
        reads = self.samfile.fetch('CONSENSUS_B_GAG_POL', start, end, until_eof=1)
        L = (end - start) / 3

        count = 0
        with open(self.pro_reads, 'w') as outf:
            for read in reads:
                # incorporate paired reads
                if read.is_proper_pair and read.is_read1:
                    pointer = self.samfile.tell() # save current position
                    try: mate = self.samfile.mate(read)
                    except ValueError: continue
                    finally:
                        self.samfile.seek(pointer)
                        read1, start1 = self.translate_read(read, start, end)
                        read2, start2 = self.translate_read(mate, start, end)
                        
                        if start2 < start1:
                            read1, read2 = read2, read1
                            start1, start2 = start2, start1

                        len1 = len(read1)
                        len2 = len(read2)
                            
                        seq = '.'*start1 + read1 + '.'*(start2-len1-start1) +\
                              read2 + '.'*(L-len2-start2)
                    
                    outf.write('%s\n'%seq)
                    count += 1
                    if count%1000==0: print count

                elif not read.is_proper_pair:
                    read1, start1 = self.translate_read(read, start, end)
                    if read1 == '': continue
                    len1 = len(read1)
                    seq = '.'*start1 + read1 + '.'*(L-len1)
                    
                    outf.write('%s\n'%seq)
                    count += 1
                    if count%1000==0: print count

    def cwr(iterable, r):
        # combinations_with_replacement (itertools 2.7 generator)
        pool = tuple(iterable)
        n = len(pool)
        for indices in product(range(n), repeat=r):
            if sorted(indices) == list(indices):
                yield typle(pool[i] for i in indices)

    def protein_pair_counts(self, reference, start, end):
        # THIS SUCKS AND IS WAY TOO SLOW
        reads = self.samfile.fetch('CONSENSUS_B_GAG_POL', start, end)

        mutations = []
        all_aas = 'ACDEFGHIKLMNPQRSTVQYX*'
        possible_combos = [''.join(aas) for aas in product(all_aas, repeat=2)]
        for aa_pair in combinations(reference, 2):
            d = dict(zip(possible_combos, [0]*len(possible_combos)))
            d['ref'] = ''.join(aa_pair)
            mutations.append(d)

        for read in reads:
            # If its read one and part of a pair, try to get bivariates from 
            # pair
            if read.is_proper_pair and read.is_read1:
                pointer = self.samfile.tell()
                try:
                    mate = self.samfile.mate(read)
                except ValueError:
                    continue
                finally:
                    self.samfile.seek(pointer)

                    read1, start1 = self.translate_read(read, start, end)
                    read2, start2 = self.translate_read(mate, start, end)
                    if read1 == '' or read2 == '':
                        pass#continue

                    # Ensure read1 starts before read2
                    if start2 < start1:
                        swpread, swpstart = read2, start2
                        read2, start2 = read1, start1
                        read1, start1 = swpread, swpstart

                    cur_pos = len([j for i in range(start1) 
                                     for j in range(len(reference)-i)])\
                                     + start2

                    for i in range(len(read1)):
                        for j in range(min(i+1, start2), len(read2)):
                            pair = read1[i] + read2[j]
                            mutations[cur_pos][pair] = mutations[cur_pos].get(
                                                        pair, 0) + 1
                            cur_pos += 1
                        cur_pos += len(reference) - i + start2 
            
            # Regardless of what read it is, we want the bivariates from just
            # the single read. The mate to this read will get its turn when
            # it fails the is_proper_pair/is_read1 if above. This catches reads
            # that are unpaired. 
            read1, start1 = self.translate_read(read, start, end)
            cur_pos = len([j for i in range(start1) 
                             for j in range(len(reference)-i)])

            for i in range(len(read1)):
                for j in range(i+1, len(read1)):
                    pair = read1[i] + read1[j]
                    mutations[cur_pos][pair] = mutations[cur_pos].get(pair,
                                                                        0) + 1
                    cur_pos += 1
                cur_pos += len(reference) - i
            
        return mutations


#------------------------------------------------------------------------------


    def get_local_codon(self, mutations, pos, mut=None):
        codon_pos = pos%3
        codon_seq = mutations[pos - codon_pos:pos + 3 - codon_pos]
        codon_seq = ''.join(map(lambda _x: _x['ref'], codon_seq))
        if mut:
            codon_seq = codon_seq[:codon_pos]+ mut + codon_seq[codon_pos + 1:]
        return translate(codon_seq)

    def export_nucleotide_frequencies(self, mutations, outfilename, start,
                                      threshold=None):
        outputfile = open(outfilename, 'w')
        writer = csv.writer(outputfile)

        sig_freqs = []

        N = len(mutations)
        for pos in range(N):
            pos_info = mutations[pos].copy()
            pos_ref = pos_info.pop('ref')
            pos_nts = pos_info.items()
            
            # synonymous mutation info
            ref_codon = self.get_local_codon(mutations, pos)
            
            total = float(sum(pos_info.values()))
            if total == 0: total = 1.
            for nt, count in pos_nts:
                freq = count/total
                wt_flag = int(nt == pos_ref)
                #if not threshold:
                writer.writerow([pos, nt, count, "%.8f"%freq, wt_flag])
                #else:
                #    if wt_flag == False and freq > threshold:
                #        mut_codon = self.get_local_codon(mutations, pos, mut=nt)
                #        syn_flag = int(mut_codon == ref_codon)
                #        writer.writerow([pos+start, nt, count, "%.8f"%freq, 
                #                        'ref:'+pos_ref,
                #                        'aa:%s->%s'%(ref_codon, mut_codon)])
                #        sig_freqs.append({'pos': pos, 'freq':freq,
                #                          'mut_nt':nt, 'mut_aa':mut_codon,
                #                          'ref_nt':pos_ref, 'ref_aa':ref_codon})
        outputfile.close()
        return sig_freqs

    def export_amino_acid_frequencies(self, mutations, outfilename, start):
        print "exporting amino acids to %s"%outfilename
        outputfile = open(outfilename, 'wb')
        writer = csv.writer(outputfile)
        bad_codons = set(['*', 'X'])

        #sig_freqs []
        N = len(mutations)
        for pos in range(N):
            pos_info = mutations[pos].copy()
            pos_ref = pos_info.pop('ref')
            pos_aas = pos_info.items()

            total = float(sum(pos_info.values()))
            if total == 0: total = 1.
            for aa, count in pos_aas:
                freq = count/total
                wt_flag = int(aa == pos_ref)
                writer.writerow([pos, aa, count, "%.8f"%freq, wt_flag])
        outputfile.close()

    def export_amino_acid_pair_frequencies(self, mutations, outfilename):
        print "exporting amino acid pair counts to %s"%outfilename
        outputfile = open(outfilename, 'wb')
        writer = csv.writer(outputfile)

        N = len(mutations)
        for pos in range(N):
            pos_info = mutations[pos].copy()
            pos_ref = pos_info.pop('ref')
            pos_ass = pos_info.items()

            total = float(sum(pos_info.values()))
            if total < 1: total = 1.
            for aa, count in pos_aas:
                freq = count/total
                wt_flag = int(aa == pos_ref)
                writer.writerow([pos, aa, count, "%.8f"%freq, wt_flag])
        outputfile.close()

    def export_protease_reads(self):
        self.export_reads(PRO_AA_REFERENCE, PRO_START, PRO_END)
    def export_protease_sequences(self):
        self.export_sequences2(PRO_AA_REFERENCE, PRO_START, PRO_END)

#------------------------------------------------------------------------------

    def plot_nucleotide_frequencies(self, outfile, mutations, start, end):
        import matplotlib.pyplot as plt
        N = len(mutations)
        nt_index = dict(zip("ACGTN", range(5)))
        nt_counts = [[] for _ind in nt_index]

        for pos in range(N):
            pos_info = mutations[pos]
            pos_ref = pos_info.pop('ref')
            pos_nts = pos_info.items()
            
            total = float(sum(pos_info.values()))
            if total == 0: total = 1.
            for nt, count in pos_nts:
                freq = count/total
                nt_counts[nt_index[nt]].append(freq)
        nt_counts = np.array(nt_counts)

        N_plts = np.round((end-start)/150.)
        for z in range(int(N_plts)):
            fig = plt.figure(figsize=(18.7,10.5))
            ax = fig.add_subplot(111)
            plt.title('Nucleotide Frequencies of Patient %s'%self.patientid)

            l_edge = z*N/N_plts
            r_edge = (z+1)*N/N_plts
            x = np.arange(N)[l_edge:r_edge]
            w = 1
            c = nt_counts[:, l_edge:r_edge]

            plt_a = plt.bar(x, c[0], w, color='g', align='center', alpha=0.5)
            plt_c = plt.bar(x, c[1], w, color='b', align='center', alpha=0.5,
                            bottom=c[:1].sum(0))
            plt_g = plt.bar(x, c[2], w, color='Gold', align='center', alpha=0.5,
                            bottom=c[:2].sum(0))
            plt_t = plt.bar(x, c[3], w, color='r', align='center', alpha=0.5,
                            bottom=c[:3].sum(0))
            plt_n = plt.bar(x, c[4], w, color='k', align='center', alpha=0.5,
                            bottom=c[:4].sum(0))

            plt.legend( (plt_a[0], plt_c[0], plt_g[0], plt_t[0], plt_n[0]), 
                        list('ACGTN'), bbox_to_anchor=(1,1), loc=2)

            plt.ylabel('Nucleotide Frequency')
            plt.xlabel('Gag-Pol sequence position')
            plt.xlim([l_edge -1, r_edge])
            plt.ylim([0, 1])


            locs, labels = plt.xticks()
            plt.xticks(locs, map(lambda x: "%i"%(x+start), locs))
            plt.xlim([l_edge -1, r_edge])

            #mng = plt.get_current_fig_manager()
            #mng.window.maximize()
            fig.savefig(outfile+'%i.png'%z, orientation='landscape')
            #plt.show()
       
       #for nt, ind in nt_index.items():
        #    plt.subplot(511+int(ind))
        #    g = plt.bar(range(N), nt_counts[ind],
        #                linewidth=0, align='center')
        #    lines.append(g)
        #    legend_nts.append(nt)
        #    plt.ylabel('%s Frequency'%nt)
        #    plt.axis([0,50, 0,0.2])
        #plt.xlabel('Gag-Pol Sequence Position')
        #plt.show()
        #fig.savefig('testfig.png', dpi=200)


    def plot_significant_mutations(self, outfile, sig_muts, start, end):
        import matplotlib.pyplot as plt
        x = np.arange(start, end)

        fig = plt.figure(figsize=(18.7,10.5))
        ax1 = fig.add_subplot(111)

        color_table = zip('ACGTN', ('g','b','Gold','r','k'))
        for color_nt, color in color_table:
            y = np.zeros(end-start)
            syn = []
            for mut in sig_muts:
                pos, freq, nt = (mut['pos'], mut['freq'], mut['mut_nt'])
                syn_mut = (mut['ref_aa'] == mut['mut_aa'])
                if nt == color_nt:
                    y[pos] = freq
                    syn.append((pos, syn_mut))

            b = ax1.bar(x, y, width=1, align='center', alpha=0.5, color=color,
                            label=color_nt)
            for i, s in syn:
                if not s:
                    b[i].set_hatch('/')

        ax1.legend(bbox_to_anchor=(1,1), loc=2)
        plt.xlim(start,end)

        plt.title('Significant Nucleotide Mutations of Patient %s'
                    %self.patientid)
        plt.xlabel('Gag-Pol sequence position')
        plt.ylabel('Mutational Frequency')
        fig.savefig(outfile + '.png', orientation='landscape')
#        plt.show()


    def analyze_protease_nucleotides(self):
        counts = self.nucleotide_counts(PRO_NT_REFERENCE, PRO_START, PRO_END)
        sig_muts = self.export_nucleotide_frequencies(counts, 
                                                      self.pro_counts+'_nt', 
                                                      PRO_START)#,threshold=.05)
        #self.plot_nucleotide_frequencies(self.pro_plts, counts, 
        #                                 PRO_START, PRO_END)
        #self.plot_significant_mutations(self.sig_pro_plts, sig_muts, 
        #                                PRO_START, PRO_END)

    def analyze_gag_nucleotides(self):
        counts = self.nucleotide_counts(GAG_NT_REFERENCE, GAG_START, GAG_END)
        sig_muts = self.export_nucleotide_frequencies(counts, 
                                                      self.gag_counts+'_nt',
                                                      GAG_START, threshold=0.05)
        #self.plot_nucleotide_frequencies(self.gag_plts, counts, 
        #                                 GAG_START, GAG_END)
        #self.plot_significant_mutations(self.sig_gag_plts, sig_muts,
        #                                GAG_START, GAG_END)

    def analyze_genome_nucleotides(self):
        self.analyze_gag_nucleotides()
        self.analyze_protease_nucleotides()

    def analyze_protease_amino_acids(self):
        counts = self.protein_counts(PRO_AA_REFERENCE, PRO_START, PRO_END)
        pair_counts = self.protein_pair_counts(PRO_AA_REFERENCE, PRO_START, 
                                               PRO_END)
        self.export_amino_acid_frequencies(counts, self.pro_counts+'_aa', 
                                           PRO_START)
        self.export_amino_acid_pair_frequencies(pair_counts, self.pro_pair_counts+'_aa', 
                                                PRO_START)

    def analyze_protease_amino_acid_pairs(self):
        print "analyzing"
        pair_counts = self.protein_pair_counts(PRO_AA_REFERENCE, PRO_START,
                                               PRO_END)
        self.export_amino_acid_pair_frequencies(pair_counts, self.pro_pair_counts+'_aa',
                                                PRO_START)
        print "done analyzing"


    def analyze_gag_amino_acids(self):
        counts = self.protein_counts(GAG_AA_REFERENCE, GAG_START, GAG_END)
        self.export_amino_acid_frequencies(counts, self.gag_counts+'_aa',
                                           GAG_START)

    def analyze_all(self):
        pro_aa = self.pro_counts+'_aa'
        pro_nt = self.pro_counts+'_nt'
        gag_aa = self.gag_counts+'_aa'
        gag_nt = self.gag_counts+'_nt'
        #pro
        counts = self.nucleotide_counts(PRO_NT_REFERENCE, PRO_START, PRO_END)
        self.export_nucleotide_frequencies(counts, pro_nt, PRO_START)
        counts = self.protein_counts(PRO_AA_REFERENCE, PRO_START, PRO_END)
        self.export_amino_acid_frequencies(counts, pro_aa, PRO_START)
        #gag
        counts = self.nucleotide_counts(GAG_NT_REFERENCE, GAG_START, GAG_END)
        self.export_nucleotide_frequencies(counts, gag_nt, GAG_START)
        counts = self.protein_counts(GAG_AA_REFERENCE, GAG_START, GAG_END)
        self.export_amino_acid_frequencies(counts, gag_aa, GAG_START)

    def analyze_genome_amino_acids(self):
        self.analyze_gag_amino_acids()
        self.anaylze_protease_amino_acids()

    def dNdS_pro(self):
        self.dNdS(PRO_NT_REFERENCE, PRO_START, PRO_END)

    def dNdS_gag(self):
        self.dNdS(GAG_NT_REFERENCE, GAG_START, GAG_END)

def parser_setup():
    import argparse

    descript =  "Process a SAM/BAM file and produce a CSV file containing "\
                "mutations in gag and protease"
    parser = argparse.ArgumentParser(description=descript)
    parser.add_argument('filename', nargs=1)
    parser.add_argument('-g', '--gag', default=False, action='store_true',
                        help="produce file containing gag mutations as well")
    parser.add_argument('--stats', default=False, action='store_true',
                        help="produce sam stats instead of mutation counts")
    parser.add_argument('--cov', default=False, action='store_true',
                        help="return coverage over whole sam file")
    parser.add_argument('--dnds', default=False, action='store_true',
                        help="returns dN/dS ratio over region")
    parser.add_argument('-a', '--aa', default=False, action='store_true',
                        help="produce amino acid counts instead of nucleotide")
    parser.add_argument('--all', default=False, action='store_true',
                        help="produce nt/aa counts for gag and pro")
    parser.add_argument('-p', '--pair', default=False, action='store_true',
                        help="produce only amino acid pair counts")
    parser.add_argument('-n', '--notemp', default=False, action='store_true',
                        help="do not use tempfiles, store .bam files")
    parser.add_argument('-r', '--reads', default=False, action='store_true')

    return parser

if __name__ == "__main__":
    parser = parser_setup()
    args = parser.parse_args()
    analyzer = Analyzer(args.filename[0], use_temp=(not args.notemp))

    if args.stats:
        analyzer.sam_stats()
    elif args.cov:
        print analyzer.sam_coverage()
    elif args.dnds:
        analyzer.dNdS_gag()
    elif args.all:
        analyzer.analyze_all()

    else:
        if args.aa:
            if args.gag:
                analyzer.analyze_gag_amino_acids()
            else:
                if args.reads:
                    analyzer.export_protease_sequences()
                elif args.pair:
                    analyzer.analyze_protease_amino_acid_pairs()
                else:
                    analyzer.analyze_protease_amino_acids()
        else:
            if args.gag:
                analyzer.analyze_gag_nucleotides()
            else:
                analyzer.analyze_protease_nucleotides()
