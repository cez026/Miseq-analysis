#!/usr/bin/env python
"""
Created on Mar 13, 2013

@author: Bill
"""

import os
import sys
from sample_analysis import Analyzer

def parser_setup():
    import argparse

    description =   "Batch analysis via sample_analysis.py"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('directory', nargs=1)
    parser.add_argument('-g', '--gag', default=False, action='store_true',
                        help="produce files containing gag mutations as well")
    parser.add_argument('-a', '--aa', default=False, action='store_true',
                        help="produce amino acid counts instead of nucleotide")
    parser.add_argument('-t', '--threads', default=1, help="Number of threads")
    parser.add_argument('--stats', default=False, action='store_true',
                        help="produce sam stats instead of mutation counts")
    parser.add_argument('-n', '--notemp', default=False, action='store_true',
                        help="Store bam files")
    parser.add_argument('--all', default=False, action='store_true',
                        help="produce aa/nt counts for gag and pro")
    return parser

def runtime_func(filename, args):
    analyzer = Analyzer(filename, (not args.notemp))

    if args.stats:
        analyzer.stats()
    elif args.all:
        analyzer.analyze_all()
    else:
        if args.aa:
            if args.gag:
                analyzer.analyze_gag_amino_acids()
            else:
                analyzer.analyze_protease_amino_acids()
        else:
            if args.gag:
                analyzer.analyze_gag_nucleotides()
            else:
                analyzer.analyze_protease_nucleotides()
    del(analyzer)

def multiprocessing_run(filenames, args):
    from multiprocessing import Pool
    pool = Pool(processes=int(args.threads))
    for filename in filenames:
        pool.apply_async(runtime_func, args=(filename, args))
    pool.close()
    pool.join()


def bwa_main():
    parser = parser_setup()
    args = parser.parse_args()
    filenames = []
    for root, _, files in os.walk(args.directory[0]):
        for f in files:
            if f.endswith('-realigned.bam'):
                filenames.append(os.path.join(root, f))

    if int(args.threads) != 1:
        from multiprocessing import Pool
        multiprocessing_run(filenames, args)
    else:
        for filename in filenames:
            runtime_func(filename, args)

def old_main():
    parser = parser_setup()
    args = parser.parse_args()

    filenames = []
    #for root, dirs, files in os.walk(args.directory[0]):
    #    for f in files:
    #        if not f.endswith('-refined.sam'):
    #            continue
    #        filenames.append(os.path.join(root, f))
    for root, dirs, files in os.walk(args.directory[0]):
        have_sam, have_bam = False, False
        last_resort = ""
        for f in files:
            if not f.endswith('-refined.sam') and not f.endswith('-refined.bam_sorted.bam'):
                continue
            if f.endswith('.bam'):
                filenames.append(os.path.join(root, f))
                have_bam = True
                break
            else:
                have_sam = True
                last_resort = os.path.join(root, f)
        if have_sam and not have_bam:
            filenames.append(last_resort)

    filenames = sorted(filenames)

    if int(args.threads) != 1:
        from multiprocessing import Pool
        multiprocessing_run(filenames, args)
    else:
        for filename in filenames:
            runtime_func(filename, args)

if __name__ == "__main__":
    bwa_main()
