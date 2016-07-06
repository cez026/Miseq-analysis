#!/usr/bin/env python
"""
Created on Mar 13, 2013

@author: Bill

NOTE: THIS WHOLE METHODOLOGY SUCKS AND IS WAAYYYY TOO SLOW.
"""

import os
import sys
from sample_analysis import Analyzer

def parser_setup():
    import argparse

    description =   "Batch analysis via sample_analysis.py"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('directory', nargs=1)
    parser.add_argument('samples', nargs='+')
    parser.add_argument('-t', '--threads', default=1, help="Number of threads")
    parser.add_argument('-n', '--notemp', default=False, action='store_true',
                        help="Store bam files")
    return parser

def runtime_func(filename, args):
    analyzer = Analyzer(filename, (not args.notemp))

    analyzer.analyze_protease_amino_acid_pairs()
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
    rootdir = args.directory[0]
    for sampledir in args.samples:
        for root, _, files in os.walk(os.path.join(rootdir, sampledir)):
            for f in files:
                if f.endswith('-realigned.bam'):
                    filenames.append(os.path.join(rootdir, os.path.join(root, f)))

    if int(args.threads) != 1:
        from multiprocessing import Pool
        multiprocessing_run(filenames, args)
    else:
        for filename in filenames:
            runtime_func(filename, args)

if __name__ == "__main__":
    bwa_main()
