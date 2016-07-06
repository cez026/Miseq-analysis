#!/usr/bin/env python
import os, sys, re
import argparse
import subprocess

SAMPLE_DIR = '/u2/scripps/samples/'
CONSENSUS_PATH = '/u2/scripps/samples/consensus.fa'

def parser_setup():
    description = "Batch alignment using \"bwa mem\" command"
    parser = argparse.ArgumentParser(description=description)
    #parser.add_argument('directory', metavar='dirname', type=str, nargs=1,
    #                    help='directory in which to find subdirectories')
    #parser.add_argument('consensus', metavar='ref.fa', type=str, nargs=1,
    #                    help='path to consensus.fa')
    return parser

def main():
    """
    This function searches subdirectories in SAMPLE_DIR for pairs of fastq files
    then aligns them with the BWA 'mem' command. The first part of the function
    is highly specific for finding fastq pairs for my dataset/labelling scheme.
    
    The first code block simply returns a list of tuples of the form
        (sample_id, fastq1_path, fastq2_path)

    The second code block iterates through this list and calls the BWA command
    via the subprocess module on each pair of fastq files.

    Aligned .sam files are placed in appropriate subdirectories with the name
        'SAMPLE_DIR/sample_id/sample_id-bwape.sam'
    """
    # Get fastq pairs
    fastq_files = []
    for root, dirs, files in os.walk(SAMPLE_DIR):
        subdir_files = []
        sample_id = root.split(os.sep)[-1]
        for f in files:
            if f.endswith('.fastq.gz') and re.compile('\d').search(sample_id):
                subdir_files.append(os.path.join(root, f))
        if subdir_files == []:
            continue

        id = [sample_id.replace('Sample_', '')]
        if len(subdir_files) == 2:
            fastq_files.append(tuple(id + sorted(subdir_files)))
        elif len(subdir_files) == 4:
            # different primers - normal sorting ok
            if all([f.endswith('001.fasta.gz') for f in subdir_files]):
                subdir_files = sorted(subdir_files)
                fastq_files.append(tuple(id + subdir_files[:2]))
                fastq_files.append(tuple(id + subdir_files[2:]))
            # same primer - need to sort by reverse
            else:
                subdir_files = sorted(subdir_files, key=lambda x: x[::-1])
                fastq_files.append(tuple(id + subdir_files[:2]))
                fastq_files.append(tuple(id + subdir_files[2:]))
    fastq_files = sorted(fastq_files)

    # Run BWA subprocess
    for id, fastq1, fastq2 in fastq_files:
        # Get sam path correctly
        sam_path = os.path.join(SAMPLE_DIR, id)
        sam_path = os.path.join(sam_path, id+'-bwape.sam')

        # Format the read group tag. Just put sample id in everywhere
        rg = '-t 4 -R \"@RG\tID:{0}\tLB:{0}\tPL:illumina\tSM:{0}\"'.format(id)

        # Format subprocess call string.
        call_string = '~/bin/bwa mem {0} {1} {2} {3} > {4}'
        call_string = call_string.format(rg, CONSENSUS_PATH, fastq1, 
                                         fastq2, sam_path)
        subprocess.call(call_string, shell=True)

if __name__ == "__main__":
    # Originally used the parser to pass directory info, but I became lazy
    # Parse is now defunct. 
    parser = parser_setup()
    args = parser.parse_args()

    main()
