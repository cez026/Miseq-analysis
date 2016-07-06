import os, sys, re
import pysam
import argparse
import subprocess
from time import time

SAMPLE_DIR = '/u2/scripps/samples/'

def parser_setup():
    description = "Various batch samtools .bam operations"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-c', '--convert', default=False, action='store_true',
                        help='convert .sam to .bam')
    parser.add_argument('-s', '--sort', default=False, action='store_true',
                        help='sort .bam')
    return parser


def generate_ids():
    ids = []
    for root, dirs, files in os.walk(SAMPLE_DIR):
        id = root.split(os.sep)[-1]
        if id is '' or not id[0].isdigit():
            continue
        ids.append(id)
    return sorted(ids)


def convert_sam_to_bam():
    """
    This method should take a newly create .sam file from alignment and
        - convert it to .bam
        - sort .bam
        - index .bam
    """
    ids = generate_ids()
    for id in ids:
        start_time = time()
        print 'converting: %s'%id
        base_path = os.path.join(SAMPLE_DIR, id)
        sam_path = os.path.join(base_path, id+'-bwape.sam')
        bam_path = os.path.join(base_path, id+'-bwape.bam')

        bam_content = pysam.view('-bS', sam_path)
        bam_file = open(bam_path, 'w+')
        bam_file.writelines(bam_content)
        bam_file.close()

        pysam.sort(bam_path, bam_path+'_sorted')
        pysam.index(bam_path+'_sorted.bam')

        # indexing creates file.bam.bam. Move it to file.bam
        bam_call = "mv {0} {1}".format(bam_path+'_sorted.bam', bam_path)
        index_call = "mv {0} {1}".format(bam_path+'_sorted.bam.bai',
                                         bam_path+'.bam.bai')
        subprocess.call(bam_call, shell=True)
        subprocess.call(index_call, shell=True)
        end_time = time()
        print 'completed: %.3fs'%(end_time-start_time)


if __name__ == "__main__":
    parser = parser_setup()
    args = parser.parse_args()

    if args.convert:
        convert_sam_to_bam()
