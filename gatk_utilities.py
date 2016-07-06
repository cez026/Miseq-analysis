import os, sys
import argparse
import subprocess
from time import time

SAMPLE_DIR = '/u2/scripps/samples'
GATK_PATH = '~/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar'
CONSENSUS_PATH = os.path.join(SAMPLE_DIR, 'consensus.fa')
TAQ_PATH = os.path.join(SAMPLE_DIR, 'taq.grp')

def parser_setup():
    description = "Performs GATK BQSR and indel realignment on all samples"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-r', '--recal', default=False, action='store_true',
                        help='recalibrate and realign .bam files')
    return parser

def generate_ids():
    ids = []
    for root, dirs, files in os.walk(SAMPLE_DIR):
        id = root.split(os.sep)[-1]
        if id is '' or not id[0].isdigit():
            continue
        ids.append(id)
    return sorted(ids)

def debug(message):
    print '\t%s'%message

def recalibrate_bam():
    ids = generate_ids()
    devnull = open('/dev/null', 'w')
    for id in ids:
        start_time = time()
        print 'recalibrating: %s'%id

        base_path = os.path.join(os.path.join(SAMPLE_DIR, id), id+'-bwape')
        bam_path = base_path + '.bam'
        bqsr_path = base_path + '.BQSR.bam'
        inter_path = base_path + '.intervals'
        recal_path = base_path + '.BQSR-realigned.bam'
        bqsr_ind_path = base_path + '.BQSR.bai'
        recal_ind_path = base_path + '.BQSR-realigned.bai'

        bqsr_call = 'java -jar {0} -T PrintReads -R {1} -I {2} -BQSR {3} -o {4}'
        bqsr_call = bqsr_call.format(GATK_PATH, CONSENSUS_PATH, bam_path, 
                                     TAQ_PATH, bqsr_path)

        inter_call = 'java -Xmx4g -jar {0} -T RealignerTargetCreator -R {1} '\
                     '-I {2} -o {3}'
        inter_call = inter_call.format(GATK_PATH, CONSENSUS_PATH, bqsr_path,
                                       inter_path)

        recal_call = 'java -Xmx4g -jar {0} -T IndelRealigner -R {1} -I {2} '\
                     '--maxReadsForRealignment 500000 -targetIntervals {3} '\
                     '-o {4}'
        recal_call = recal_call.format(GATK_PATH, CONSENSUS_PATH, bqsr_path,
                                       inter_path, recal_path)

        rm_call = 'rm {0} {1} {2}'.format(bqsr_path, bqsr_ind_path, inter_path)

        subprocess.call(bqsr_call, shell=True, stdout=devnull)
        debug('bqsr complete')
        subprocess.call(inter_call, shell=True, stdout=devnull)
        debug('interval calling complete')
        subprocess.call(recal_call, shell=True, stdout=devnull)
        debug('recalibration complete')
        subprocess.call(rm_call, shell=True, stdout=devnull)


if __name__ == "__main__":
    parser = parser_setup()
    args = parser.parse_args()

    if args.recal:
        recalibrate_bam()
