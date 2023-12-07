import os
import argparse
import gzip
import logging
import subprocess

#this script requires individual VCFs
parser = argparse.ArgumentParser()
parser.add_argument('-vd', '--VCF_directory', required=True, type=str,help='path to directory of single-sample VCFs')
parser.add_argument('-d', '--working_directory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-smf', '--species_maskfile', required=False, type=str, help='path to bed file of commonly masked regions of genome')
parser.add_argument('-bed', '--bedgraph', required=False, type=str, help="path to bed coverage file (bedgraph) for vcf")
parser.add_argument('-cd', '--coverage_depth', required=False, default=10, type=int, help="minimum coverage depth for any given call before that call is considered dubious")
parser.add_argument('-l', '--logging', required=False, default=True, type=bool, help="if True, logging.debug verbose logging to diff.log, else suppress most logging")

args = parser.parse_args()
vd = args.VCF_directory
wd = args.working_directory
#smf = args.species_maskfile
#bed = args.bedgraph
#min_coverage = args.coverage_depth

for v in os.listdir(vd):
    print(v)
    arg = f'python3 ~/scripts/parsevcf/vcf_to_diff_script.py -v {vd}{v} -d {wd}'
    #arg = arg.split(' ')
    print('arg', arg)
    subprocess.run(arg, shell=True, check=True)