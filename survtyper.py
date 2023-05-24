#!/usr/bin/env python

import argparse, os, time, pysam, pyfaidx
import numpy as np
from random_pos_generator import RandomPositionGenerator

cmd_parser = argparse.ArgumentParser(description='SurVTyper, an SV genotyper.')
cmd_parser.add_argument('vcf_file', help='Input vcf file.')
cmd_parser.add_argument('bam_file', help='Input bam file.')
cmd_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
cmd_parser.add_argument('reference', help='Reference genome in FASTA format.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
cmd_parser.add_argument('--seed', type=int, default=0, help='Seed for random sampling of genomic positions.')
cmd_parser.add_argument('--sampling-regions', help='File in BED format containing a list of regions to be used to estimate'
                                                   'statistics such as depth and insert size distributions.')
cmd_parser.add_argument('--samplename', default='', help='Name of the sample to be used in the VCF output.'
                                                         'If not provided, the basename of the bam/cram file will be used,'
                                                         'up until the first \'.\'')
cmd_parser.add_argument('--max-size-diff', default=100, help='If the event being genotyped exists, but its size'
                                                             'is off by max-size-diff bp, it will be marked as a'
                                                             'FT fail.')
cmd_parser.add_argument('--save-evidence', action='store_true', help='Output evidence for and against each SV.')
cmd_args = cmd_parser.parse_args()

EXEC_PATH = os.path.dirname(os.path.realpath(__file__))

MAX_READS = 1000
GEN_DIST_SIZE = 100000
MAX_ACCEPTABLE_IS = 20000

config_file = open(cmd_args.workdir + "/config.txt", "w")
config_file.write("bam_file %s\n" % cmd_args.bam_file)
config_file.write("workdir %s\n" % cmd_args.workdir)
config_file.write("reference %s\n" % cmd_args.reference)
config_file.write("threads %d\n" % cmd_args.threads)
config_file.write("seed %d\n" % cmd_args.seed)
config_file.write("max_size_diff %d\n" % cmd_args.max_size_diff)
config_file.write("save_evidence %d\n" % cmd_args.save_evidence)

# Find read length
read_len = 0
bam_file = pysam.AlignmentFile(cmd_args.bam_file, reference_filename=cmd_args.reference)
for i, read in enumerate(bam_file.fetch(until_eof=True)):
    if i > MAX_READS: break
    read_len = max(read_len, read.query_length)
config_file.write("read_len %d\n" % read_len)

reference_fa = pyfaidx.Fasta(cmd_args.reference)

rand_pos_gen = RandomPositionGenerator(reference_fa, cmd_args.seed, cmd_args.sampling_regions)
random_positions = []
for i in range(1,1000001):
    if i % 100000 == 0: print(i, "random positions generated.")
    chr, pos = rand_pos_gen.next()
    random_positions.append((chr, pos))
rand_pos_gen = None

with open("%s/random_pos.txt" % cmd_args.workdir, "w") as random_pos_file:
    for random_pos in random_positions:
        random_pos_file.write("%s %d\n" % random_pos)

general_dist = []
rnd_i = 0
while len(general_dist) < GEN_DIST_SIZE:
    chr, pos = random_positions[rnd_i]
    rnd_i += 1

    if pos > len(reference_fa[chr]) - 10000:
        continue

    i = 0
    for read in bam_file.fetch(contig=chr, start=pos, end=pos + 10000):
        if read.is_proper_pair and not read.is_secondary and not read.is_supplementary and \
                0 < read.template_length < MAX_ACCEPTABLE_IS:
            if i > 100:
                break
            i += 1
            general_dist.append(read.template_length)

reference_fa = None

mean_is = np.mean(general_dist)
stddev_is = np.std(general_dist)

general_dist = [x for x in general_dist if abs(x - mean_is) < 5 * stddev_is]

mean_is = int(np.mean(general_dist))
lower_stddev_is = int(np.sqrt(np.mean([(mean_is - x) ** 2 for x in general_dist if x < mean_is])))
higher_stddev_is = int(np.sqrt(np.mean([(x - mean_is) ** 2 for x in general_dist if x > mean_is])))
general_dist = None

min_is, max_is = mean_is - 3 * lower_stddev_is, mean_is + 3.5 * higher_stddev_is

config_file.write("min_is %d\n" % min_is)
config_file.write("avg_is %d\n" % mean_is)
config_file.write("max_is %d\n" % max_is)

if not os.path.exists(cmd_args.workdir + "/mateseqs/"):
    os.makedirs(cmd_args.workdir + "/mateseqs/")
if cmd_args.save_evidence and not os.path.exists(cmd_args.workdir + "/evidence/"):
    os.makedirs(cmd_args.workdir + "/evidence/")

if cmd_args.samplename:
    sample_name = cmd_args.samplename
else:
    sample_name = os.path.basename(cmd_args.bam_file).split(".")[0]

config_file.write("sample_name %s\n" % sample_name)
config_file.close()

start = time.time()
genotype_cmd = EXEC_PATH + "/genotype_svs %s %s" % (cmd_args.vcf_file, cmd_args.workdir)
print("Executing:", genotype_cmd)
os.system(genotype_cmd)
end = time.time()
print("SVs genotyped in %d [s]" % (end-start))
