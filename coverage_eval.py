#!/usr/bin/env python

import os
import sys
import argparse
import csv
import subprocess
#import sklearn.cluster as cluster
#import numpy as np

# needs documented and rewritten

class contig:
    def __init__(self, name, length, seq):
        self.depth_position = []
        self.name = name
        self.length = length
        self.windowsize = 49
        self.gccontent = []
        self.seq = seq
        self.label = 1

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

# def calc_clusters(contig_list):
#     [print(c.avg_depth) for c in contig_list]
#     points = np.array([c.avg_depth for c in contig_list])
#     cent,labels,inertia = cluster.k_means(points.reshape(-1, 1), n_clusters=2)
#     [print(l) for l in labels]

#     for i, c in enumerate(contig_list):
#         c.label = labels[i]
    
#     return labels




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculate average coverage across each contig")
    parser.add_argument("bam", type=argparse.FileType('r'))
    parser.add_argument("asm", type=argparse.FileType('r'))
    parser.add_argument("out_dir", type=str, help="output directory to write")
    args = parser.parse_args()
    contigs = {}
    contig_list = []


    if os.path.isdir("out_dir"):
        os.chdir("out_dir")

    for name, seq, qual in readfq(args.asm):
        contigs[name] = contig(name, len(seq), seq.upper())
        contig_list.append(contigs[name])

    depth_cmd = "samtools depth -a {} > {}/depth.out".format(args.bam.name, args.out_dir)
    cp = subprocess.run([ depth_cmd ], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    if cp.returncode:
        print(cp.stdout)
        exit(cp.returncode)
    else:
        with open("{}/depth.out".format(args.out_dir), 'r') as DEPTH:
            for line in DEPTH:
                fields = line.strip().split()
                contigs[fields[0]].depth_position.append(int(fields[2]))

    cvg_file = args.asm.name.replace(".fasta", ".all.cvg")
    with open(cvg_file, 'w') as cvg:
        for name, contig in contigs.items():
            contig.avg_depth = sum(contig.depth_position)/contig.length
            cvg.write("{}\t{}\t{:0.2f}\n".format(name, contig.length, contig.avg_depth))

    window_cmd = "bedtools makewindows -s 50 -g {} -w 100 > {}/windows.bed".format(cvg_file, args.out_dir)
    cp = subprocess.run([ window_cmd ], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    if cp.returncode:
        print(cp.stdout)
        exit(cp.returncode)

    exit(0)

    gc_cmd = "bedtools nuc -fi {} -bed {}/windows.bed  > {}/gcontent.tab".format(args.asm.name, args.out_dir, args.out_dir)
    cp = subprocess.run([ gc_cmd ], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    if cp.returncode:
        print(cp.stdout)
        exit(cp.returncode)
    else:
        with open("{}/gcontent.tab".format(args.out_dir), 'r') as gc:
            for line in gc:
                fields = line.strip().split()
                if len(fields) > 3:
                    if fields[0] != "#1_usercol":
                        contigs[fields[0]].gccontent.append(float(fields[4]))


#    if len(contig_list) > 1:
#        contig_list.sort(key = lambda c: c.avg_depth)
#        calc_clusters(contig_list)


    # keep all contigs that have same label as best contig
#    contig_list.sort(key = lambda c: c.avg_depth, reverse=True)
#    keep_label = contig_list[0].label
    
#    filter_file = "{}/filtered.fasta".format(args.out_dir)
#    filter_cvg =  args.asm.name.replace(".contig", ".cvg")
#    with open(filter_file, 'w') as filter:
#        with open(filter_cvg, 'w') as fcvg:
#            for contig in contig_list:
#                if contig.label == keep_label:
#                    filter.write(">{}\n".format(contig.name))
#                    filter.write("{}\n".format(contig.seq))
#                    fcvg.write("{}\t{}\t{:0.2f}\n".format(name, contig.length, contig.avg_depth))




    with open("{}/gc_cov.tab".format(args.out_dir), 'w') as gc_cov:
        for contig in contig_list:
            if contig.label == keep_label:
                for start_win in range(0, (contig.length - contig.windowsize)):
                    end_win = start_win + contig.windowsize - 1
                    gc_win = contig.seq[ start_win : end_win ].count("C") + contig.seq[ start_win : end_win ].count("G")
                    depth_win = sum( contig.depth_position[start_win:end_win] ) / contig.windowsize
                    gc_cov.write("{}\t{:0.1f}\t{}\n".format(contig.name, depth_win, gc_win * 2))

