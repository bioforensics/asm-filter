#!/usr/bin/env python

import os
import sys
import argparse
import csv
import subprocess
import tempfile

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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="fasta file to check for circles")
    parser.add_argument("asm", type=argparse.FileType('r'))
    parser.add_argument("out_dir", type=str)
    args = parser.parse_args()
    
    if not os.path.isdir(args.out_dir):
        print(f"{args.out_dir} is not a directory")
        exit(1)

    # minimum size of sequence to check for circle
    MINSIZE = 500
    WINDOW = 15

    with open(f"{args.out_dir}/{args.asm.name}_circle.fasta", 'w') as CIRCLE:

        for name, seq, qual in readfq(args.asm):
            
            # seq is over 1k split in half and align to self
            # split seq 
            if len(seq) < MINSIZE:
                continue
            first, second = seq[:len(seq)//2], seq[len(seq)//2:]
            node = "_".join(name.split("_")[:2])
            with tempfile.NamedTemporaryFile('w') as front, tempfile.NamedTemporaryFile('w') as back:
                front.write(f">{node}_front\n")
                front.write(first)
                front.flush()
                os.fsync(front.fileno())
                back.write(f">{node}_back\n")
                back.write(second)
                back.flush()
                os.fsync(back.fileno())

                align_cmd = f"nucmer --prefix {node} -f --nosimplify -maxmatch -l 10 -c 15 {front.name} {back.name} && \
                          show-coords -c -l -o -r -T -H -I 95 {node}.delta > {args.out_dir}/{node}.coords"
                cp = subprocess.run([ align_cmd ], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
                os.remove(f"{node}.delta")
                # if return is good then open coords file
                circlep = False
                with open(f"{args.out_dir}/{node}.coords", 'r') as COORDS:
                    aligns = csv.reader(COORDS, delimiter='\t')
                    for align in aligns:
                        if int(align[0]) < WINDOW:
                            if int(align[8]) - int(align[3]) < WINDOW:
                                remove = int(align[8]) - int(align[2])
                                CIRCLE.write(f">{name}_circle_length_{len(seq[:-remove])}\n")
                                CIRCLE.write(seq[:-remove] + "\n")
                                CIRCLE.flush()
                                os.fsync(CIRCLE.fileno())
                                circlep = True

                if circlep == False:            
                    os.remove(f"{args.out_dir}/{node}.coords")
            

