#!/usr/bin/env python

import os
import subprocess
import argparse
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fast1", type=argparse.FileType('r'))
    parser.add_argument("fast2", type=argparse.FileType('r'))
    parser.add_argument("sketch", type=str)
    
    parser.add_argument("--genome_size", type=str)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min", type=int, default=2)
#    parser.add_argument("--jellyfish", type=bool, default=False)
#  jellyfish no longer supported
    
    args = parser.parse_args()

    mash_command = f"mash sketch -m {args.min} -p {args.threads} {args.fast1.name} -o {args.sketch}"
    cp = subprocess.run(mash_command, shell=True, stderr=subprocess.PIPE, universal_newlines=True)
    if cp.returncode:
        print(cp.stderr)
        exit(cp.returncode)

    s = re.search(r'Estimated genome size: ([e\d\.\+]+)', cp.stderr)
    if s:
        genome_size = int(float(s.group(1)))
    else:
        print("genome size not found")
        exit(1)
            
    if args.genome_size:
        print(f"{args.genome_size}")
    else:
        print(f"{genome_size}")
