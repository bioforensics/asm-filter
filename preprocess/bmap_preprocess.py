#!/usr/bin/env python

import os
import subprocess
import argparse
from pathlib import Path
import snakemake

# default installation directory
INSTALL_DIR = "./"
config_name = "config.yml"
default_config = INSTALL_DIR + config_name
PREFIX = "bmap_"

# return path to config file
# search sample directory -> cwd -> default
def findConfig(sample_dir):
    config = sample_dir / config_name
    if config.exists():
        return config

    config = Path.cwd() / config_name
    if config.exists():
        return config

    return default_config



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", "-r1", type=argparse.FileType('r'), help="fastq file containing read pair 1")
    parser.add_argument("--r2", "-r2", type=argparse.FileType('r'), help="fastq file containing read pair 2")
    # if not run in main project directory the cd to proj_dir and then run
    parser.add_argument("-s", "--sample", type=str, default=".")
    parser.add_argument("--configfile",  type=argparse.FileType('r'), help="YAML config file")
    parser.add_argument("--cores", type=int, default=32)
    parser.add_argument("--qual", type=int, default=20)
    parser.add_argument("--length", type=int, default=100)
    args = parser.parse_args()

    sample_dir = Path(args.sample)
    seq_dir = sample_dir / "seq"

    config = { "samples" : args.sample }
    if args.qual:
        config["qual"] = args.qual
    if args.length:
        config["minlength"] = args.length

    configfile = ""
    if args.configfile:
        configfile = args.configfile.name
    else:
        configfile = findConfig(sample_dir)

    print(f"config file  is {configfile}")


    if bool(args.r1) ^ bool(args.r2):
        parser.error('--r1 and --r2 both must be given together')
        exit(1)
        
    # if reads are passed make links to then in sample
    if args.r1 and args.r2:

        if not sample_dir.is_dir():
            print(f"{sample_dir} is not already a directory")
            print(" Making sample directory ")
            sample_dir.mkdir(parents=True, exist_ok=True)

        if not seq_dir.is_dir():
            print(f"making seq directories ")
            seq_dir.mkdir(parents=True, exist_ok=True)

        # sample and seq directories should exist at the point
        # now make links into seq dir for the reads
        r1 = Path(args.r1.name)
        r2 = Path(args.r2.name)
        r1_target = (seq_dir / (args.sample + "_R1.fastq.gz"))
        if r1_target.exists():
            print(f" {r1_target} already exists, not linking and using existing file")
        else:
            r1_target.symlink_to(r1.resolve())
                
        r2_target = (seq_dir / (args.sample + "_R2.fastq.gz"))
        if r2_target.exists():
            print(f" {r2_target} already exists, not linking and using existing file")
        else:
            r2_target.symlink_to(r2.resolve())

        workflow = INSTALL_DIR + "preprocess.smk"

        snakemake.snakemake(workflow,  config=config, configfile=configfile, printshellcmds=True, cores=args.cores)
                        #, detailed_summary=True)
                        # workdir=sample_dir,

