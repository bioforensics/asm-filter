# Assembly Tools
 Snakemake workflows used to assemble bacterial isolates sequenced on illumina paired end sequencers.

These workflows have been used to assemble five historical Bacillus anthracis assemblies soon to be published in Microbiology Resource Annoucements (http://mra.asm.org).  

The assemblies have been deposited in DDBJ/ENA/GenBank under the accession numbers (SAMN12620928, SAMN12620929, SAMN12620930, SAMN12620931, SAMN12620932).  The raw Illumina paired end sequencing reads are archived in the SRA under the accession numbers (SRR10019497, SRR10019498, SRR10019499, SRR10019500, SRR10019501).

# Preprocessing Pair End Reads

The preprocessing.smk snakemake workflow is prepares illumina reads to be assembled.
1. Run fastp to trim adapter sequence, low quality bases, and very short reads.  By default bases below Q20 at ends of reads will be trimmed. Any reads below length 75 and/or containing Ns will be removed.  
2. Run "mash screen" against Refseq to check for contamenent.
3. Estimate Genome size by building a k-mer profile on the reads.
4. Randomly downsample reads to 150x coverage of the estimated genome size using sample-reads program.

# Snakemake GraphViz Image Here
