#!/usr/bin/python
# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
# snakemake -p --profile <your profile> --keep-going --immediate-submit --use-conda
#

import os
import glob
jn = os.path.join

shell.prefix("source /data/antwerpen/grp/asvardal/share/hscon5_setup.sh; ")

## can also be pulled here from yaml file with
# configfile: "configs/config_main.yaml"

config={'fasta_dir':'PATH_TO_YOUR_DATA',
        'out_dir':'OUTPUT_DIR'  }

rule all:
    input:
       dummy = jn(config['out_dir'], 'dummy')

# file name is hardcoded as "assembly.fasta"
# better to move it to config
rule index_fas:
    input:
        fas = jn(config['fasta_dir'], 'assembly.fasta')
    output:
        faidx = jn(config['fasta_dir'], 'assembly.fasta.fai')
    resources:
        mem_gb=4,
        walltime=1
    shell:
        'samtools faidx {input.fas}'

rule prepare_seq_lists:
    input:
        faidx = jn(config['fasta_dir'], 'assembly.fasta.fai')
    output:
        dynamic(jn(config['out_dir'], 'seq_names_part.{part}'))
    resources:
        mem_gb=4,
        walltime=1
    params:
        seqlist = jn(config['fasta_dir'], 'seq_names.txt'),
        seqlist_parts = jn(config['out_dir'], 'seq_names_part.')
    shell:
        """cut -f1 {input.faidx}  > {params.seqlist} ; 
           split -l 1000 {params.seqlist} {params.seqlist_parts} """

rule split_fas:
    input:
        seqlist_parts = jn(config['out_dir'], 'seq_names_part.{part}'),
        fas = jn(config['fasta_dir'], 'assembly.fasta'),
        faidx = jn(config['fasta_dir'], 'assembly.fasta.fai')
    output:
        fas_parts = jn(config['out_dir'], 'assembly_part.{part}.fas')
    resources:
        mem_gb=4,
        walltime=1
    shell:
       'samtools faidx -r {input.seqlist_parts} {input.fas} > {output.fas_parts}'

rule blast_refseq_euk:
    input:
        fas_parts = jn(config['out_dir'], 'assembly_part.{part}.fas')
    output:
        fas_parts_m6 = jn(config['out_dir'], 'assembly_part.{part}.m6')
    threads: 63
    resources:
        mem_mb=480000,
        mem_gb=480,
        walltime=71,
        # specific to your cluser, better to move to config
        partition="zen3_512"
    shell:
       """module load BLAST+/2.13.0-gompi-2022a;
          blastn \
          -query {input.fas_parts} \
          -db /scratch/antwerpen/grp/asvardal/projects/amphipods/NCBI/ref_euk_rep_genomes \
          -outfmt '6 qseqid staxids bitscore std' \
          -max_target_seqs 1 \
          -max_hsps 1 \
          -num_threads 63 \
          -evalue 1e-25 > {output.fas_parts_m6}"""

rule blast_refseq_prok:
    input:
        fas_parts = jn(config['out_dir'], 'assembly_part.{part}.fas')
    output:
        fas_parts_m6 = jn(config['out_dir'], "refseq_prok", 'assembly_part.{part}.m6')
    threads: 63
    resources:
        mem_mb=50000,
        mem_gb=50,
        walltime=71,
    shell:
       """module load BLAST+/2.13.0-gompi-2022a;
          blastn \
          -query {input.fas_parts} \
          -db /scratch/antwerpen/grp/asvardal/projects/amphipods/NCBI/ref_prok_rep_genomes \
          -outfmt '6 qseqid staxids bitscore std' \
          -max_target_seqs 1 \
          -max_hsps 1 \
          -num_threads 63 \
          -evalue 1e-25 > {output.fas_parts_m6}"""


rule collect:
    input:
       fas_parts_euk = dynamic(jn(config['out_dir'], 'assembly_part.{part}.m6')),
       fas_parts_prok = dynamic(jn(config['out_dir'], "refseq_prok", 'assembly_part.{part}.m6'))
    output:
       dummy = jn(config['out_dir'], 'dummy')
    shell:
       'touch {output.dummy}'