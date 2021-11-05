#!/usr/bin/env python3

import pathlib2
import pandas
import os

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

def find_read_files(read_dir):
#Make list of files
    path_generator = os.walk(read_dir, followlinks = True)
    my_files = list((dirpath, filenames)
        for (dirpath, dirname, filenames)
        in path_generator)
#Make new dictionary & populate with files (flowcell = key)
    my_fastq_files = {}
    for dirpath, filenames in my_files:
        for filename in filenames:
            if filename.endswith('.fastq.gz'):
                my_flowcell = pathlib2.Path(dirpath).name
                my_fastq = str(pathlib2.Path(dirpath,filename))
                if my_flowcell in my_fastq_files:
                    my_fastq_files[my_flowcell].append(my_fastq)
                else:
                    my_fastq_files[my_flowcell]= []
                    my_fastq_files[my_flowcell].append(my_fastq)
    return(my_fastq_files)

def sample_name_to_fastq(wildcards):
    sample_row = sample_key[sample_key['Sample_name'] == wildcards.sample]
    sample_id = sample_row.iloc[-1]['OGF_sample_ID']
    sample_flowcell = sample_row.iloc[-1]['Flow_cell']
    sample_all_fastq = [x for x in all_fastq[sample_flowcell]
                        if '-{}-'.format(sample_id) in x]
    sample_r1 = sorted(list(x for x in sample_all_fastq
                            if '_R1_' in os.path.basename(x)))
    sample_r2 = sorted(list(x for x in sample_all_fastq
                            if '_R2_' in os.path.basename(x)))
    return({'r1': sample_r1, 'r2': sample_r2})

###########
# GLOBALS #
###########

star_reference_folder = 'output/star/star_reference'

read_dir = 'data/reads'
sample_key_file = 'data/reads/rnaseq_sample_key.csv'
sample_key = pandas.read_csv(sample_key_file)

bbduk_adapters = '/adapters.fa'
bbduk_ref = '/phix174_ill.ref.fa.gz'

##############
# CONTAINERS #
##############

bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
star_container = 'shub://TomHarrop/singularity-containers:star_2.7.0c'
tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'
interproscan_container = 'docker://agbase/interproscan:5.45-80.0_0'

#########
# SETUP #
#########

# generate name to filename dictionary
all_fastq = find_read_files(read_dir)
all_samples = sorted(set(sample_key['Sample_name']))

#########
# RULES #
#########

rule target:
    input:
        expand('output/star/star_pass2/{sample}.Aligned.sortedByCoord.out.bam.bai', sample=all_samples),
        'output/blastn_transcriptome/blastn_transcriptome.outfmt6'

######################################################
## Blast viral scaffold genes against transcriptome ##
######################################################

##look for transcripts in transcriptome with hits to viral peptides
rule blastn_prodigal_preds_mh_transcriptome:
    input:
        prodigal_nt_preds = 'output/prodigal/nucleotide_seq.fasta',
        blast_db = 'output/blastdb/mh_genome.nhr'
    output:
        blastn_res = 'output/blastn_transcriptome/blastn_transcriptome.outfmt6'
    params:
        hyperodae_transcriptome_db = 'output/blastdb/mh_transcriptome'
    threads:
        10
    log:
        'output/logs/blastn_transcriptome_prodigal_preds.log'
    shell:
        'blastn '
        '-query {input.prodigal_nt_preds} '
        '-db {params.hyperodae_transcriptome_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastn_res} '
        '2>{log}'

rule mh_transcriptome_blast_db:
    input:
        mh_transcriptome = 'data/mh_transcriptome/mh_transcriptome_length_filtered.fasta',
        blast_db = 'output/blastdb/mh_genome.nhr'
    output:
        blast_db = 'output/blastdb/mh_transcriptome.nhr'
    params:
        db_name = 'mh_transcriptome',
        db_dir = 'output/blastdb/mh_transcriptome'
    threads:
        10
    log:
        'output/logs/hyperodae_transcriptome_blast_db.log'
    shell:
        'makeblastdb '
        '-in {input.mh_transcriptome} '
        '-dbtype nucl '
        '-title {params.db_name} '
        '-out {params.db_dir} '
        '-parse_seqids '
        '2> {log}'

rule hyperodae_transcriptome_blast_db:
    input:
        mh_genome = 'data/Mh_assembly.fa'
    output:
        blast_db = 'output/blastdb/mh_genome.nhr'
    params:
        db_name = 'mh_genome',
        db_dir = 'output/blastdb/mh_genome'
    threads:
        20
    log:
        'output/logs/hyperodae_blast_db.log'
    shell:
        'makeblastdb '
        '-in {input.mh_genome} '
        '-dbtype nucl '
        '-title {params.db_name} '
        '-out {params.db_dir} '
        '-parse_seqids '
        '2> {log}'     

##########################
## mapping RNAseq reads ##
##########################

##index so can view in igv
rule index_star_bam:
    input:
        bam = 'output/star/star_pass2/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        bai = 'output/star/star_pass2/{sample}.Aligned.sortedByCoord.out.bam.bai'
    threads:
        20
    singularity:
        samtools_container
    log:
        'output/logs/{sample}_index.log'
    shell:
        'samtools '
        'index '
        '{input.bam} '

###map RNAseq onto genome & see whether viral genes are expressed, in what tissues & whether reads map over introns
rule star_second_pass:
    input:
        left = 'output/bbduk_trim_rna/{sample}_r1.fq.gz',
        right = 'output/bbduk_trim_rna/{sample}_r2.fq.gz',
        junctions = expand('output/star/star_pass1/{sample}.SJ.out.tab', sample=all_samples)
    output:
        bam = 'output/star/star_pass2/{sample}.Aligned.sortedByCoord.out.bam'
    threads:
        30
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/star/star_pass2/{sample}.'
    log:
        'output/logs/star/star_pass2_{sample}.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--sjdbFileChrStartEnd {input.junctions} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outBAMcompression 10 '
        '--readFilesIn {input.left} {input.right} '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_first_pass:
    input:
        left = 'output/bbduk_trim_rna/{sample}_r1.fq.gz',
        right = 'output/bbduk_trim_rna/{sample}_r2.fq.gz',
        star_reference = 'output/star/star_reference/Genome'
    output:
        sjdb = 'output/star/star_pass1/{sample}.SJ.out.tab'
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/star/star_pass1/{sample}.'
    threads:
        30
    log:
        'output/logs/star/star_pass1_{sample}.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--outSAMtype None '
        '--readFilesIn {input.left} {input.right} '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_reference:
    input:
        mh_genome = 'data/hi-c_genome/Mh_Hi-C_PGA_assembly.fasta'
    output:
        'output/star/star_reference/Genome'
    params:
        genome_dir = star_reference_folder
    threads:
        20
    log:
        'output/logs/star/star_reference.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.genome_dir} '
        '--genomeFastaFiles {input.mh_genome} '
        '&> {log} '

rule bbduk_trim_rna:
    input:
        r1 = 'output/joined/{sample}_r1.fq.gz',
        r2 = 'output/joined/{sample}_r2.fq.gz'
    output:
        r1 = 'output/bbduk_trim_rna/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim_rna/{sample}_r2.fq.gz'
    params:
        adapters = bbduk_adapters
    log:
        'output/logs/bbduk_trim/{sample}.log'
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'

##won't work on ovary samples - have had to do this manually myself for now
rule cat_reads:
    input:
        unpack(sample_name_to_fastq)
    output: 
        r1 = temp('output/joined/{sample}_r1.fq.gz'),
        r2 = temp('output/joined/{sample}_r2.fq.gz')
    threads:
        1
    shell:
        'cat {input.r1} > {output.r1} & '
        'cat {input.r2} > {output.r2} & '
        'wait'
