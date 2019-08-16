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
sample_key_file = 'data/rnaseq_sample_key.csv'
sample_key = pandas.read_csv(sample_key_file)

bbduk_adapters = '/adapters.fa'
bbduk_ref = '/phix174_ill.ref.fa.gz'

##############
# CONTAINERS #
##############

bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
star_container = 'shub://TomHarrop/singularity-containers:star_2.7.0c'
tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'

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
    	##do rnaseq reads map onto viral scaffolds & if so do they map over introns?
        expand('output/star/star_pass2/{sample}.Aligned.sortedByCoord.out.bam.bai', sample=all_samples),
        'output/reapr/pipeline/05.summary.report.txt',
        'output/samtools/mean_depth_table.csv',
        'output/bb_stats/gc.txt'

rule bb_stats:
    input:
        genome = 'data/Mh_assembly.fa'
    output:
        gc = 'output/bb_stats/gc.txt',
        stats = 'output/bb_stats/bb_stats.out',
        gc_hist = 'output/bb_stats/gc_hist.out'
    log:
    singularity:
        bbduk_container
    shell:
        'stats.sh '
        'in={input.genome} '
        'out={output.stats} '
        'gc={output.gc} '
        'gcformat=4 '
        'gchist={output.gc_hist} '

##end reapr smalt bam depth

rule samtools_depth_smalt:
    input:
        smalt_bam = 'output/reapr/smalt/smalt.bam'
    output:
        depth_out = 'output/samtools/smalt_depth.out'
    log:
        'output/logs/samtools_depth_smalt.log'
    threads:
        20
    shell:
        'samtools depth '
        '{input.smalt_bam} '
        '-o {output.depth_out} '
        '-H '
        '2> {log}'

###try for reapr smalt bam output too

rule calc_mean_depth:
    input:  
        depth = 'output/samtools/depth.out'
    output:
        mean_depth_table = 'output/samtools/mean_depth_table.csv'
    singularity:
        tidyverse_container
    threads:
        20
    log:
        'output/logs/calc_mean_depth.log'
    script:
        'src/calc_mean_depth.R'


rule samtools_depth:
    input:
        sorted_bam = 'output/samtools/sorted.bam'
    output:
        depth_out = 'output/samtools/depth.out'
    log:
        'output/logs/samtools_depth.log'
    threads:
        20
    shell:
        'samtools depth '
        '{input.sorted_bam} '
        '-o {output.depth_out} '
        '-H '
        '2> {log}'

rule samtools_sort:
    input:
        sam = 'output/bwa/bwa_mem.sam'
    output:
        sorted_bam = 'output/samtools/sorted.bam'
    log:
        'output/logs/samtools_sort.log'
    threads:
        20
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output.sorted_bam} '
        '2> {log}'

rule bwa_mem:
    input:
        index = 'output/bwa/index/index.bwt',
        trimr1 = 'output/bbduk_trim_dna/mhyp_trimr1.fq.gz',
        trimr2 = 'output/bbduk_trim_dna/mhyp_trimr2.fq.gz'
    output:
        sam = 'output/bwa/bwa_mem.sam'
    params:
        index_dir = 'output/bwa/index/index'
    threads:
        50
    log:
        'output/logs/bwa_mem.log'
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.index_dir} '
        '{input.trimr1} '
        '{input.trimr2} '
        '> {output.sam} '
        '2> {log}'

##bwa-mem for coverage of scaffolds - compare to reapr smalt?
rule bwa_index:
    input:
        genome = 'data/Mh_assembly.fa'
    output:
        index = 'output/bwa/index/index.bwt'
    params:
        outdir = 'output/bwa/index/index'
    threads:
        20
    log:
        'output/logs/bwa_index.log'
    shell:
        'bwa index '
        '{input.genome} '
        '-p {params.outdir} '
        '2> {log} '

##run reapr - does this give the same output as that in the LbFV genome paper?
rule reapr_pipeline:
    input:
        mh_genome = 'data/Mh_assembly.fa',
        bam = 'output/reapr/smalt/smalt.bam',
        perfectmap = 'output/reapr/perfectmap/perfect.perfect_cov.gz'
    output:
        summary = 'output/reapr/pipeline/05.summary.report.txt'
    params:
        perfectmap_prefix = 'output/reapr/perfectmap/perfect',
        outdir = 'output/reapr/pipeline'
    threads:
        50
    log:
        'output/logs/reapr/pipeline.log'
    shell:
        'bin/Reapr_1.0.18/reapr pipeline '
        '{input.mh_genome} '
        '{input.bam} '
        '{params.outdir} '
        '{params.perfectmap_prefix} '
        '2> {log} '

rule reapr_perfectmap:
    input:
        mh_genome = 'data/Mh_assembly.fa',
        ##perfectmap needs full path for symlinks or fails & also needs all reads same length so not adapter/quality trimmed
        filr1 = '/Volumes/userdata/student_users/sarahinwood/Projects/mh-projects/mh-genome/output/bbduk_trim_dna/mhyp_filr1.fq.gz',
        filr2 = '/Volumes/userdata/student_users/sarahinwood/Projects/mh-projects/mh-genome/output/bbduk_trim_dna/mhyp_filr2.fq.gz'
    output:
        coverage_out = 'output/reapr/perfectmap/out.perfect_cov.gz'
    params:
        prefix = 'output/reapr/perfectmap/perfect'
    threads:
        50
    log:
        'output/logs/reapr/perfectmap.log'
    shell:
        'bin/Reapr_1.0.18/reapr perfectmap '
        '{input.mh_genome} '
        '{input.filr1} '
        '{input.filr2} '
        '101 '
        '{params.prefix} '
        '2> {log} '

rule reapr_smalt:
    input:
        mh_genome = 'data/Mh_assembly.fa',
        trimr1 = 'output/bbduk_trim_dna/mhyp_trimr1.fq',
        trimr2 = 'output/bbduk_trim_dna/mhyp_trimr2.fq'
    output:
        bam = 'output/reapr/smalt/smalt.bam'
    threads:
        50
    log:
        'output/logs/reapr/smalt.log'
    shell:
        'bin/Reapr_1.0.18/reapr smaltmap '
        '{input.mh_genome} '
        '{input.trimr1} '
        '{input.trimr2} '
        '{output.bam} '
        '-n {threads} '
        '2> {log} '

##temp gunzip for smalt
rule unzip_bbduk:
    input:
        trimr1 = 'output/bbduk_trim_dna/mhyp_trimr1.fq.gz',
        trimr2 = 'output/bbduk_trim_dna/mhyp_trimr2.fq.gz'
    output:
        unzipr1 = temp('output/bbduk_trim_dna/uz_mhyp_trimr1.fq'),
        unzipr2 = temp('output/bbduk_trim_dna/uz_mhyp_trimr2.fq')
    threads:
        10
    shell:
        'gunzip {input.trimr1} > {output.unzipr1} & '
        'gunzip {input.trimr2} > {output.unzipr2} & '
        'wait '

##trim and decontaminate DNA reads to map onto genome
rule bbduk_trim_dna:
    input:
        filr1 = 'output/bbduk_trim_dna/mhyp_filr1.fq.gz',
        filr2 = 'output/bbduk_trim_dna/mhyp_filr2.fq.gz'
    output:
        trimr1 = 'output/bbduk_trim_dna/mhyp_trimr1.fq.gz',
        trimr2 = 'output/bbduk_trim_dna/mhyp_trimr2.fq.gz',
        t_stats = 'output/bbduk_trim_dna/trim-stats.txt'
    log:
        trim = 'output/logs/bbduk_trim_dna/trim.log'
    params:
        adapters = bbduk_adapters
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in1={input.filr1} '
        'in2={input.filr2} '
        'int=f '
        'out1={output.trimr1} '
        'out2={output.trimr2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '

rule bbduk_filter_dna:  
    input:
        r1 = 'data/1704KHP_macrogen_M_hyperodae/Mhyp_1.fastq.gz',
        r2 = 'data/1704KHP_macrogen_M_hyperodae/Mhyp_2.fastq.gz'
    output:
        filr1 = 'output/bbduk_trim_dna/mhyp_filr1.fq.gz',
        filr2 = 'output/bbduk_trim_dna/mhyp_filr2.fq.gz',
        f_stats = 'output/bbduk_trim_dna/filter-stats.txt'
    log:
        filter = 'output/logs/bbduk_trim_dna/filter.log'
    params:
        ref = bbduk_ref
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in1={input.r1} '
        'in2={input.r2} '
        'out1={output.filr1} '
        'out2={output.filr2} '
        'ref={params.ref} '
        'hdist=1 '
        'stats={output.f_stats} '       
        '2> {log.filter} '         

##########working with RNAseq reads from here down##########

##index so can view in igv
rule index_star_bam:
	input:
		bam = 'output/star/star_pass2/{sample}.Aligned.sortedByCoord.out.bam'
	output:
		bai = 'output/star/star_pass2/{sample}.Aligned.sortedByCoord.out.bam.bai'
	threads:
		20
	log:
		'output/logs/{sample}_index.log'
	shell:
		'samtools '
		'index '
		'{input.bam} '

###map RNAseq onto genome & see whether viral genes are expressed, in what tissues & whether reads map over introns?
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
        mh_genome = 'data/Mh_assembly.fa'
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
