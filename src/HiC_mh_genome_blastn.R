library("data.table")
library("dplyr")
library("stats")

##Hi-C genome has 12 scaffolds + unscaffolded -> where are viral genes?

prodigal_hic_blastn <- fread("output/blastn_hi_c_genome/blastn_hi_c.outfmt6")
setnames(prodigal_hic_blastn, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new=c("prodigal_nt_id", "hi-c_genome_location", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(prodigal_hic_blastn, prodigal_nt_id, evalue, -bit_score)

##extract result with lowest evalue for each peptide
HiC_blastn <- prodigal_hic_blastn[,.SD[which.min(evalue)], by=prodigal_nt_id]
HiC_vs_old_scaffold_ids <- HiC_blastn[,c(1,2)]

##previous genome viral genes
prodigal_blast_table <- fread("output/prodigal_blastp/blast_gff_coords.csv")
##merge with new genome location
HiC_vs_old_viral_gene_locations <- merge(HiC_vs_old_scaffold_ids, prodigal_blast_table, by.x="prodigal_nt_id", by.y="peptide_id")
fwrite(HiC_vs_old_viral_gene_locations, "output/blastn_hi_c_genome/HiC_vs_oldassembly_viral_genes.csv")
