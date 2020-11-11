library("data.table")
library("dplyr")
library("stats")

prodigal_transcriptome_blastn <- fread("output/blastn_transcriptome/blastn_transcriptome.outfmt6")
setnames(prodigal_transcriptome_blastn, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new=c("prodigal_nt_id", "transcript_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(prodigal_transcriptome_blastn, prodigal_nt_id, evalue, -bit_score)

##extract result with lowest evalue for each transcript to find which transcripts best match viral genes in genome
transcriptome_blastn <- prodigal_transcriptome_blastn[,.SD[which.min(evalue)], by=prodigal_nt_id]
fwrite(transcriptome_blastn, "output/blastn_transcriptome/blastn_best_hits.csv")
hit_transcripts <- unique(transcriptome_blastn[,list(transcript_id)])
fwrite(hit_transcripts, "output/blastn_transcriptome/transcripts_genome_viral_hits.txt", col.names=FALSE)

##what trinotate hits did those peptides get?
prodigal_genes_transcripts <- transcriptome_blastn[,c(1,2,11)]
trinotate_report <- fread("data/mh_transcriptome/trinotate_annotation_report.txt")
id_vs_annot <- trinotate_report[,c(1,2,3)]
viral_gene_table <- merge(prodigal_genes_transcripts, id_vs_annot, by="transcript_id", all.x=TRUE)
fwrite(viral_gene_table, "output/blastn_transcriptome/transcripts_viral_contigs_trinotate_annots.csv")

##table of transcript & genome blast annot
viral_table_transcripts <- merge(hit_transcripts, viral_gene_table, by="transcript_id", all.x=TRUE)
viral_table_transcripts <- viral_table_transcripts[,c("transcript_id", "sprot_Top_BLASTX_hit")]
fwrite(viral_table_transcripts, "output/blastn_transcriptome/transcripts_annots.csv")

##add back to big viral gene table
prodigal_id_vs_transcriptome <- transcriptome_blastn[,c(1,2,11)]
prodigal_blast_table <- fread("output/prodigal_blastp/blast_gff_coords.csv")
prodigal_blast_table_transcript_ids <- merge(prodigal_blast_table, viral_gene_table, by="prodigal_nt_id", all.x=TRUE)
fwrite(prodigal_blast_table_transcript_ids, "output/blastn_transcriptome/viral_genes_vs_hicgenome_transcriptome.csv")
fwrite(prodigal_blast_table_transcript_ids, "output/viral_gene_table.csv")


