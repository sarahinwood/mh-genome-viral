library("data.table")
library("dplyr")
library("rtracklayer")

viral_scaffold_nr_blastp <- fread("output/prodigal_blastp/nr_blastp.outfmt6")
setnames(viral_scaffold_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(viral_scaffold_nr_blastp, peptide_id, evalue, -bit_score)
fwrite(viral_scaffold_nr_blastp, "output/prodigal_blastp/blast_all_res.csv")

##extract result with lowest evalue for each peptide
blast_min_evalues <- viral_scaffold_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]
blast_min_evalues$annotation <- tstrsplit(blast_min_evalues$annotation, "<>", keep=c(1))
fwrite(blast_min_evalues, "output/prodigal_blastp/blast_min_evalues.csv")

##LbFV
##hits before filtering for highest e-value hit
mh_all_leptopilina <- dplyr::filter(viral_scaffold_nr_blastp, grepl('Leptopilina boulardi filamentous virus', annotation))
##no. unique LbFV annots
length(unique(mh_all_leptopilina$annotation))
##highest e-value hit
mh_leptopilina <- dplyr::filter(blast_min_evalues, grepl('Leptopilina boulardi filamentous virus', annotation))
fwrite(mh_leptopilina, "output/prodigal_blastp/lbfv_annots.csv")
##no. unique LbFV annots
length(unique(mh_leptopilina$annotation))

##see which LbFV hits are lost due to lower e-values
lost_lbfv_hits <- data.table(setdiff(mh_all_leptopilina$nr_db_id, mh_leptopilina$nr_db_id))
lost_lbfv_annots <- merge(mh_all_leptopilina, lost_lbfv_hits, by.x="nr_db_id", by.y="V1")
##compare to LbFV hits from prior peptide predictions
pepdb_lbfv <- fread("data/peptide_db_lbfv_annots.csv")
prodigal_only_lbfv <- data.table(setdiff(mh_leptopilina$nr_db_id, pepdb_lbfv$nr_db_id))
prodigal_only_lbfv_annots <- merge(prodigal_only_lbfv, mh_leptopilina, by.x="V1", by.y="nr_db_id")

##read GFF
prodigal_gff <- readGFF("output/prodigal/gene_predictions.gff")
gene_coords <- subset(prodigal_gff, select=c(start, end, ID, conf))
##scaffold id to gene id
scaffold_to_geneid <- data.table(prodigal_gff$seqid, prodigal_gff$ID)
setnames(scaffold_to_geneid, old=c("V1", "V2"), new=c("scaffold_id", "gene_id"))
scaffold_to_geneid$gene_no <- tstrsplit(scaffold_to_geneid$gene_id, "_", keep=c(2))
scaffold_to_geneid$peptide_id <- data.table(paste(scaffold_to_geneid$scaffold_id,"_",scaffold_to_geneid$gene_no, sep=""))
blast_gene_ids<- merge(blast_min_evalues, scaffold_to_geneid, by="peptide_id", all=TRUE)
blast_gff <- merge(blast_gene_ids, gene_coords, by.x="gene_id", by.y="ID", all=TRUE)
fwrite(blast_gff, "output/prodigal_blastp/blast_gff_coords.csv")