library("data.table")
library("dplyr")
library("rtracklayer")

prodigal_gff <- readGFF("output/prodigal/gene_predictions.gff")
prodigal_gff$conf <-as.numeric(prodigal_gff$conf)
setorder(prodigal_gff, -conf)

##complete genes
sum(prodigal_gff$partial=="00")
##incomplete at right edge
sum(prodigal_gff$partial=="10")
##incomplete at left edge
sum(prodigal_gff$partial=="01")
##incomplete both edges
sum(prodigal_gff$partial=="11")


##conf: A confidence score for this gene, representing the probability that this gene is real, i.e. 78.3% means Prodigal believes that gene is real 78.3% of the time and a false positive 21.7% of the time.
##score: The total score for this gene.
##cscore: The hexamer coding portion of the score, i.e. how much this gene looks like a true protein.
##sscore: A score for the translation initiation site for this gene; it is the sum of the following three fields.
##rscore: A score for the RBS motif of this gene.
##uscore: A score for the sequence surrounding the start codon.
##tscore: A score for the start codon type (ATG vs. GTG vs. TTG vs. Nonstandard).
##mscore: A score for the remaining signals (stop codon type and leading/lagging strand information).
