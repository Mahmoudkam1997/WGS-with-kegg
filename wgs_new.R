rm(list=ls())
library(ShortRead)       # For handling FASTQ files
library(Biostrings)      # For sequence operations
library(Rsubread)        # For alignment
library(VariantAnnotation) # For variant calling
library(GenomicAlignments) # For BAM file manipulation
library(tidyverse)
library(AnnotationHub)
###############################################################
#################################################################
fastq_wgs <- "E:/wgs_tut/"
reads <- readFastq(fastq_wgs)
############################################################
qu <- quality(reads)
summary(qu)
####################################################
qa <- qa(fastq_wgs, type = "fastq")
report(qa, dest="quality_report")
#########################################################
trimmed_reads <- trimTails(reads, 5, "1", consecutive = TRUE)
filtered_reads <- reads[width(reads) > 50]
####################################################################
writeFastq(filtered_reads, "filtered_reads.fastq.gz", compress=TRUE)
#####################################################################
reference_genome <- "E:/GCF_000005845.2_ASM584v2_genomic.fna"
buildindex(basename = "ref_index", reference = reference_genome)
######################################################################
# Align reads
align(index = "ref_index",
      readfile1 = "filtered_reads.fastq.gz",
      output_file = "aligned_reads.bam",
      nthreads = 4)
##############################################################
library(Rsamtools)

# Sort BAM file
sorted_bam <- "E:/aligned_reads.sorted"
sortBam("E:/aligned_reads.bam", sorted_bam)

# Index BAM file
indexBam(sorted_bam)
################################################################
#variant calling in Galaxy by bcf tools
############################################################
#variant annotation
vcf <- readVcf("E:/Galaxy5-[bcftools mpileup on data 3 and data 4].vcf", 
               genome = "E:/GCF_000005845.2_ASM584v2_genomic.fna")
info(vcf)
#################################################################
#Annotate variants another way 
# Load the genome annotation
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("E:/GCF_000005845.2_ASM584v2_genomic.gff")

# Extract genes
genes <- genes(txdb)

# Find overlaps between variants and genes
overlaps <- findOverlaps(rowRanges(vcf), genes)

# View annotated variants
annotated_variants <- genes[subjectHits(overlaps)]
annotated_variants
##############################################################
BiocManager::install("Gviz")
library(Gviz)
# Load BAM file
bamTrack <- AlignmentsTrack("E:/aligned_reads.sorted.bam", 
                            genome = "txdb")

# Create a genome axis track
genomeAxis <- GenomeAxisTrack()

# Plot the tracks
plotTracks(list(genomeAxis, bamTrack), from = 100000, to = 101000)
#####################################################################
library(clusterProfiler)
BiocManager::install("org.EcK12.eg.db")  # For *E. coli*
# Convert gene IDs to KEGG
library(org.EcK12.eg.db)
kegg_genes <- bitr_kegg(geneID = annotated_variants$gene_id,
                        fromType = "ncbi-geneid",  # KEGG expects ENTREZ IDs
                        toType = "Path",
                        organism = "eco")    

# Get gene details
genes <- select(org.EcK12.eg.db, 
       keys = "b2378", 
       keytype = "ALIAS", 
       columns = c("GENENAME", "ENTREZID", "SYMBOL", "PATH"))
#####################################################################
entrez_ids <- bitr(annotated_variants$gene_id, fromType = "GENENAME",
                   toType = "ENTREZID", OrgDb = "org.EcK12.eg.db")
# Perform KEGG enrichment analysis
kegg_enrich <- enrichKEGG(gene = annotated_variants$gene_id, 
                          organism = "eco")
head(kegg_enrich)
dotplot(kegg_enrich, showCategory = 20)
emapplot(kegg_enrich)
cnetplot(kegg_enrich, showCategory = 10, circular = TRUE, colorEdge = TRUE)



