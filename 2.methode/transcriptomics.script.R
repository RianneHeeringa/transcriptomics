#instellen working direction
setwd("C:/Users/riann/OneDrive - NHL Stenden/BML/transciptomics")
getwd()
# downloaden benodigde packages
install.packages('BiocManager')
BiocManager::install('Rsubread')
library('Rsubread')
browseVignettes('Rsubread')
#reads mappen
setwd("C:/Users/riann/OneDrive - NHL Stenden/BML/transciptomics casus/")
#indexeren
buildindex(
  basename = 'human',
  reference = 'ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)
# humane monsters mapping
align.human <- align(index = "human" , readfile1 = "Data_RA_raw/SRR4785819_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785819_2_subset40k.fastq" , output_file = "SRR4785819")
align.human <- align(index = "human" , readfile1 = "Data_RA_raw/SRR4785820_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785820_2_subset40k.fastq" , output_file = "SRR4785820")
align.human <- align(index = "human" , readfile1 = "Data_RA_raw/SRR4785828_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785828_2_subset40k.fastq" , output_file = "SRR4785828")
align.human <- align(index = "human" , readfile1 = "Data_RA_raw/SRR4785831_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785831_2_subset40k.fastq" , output_file = "SRR4785831")
align.human <- align(index = "human" , readfile1 = "Data_RA_raw/SRR4785979_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785979_2_subset40k.fastq" , output_file = "SRR4785979")
align.human <- align(index = "human" , readfile1 = "Data_RA_raw/SRR4785980_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785980_2_subset40k.fastq" , output_file = "SRR4785980")
align.human <- align(index = "human" , readfile1 = "Data_RA_raw/SRR4785986_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785986_2_subset40k.fastq" , output_file = "SRR4785986")
align.human <- align(index = "human" , readfile1 = "Data_RA_raw/SRR4785988_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785988_2_subset40k.fastq" , output_file = "SRR4785988")
# BAM bestanden generen
library(Rsamtools)
sample <- c('SRR4785819', 'SRR4785820', 'SRR4785828', 'SRR4785831', 'SRR4785979', 'SRR4785980', 'SRR4785986', 'SRR4785988')
# count matrix
library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)
# Kolomnamen toevoegen
colnames() <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# feuturecounts uivoeren
# Je definieert een vector met namen van BAM-bestanden. Elke BAM bevat reads van een RNA-seq-experiment (bijv. behandeld vs. controle).

allsamples <- c("SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831", "SRR4785979", "SRR4785980", "SRR4785980", "SRR4785986", "SRR4785988")

count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "GCF_000001405.40_GRCh38.p14_genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)
#resultaten bekijken
head(count_matrix$annotation)
head(count_matrix$counts)
str(count_matrix)
counts <- count_matrix$counts
colnames(counts) <- c("GeneID", "SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831", "SRR4785979", "SRR4785980", "SRR4785986", "SRR4785988")
rownames(counts) <- counts[,1]
counts <- counts[, -1]
write.csv(counts, "bewerkt_countmatrix.csv")
head(counts)
counts <- read.csv("bewerkt_countmatrix.csv", row.names = NULL)
#behandelingstabel maken
treatment <- c("SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831", "SRR4785979", "SRR4785980", "SRR4785986", "SRR4785988")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c("SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831", "SRR4785979", "SRR4785980", "SRR4785986", "SRR4785988"
BiocManager::install("DESeq2") 
BiocManager::install ("KEGGREST")
library(DESeq2)
library(KEGGREST)
