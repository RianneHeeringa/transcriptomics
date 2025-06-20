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
group <- c("normaal", "normaal", "normaal", "normaal", "reuma", "reuma", "reuma", "reuma")
treatment_table <- data.frame(treatment, group)
rownames(treatment_table) <- c("SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831", "SRR4785979", "SRR4785980", "SRR4785986", "SRR4785988")
#packages laden
BiocManager::install("DESeq2") 
BiocManager::install ("KEGGREST")
library(DESeq2)
library(KEGGREST)
#statestiek
#grote matrix inladen
grote_matrix <- read.table("c:/Users/riann/OneDrive - NHL Stenden/BML/count_matrix.txt")
dds <- DESeqDataSetFromMatrix(countData = grote_matrix, colData = treatment_table, design = ~group)
grote_matrix <- round(grote_matrix)
dds <- DESeq(dds)
resultaten <- results(dds)
write.table(resultaten, file = 'res', row.names = TRUE, col.names = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
head(laagste_p_waarde)
#volcano plot
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)
dev.copy(png, 'VolcanoplotWC.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()
# pathview
#packages instaleren
BiocManager::install("goseq")
BiocManager::install("org.Hs.eg.db")  # ⬅ menselijke annotatie
BiocManager::install("geneLenDataBase")  # lengte-informatie
BiocManager::install("GO.db")  
library(goseq)
library(org.Hs.eg.db)
library(geneLenDataBase)
library(GO.db)
resu <- read.table(resultaten)


all <- rownames(resultaten)
deg <- resultaten %>%
  filter(padj < 0.05)
deg <- rownames(deg)

class(deg)
head(class(deg))
deg.vector <- c(t(deg))
all.vector <- c(t(all))
gene.vector=as.integer(all.vector%in%deg.vector)
names(gene.vector)=all.vector
head(gene.vector)
tail(gene.vector)

pwf=nullp(gene.vector, "hg19", "geneSymbol")
GO.wall=goseq(pwf, "hg19", "geneSymbol")
class(GO.wall)
head(GO.wall)
nrow(GO.wall)

enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
class(enriched.GO)
head(enriched.GO)
length(enriched.GO)

library(GO.db)
capture.output(for(go in enriched.GO[1:258]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}, file="SigGo.txt")

#visualiseren
install.packages("ggplot2")
library(ggplot2)

dim(GO.wall)
mode(GO.wall)
GO.wall_mat=as.matrix(GO.wall)
mode(GO.wall_mat)="matrix"
class(GO.wall_mat)

termbeschrijving <- sapply(enriched.GO, function(x) {
  if (!is.null(GOTERM[[x]])) Term(GOTERM[[x]]) else NA
})
enriched_df <- GO.wall[GO.wall$category %in% enriched.GO, ]
enriched_df$Term <- termbeschrijving[match(enriched_df$category, names(termbeschrijving))]
enriched_df$logp <- -log10(enriched_df$over_represented_pvalue)
top_terms <- enriched_df[order(enriched_df$over_represented_pvalue), ][1:10, ]
ggplot(top_terms, aes(x = reorder(Term, logp), y = logp)) +
  geom_bar(stat = "identity", fill = "pink") +
  coord_flip() +
  labs(title = "Top 10 verrijkte GO-termen",
       x = "GO-term beschrijving",
       y = "-log10(p-waarde)") +
  theme_minimal()
#KEGGgraph
BiocManager::install("clusterProfiler")
BiocManager::install("pathview") 
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
resultaten[1] <- NULL
resultaten[2:5] <- NULL
pathview (gene.data = gene.vector, pathway.id = "hsa05323", species = "hsa", gene.idtype = "SYMBOL", limit = (gene = 5))
