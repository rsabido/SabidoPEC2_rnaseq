workingDir <-getwd()

set.seed(721917)

library(dplyr)
library(edgeR)
library(GenomicAlignments)
library(BiocParallel)
library(DESeq2)
library(ggbeeswarm)
library(DOSE)
library(enrichplot)
library(pheatmap)
library(RColorBrewer)
library(ReportingTools)


#Introducimos los datos
targets <- read.table(file.path(workingDir,"targets.csv"), header = TRUE, sep = ",")
counts <- read.table(file.path(workingDir,"counts.csv"), header = TRUE, sep = ";")
counts$X <- gsub("\\..*", "", counts$X, fixed = FALSE)

#Escogemos 10 de cada categoria
NIT_table <- filter(targets, Group == "NIT")
ELI_table <- filter(targets, Group == "ELI")
SFI_table <- filter(targets, Group == "SFI")
NIT.random <- NIT_table[sample(1:nrow(NIT_table), 10),]
ELI.random <- ELI_table[sample(1:nrow(ELI_table), 10),]
SFI.random <- SFI_table[sample(1:nrow(SFI_table), 10),]
selected_table <- rbind(NIT.random,ELI.random,SFI.random)

selected_counts <- as.matrix(counts %>% select(selected_table$Sample_Name))
rownames(selected_counts) <- counts$X

dds <- DESeqDataSetFromMatrix(countData = selected_counts,
                              colData = selected_table,
                              design = ~ Group)

#Hacemos un prefiltraje para eliminar todos aquellos genes con menos de 1 read. Nos quedamos con 46012 genes
quantile(rowSums(counts(dds)), c(0.33))
keep <- rowSums(counts(dds)) >= 7
dds <- dds[keep,]

dds <- DESeq(dds)

#Analisis de expresion diferencial

res_SFIvsELI <- results(dds, contrast = c("Group", "SFI", "ELI"))
res_NITvsELI <- results(dds, contrast = c("Group", "NIT", "ELI"))
res_NITvsSFI <- results(dds, contrast = c("Group", "NIT", "SFI"))

res05_SFIvsELI <- results(dds, alpha = 0.05, contrast = c("Group", "SFI", "ELI"))
table(res05_SFIvsELI$padj < 0.05)
summary(res05_SFIvsELI)

res05_NITvsELI <- results(dds, alpha=0.05, contrast = c("Group", "NIT", "ELI"))
table(res05_NITvsELI$padj < 0.05)
summary(res05_NITvsELI)

res05_NITvsSFI <- results(dds, alpha=0.05, contrast = c("Group", "NIT", "SFI"))
table(res05_NITvsSFI$padj < 0.05)
summary(res05_NITvsSFI)

#Anotación y exportación de resultados
library(AnnotationDbi)
library("org.Hs.eg.db")

res_SFIvsELI$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_SFIvsELI),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NITvsELI$symbol <- mapIds(org.Hs.eg.db,
                                keys=row.names(res_NITvsELI),
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first")

res_NITvsSFI$symbol <- mapIds(org.Hs.eg.db,
                                keys=row.names(res_NITvsSFI),
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first")


res_SFIvsELI$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_SFIvsELI),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NITvsELI$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_NITvsELI),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NITvsSFI$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_NITvsSFI),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")



#Exportamos
res_SFIvsELIDF <- as.data.frame(res_SFIvsELI)
res_NITvsELIDF <- as.data.frame(res_NITvsELI)
res_NITvsSFIDF <- as.data.frame(res_NITvsSFI)

htmlRep_SFIvsELI <- HTMLReport(shortName="report_SFIvsELI", title="My report SFIvsELI",
                               reportDirectory=".")
publish(res_SFIvsELIDF, htmlRep_SFIvsELI)
url_SFIvsELI <- finish(htmlRep_SFIvsELI)
browseURL(url_SFIvsELI)

htmlRep_NITvsELI <- HTMLReport(shortName="report_NITvsELI", title="My report NITvsELI",
                              reportDirectory=".")
publish(res_NITvsELIOrderedDF, htmlRep_NITvsELI)
url_NITvsELI <- finish(htmlRep_NITvsELI)
browseURL(url_NITvsELI)

htmlRep_NITvsSFI <- HTMLReport(shortName="report_NITvsSFI", title="My report NITvsSFI",
                               reportDirectory=".")
publish(res_NITvsSFIOrderedDF, htmlRep_NITvsSFI)
url_NITvsSFI <- finish(htmlRep_NITvsSFI)
browseURL(url_NITvsSFI)

#Comparaciones múltiples

res05_NITvsSFIDF <- as.data.frame(res05_NITvsSFI)
res05_NITvsSFIDF<- res05_NITvsSFIDF %>%
  rownames_to_column('ENSG') %>%
  filter(padj < 0.05)
genes_NITvsSFI <- res05_NITvsSFIDF$ENSG

res05_SFIvsELIDF <- as.data.frame(res05_SFIvsELI)
res05_SFIvsELIDF<- res05_SFIvsELIDF %>%
  rownames_to_column('ENSG') %>%
  filter(padj < 0.05)
genes_SFIvsELI <- res05_SFIvsELIDF$ENSG

res05_NITvsELIDF <- as.data.frame(res05_NITvsELI)
res05_NITvsELIDF<- res05_NITvsELIDF %>%
  rownames_to_column('ENSG') %>%
  filter(padj < 0.05)
genes_NITvsELI <- res05_NITvsELIDF$ENSG


venn.plot <- draw.triple.venn(
  area1 = length(genes_NITvsSFI),
  area2 = length(genes_SFIvsELI),
  area3 = length(genes_NITvsELI),
  n12 = length(intersect(genes_NITvsSFI, genes_SFIvsELI)),
  n23 = length(intersect(genes_NITvsELI, genes_SFIvsELI)),
  n13 = length(intersect(genes_NITvsSFI, genes_NITvsELI)),
  n123 = length(intersect(intersect(genes_NITvsSFI,genes_SFIvsELI),genes_NITvsELI)),
  category = c("NIT vs SFI", "SFI vs ELI", "NIT vs ELI"),
  fill = c("blue", "red", "green"),
  scaled=TRUE)

#Análisis de significación biológica

datasetNITvsSFI <- res_NITvsSFIDF$lfcSE
names(datasetNITvsSFI) <- as.character(res_NITvsSFIDF$entrez)
datasetNITvsSFI <- sort(datasetNITvsSFI, decreasing = TRUE)

datasetSFIvsELI <- res_SFIvsELIDF$lfcSE
names(datasetSFIvsELI) <- as.character(res_SFIvsELIDF$entrez)
datasetSFIvsELI <- sort(datasetSFIvsELI, decreasing = TRUE)

datasetNITvsELI <- res_NITvsELIDF$lfcSE
names(datasetNITvsELI) <- as.character(res_NITvsELIDF$entrez)
datasetNITvsELI <- sort(datasetNITvsELI, decreasing = TRUE)


geneNITvsSFI <- names(datasetNITvsSFI)[abs(datasetNITvsSFI) > 1]
geneNITvsSFI<-geneNITvsSFI[!is.na(geneNITvsSFI)]

geneSFIvsELI <- names(datasetSFIvsELI)[abs(datasetSFIvsELI) > 1]
geneSFIvsELI<-geneSFIvsELI[!is.na(geneSFIvsELI)]

geneNITvsELI <- names(datasetNITvsELI)[abs(datasetNITvsELI) > 1]
geneNITvsELI<-geneNITvsELI[!is.na(geneNITvsELI)]


edoNITvsSFI <- enrichDGN(geneNITvsSFI)
barplot(edoNITvsSFI, showCategory=10)

edoSFIvsELI <- enrichDGN(geneSFIvsELI)
barplot(edoSFIvsELI, showCategory=10)

edoNITvsELI <- enrichDGN(geneNITvsELI)
barplot(edoNITvsELI, showCategory=10)

edoxNITvsSFI <- setReadable(edoNITvsSFI, 'org.Hs.eg.db', 'ENTREZID')
heatplot(edoxNITvsSFI, foldChange=datasetNITvsSFI, showCategory = 10)

edoxSFIvsELI <- setReadable(edoSFIvsELI, 'org.Hs.eg.db', 'ENTREZID')
heatplot(edoxSFIvsELI, foldChange=datasetSFIvsELI, showCategory = 10)

edoxNITvsELI <- setReadable(edoNITvsELI, 'org.Hs.eg.db', 'ENTREZID')
heatplot(edoxNITvsELI, foldChange=datasetNITvsELI, showCategory = 10)

wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewpNITvsSFI <- enricher(geneNITvsSFI, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewpNITvsSFI[1:5,2:5])

ewpSFIvsELI <- enricher(geneSFIvsELI, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewpSFIvsELI[1:5,2:5])

ewpNITvsELI <- enricher(geneNITvsELI, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewpNITvsELI[1:5,2:5])