#Bulk RNA seq data analysis 
## Exercise
##Th0 vs Th2 separately for naive cells / memory cell
##compare DEGs for both cell types to find overlap and difference
## The biology behind it is to find difference in expression of multiple 
##genes between Naïve CD4 t-cells and memory t-cells after stimulating
##the cells to different types of cytokines 


##we started from a mapped data 
##loading libraries
library(DESeq2)
library(vsn)
## Loading the data
ge_matrix <- read.table('./bulk_RNAseq_raw_counts.txt', 
                        header = TRUE, sep = '\t')
dim(ge_matrix)
ge_matrix[1:4, 1:4]

## Loading the meta-data
pheno_matrix <- read.table("./bulk_RNAseq_metadata.txt", 
                           header = TRUE, sep = '\t', stringsAsFactors = TRUE)
pheno_matrix[1:6, 1:5]

## Organize the data
rownames(pheno_matrix) <- pheno_matrix$sample_id
#making the first column the raw names which is sample_id  
dim(pheno_matrix)
head(pheno_matrix)

#making sure the sample id in raw names of the pheno_matix is same
##and same order as the column names of ge_matix
all(rownames(pheno_matrix) == colnames(ge_matrix))

#selecting the sample we want to analyze
#select samples corresponding to the cell type and treatment that
#we want to focus on, in this case the CD4+ Memory cells
#after 5 days of treatment vs. control.
stimTime    <- '5d'
conditions  <- c('Th2', 'Th0')
celltype    <- 'CD4_Memory'

toSelect <- pheno_matrix$stimulation_time == stimTime & 
  pheno_matrix$cytokine_condition %in% conditions &
  pheno_matrix$cell_type == celltype
pheno_matrix.subset <- pheno_matrix[toSelect, ]
ge_matrix.subset <- ge_matrix[ , toSelect]
all(rownames(pheno_matrix.subset)==colnames(ge_matrix.subset))
# Create a DESeq2 Object
#DESeqDataSetFromMatrix is used in creating a DESeqDataSet object,
#which is the main data structure used in DESeq2 for storing count data 
#and sample metadata.
DESeq2_ds <- DESeqDataSetFromMatrix(countData = ge_matrix.subset, 
                              colData = pheno_matrix.subset,
                              design = ~ cytokine_condition)
# filtering of genes, which have too few
# counts, e.g. less than 10 reads over all samples.
keep <- rowSums(counts(DESeq2_ds)) >= 10
DESeq2_ds <- DESeq2_ds[keep,]
DESeq2_ds

# Investigate the data after removing low count

normtfd <- normTransform(DESeq2_ds)
# Compare mean to sd
meanSdPlot(assay(normtfd))
# calculate rlog values which stands for regulized log of transformed DESeq2_ds
##This transformation stabilizes the variance across the dataset and makes the 
##data more amenable to downstream analyses, including PCA.
rltfd <- rlog(DESeq2_ds, blind=FALSE)
meanSdPlot(assay(rltfd))

# Normalization
DESeq2_ds <- estimateSizeFactors(DESeq2_ds)
sizeFactors(DESeq2_ds)
plot(sizeFactors(DESeq2_ds), 
     colSums(counts(DESeq2_ds, normalized=F)), 
     xlab = 'Size factor',
     ylab = 'Total number of reads', 
     pch = 19)

#computing PCA
##PCA requires normal-distributed data. We can use the rlog function to
##transform our data by library size and apply log2 transformation. One
##way to perform the PCA is then using the function prcomp.
##PCA 
library(factoextra)
rltfd.pca <- prcomp(t(assay(rltfd)), scale = TRUE)

## PCA analysis
##evaluate the complexity of the data, how
##much of the variance is explained by the first component, the second,
##We can use a scree plot to visualize that
fviz_eig(rltfd.pca) #scree plot

fviz_pca_ind(rltfd.pca)
##plotting the PCA dimension

plotPCA(rltfd, intgroup = 'sequencing_batch',ntop=26656) 
##ploting the PCa 1 and 2 designed for RNA seq which can
##interpret the grouping variable witch is sequencing patch to assess and control 
##the batch effect
plotPCA(rltfd, intgroup = 'cytokine_condition')

#Differential Expression analysis
library(EnhancedVolcano)
library(pheatmap)
DESeq2_ds <- DESeq(DESeq2_ds)
res <- results(DESeq2_ds)#print result of Differential expression analysis
dim(res)
res
##How many significantly differentially expressed genes do we find for the
##current contrast Th2 vs Th0 in CD4+ memory cells?
sum(res$padj <= 0.01 & 
      abs(res$log2FoldChange) > 1, na.rm = TRUE)# i used pvalue of 0.01 
                                                #and log2FoldChange of 1


# Visualization 1: Volcano Plot
EV<- EnhancedVolcano(res, lab = rownames(res), 
                     x = 'log2FoldChange', y = 'padj', 
                     subtitle = 'Th2 vs Th0', labSize = 3, 
                     pCutoff = 0.01,
                     FCcutoff = 1,
                     drawConnectors = TRUE)
ggsave("enhanced_volcano_plot2.png", EV, width = 15, height = 10, units = "in", dpi = 400)
#saving png of the graph by using high dimension to display all labels

# Visualization 2: Heatmap
##first we select the high differentially expressed 
DEG.idx <- which(res$padj <= 0.01 & 
                   abs(res$log2FoldChange) > 1)
df <- as.data.frame(colData(DESeq2_ds)[,c("cytokine_condition","donor_id", "sequencing_batch")])

pheatmap(assay(rltfd)[DEG.idx,], annotation_col=df,
         treeheight_row = 0, treeheight_col = 0, scale = "row")
