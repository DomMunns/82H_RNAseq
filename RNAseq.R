


#download BiocManager and devtools
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version="3.12")
BiocManager::install("rhdf5")

if(!require(devtools)){
  install.packages("devtools")
}
library("devtools")

#install and load packages
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("Biobase")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("airway")

library(limma)
library(edgeR)
library(Biobase)
library(org.Hs.eg.db)
library(ggrepel)
library(tidyverse)
library(tidybulk)
library(ggplot2)
library(gplots)

#read in static vs oss data
df <- read.csv("data_raw/data_01rpm.csv")
colnames(df)
head(df)

# Normalise the dataset
norm_df <- normalizeBetweenArrays(df[,c(4:11)])
names_df <- df[-c(4:11)]
norm_df1 <- cbind(names_df, norm_df)

#add in row for group
df1 <- list(0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2)

df <- rbind(df, df1)

#specify condition group (static = 1, osc = 2)
dge <- DGEList(counts = df[, c(4:7, 8:11)], genes = df[, 2], group = df[29153,4:11 ])
dge$samples

#map gene symbols and chromosomes to DGElist using org.Hs.eg.db package
columns(org.Hs.eg.db)

key <- dge$genes$genes
head(dge$genes)
#gene symbol
dge$genes$SYMBOL <- mapIds(org.Hs.eg.db, row.names(dge), keys = key,
                           keytype = "SYMBOL", column = "SYMBOL")
#entrez ID
dge$genes$entrez <- mapIds(org.Hs.eg.db, row.names(dge), keys = key,
                           keytype = "SYMBOL", column = "ENTREZID")
#gene function
dge$genes$genefunction <- mapIds(org.Hs.eg.db, row.names(dge),keys = key,
                                 keytype = "SYMBOL", column = "GENENAME")

head(dge$genes)

####################
#set up the expression set
####################

####expression matrix (use original df)
ex <- subset(df, select=-c(�..Ensembl.ID, Gene.type))
head(ex)

#make genes the row names
rownames(ex) <- make.names(ex$Gene.symbol, unique = TRUE)
ex <- ex[, -1]
head(ex)

#name data type numeric and matrix
ex.mat = as.matrix(as.data.frame(lapply(ex, as.numeric)))
ex.mat <- ex.mat[-29153, ]

#check transform has worked
str(ex.mat)
colnames(ex.mat)
class(ex.mat)


####feature data (use DGEList)
feat <- AnnotatedDataFrame(dge$genes)

#make genes the row names
rownames(feat) <- make.names(feat$genes, unique = TRUE)
feat <- feat[, -1]
feat <- feat[-29153,]


#check data frame worked
colnames(feat)
head(feat)
class(feat)



####phenotype data (use DGEList)
#add conditions
dge$samples$condition <- c(rep("Oscillatory Shear Stres",4), rep("Laminar Shear Stress", 4))

pheno <- subset(dge$samples, select=-c(lib.size, norm.factors))

pheno <- AnnotatedDataFrame(pheno)
#check data frame worked
head(pheno)
colnames(pheno)
class(pheno)


########################
#make the expression set
########################

eset <- ExpressionSet(assayData=ex.mat,
              phenoData = pheno,
              featureData = feat)


#log transform
exprs(eset) <- log(exprs(eset))
plotDensities(eset, group = pData(eset)[, "condition"], legend = "topright")

#Quantile normalise
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset, group = pData(eset)[, "condition"], legend = "topright")

#determine the genes with mean expression level greater than zero
keep <- rowMeans(exprs(eset)) > 0
sum(keep)

#filter the genes
eset <- eset[keep, ]
plotDensities(eset, group = pData(eset)[, "condition"], legend = "topright")

###now to do the analysis for my genes (VCAM1 and CCL2)
###VCAM1
#find row containing VCAM1
VCAM1 <- which(fData(eset)[, "SYMBOL"] == "VCAM1")
VCAM1

#plot VCAM1 expression vs condition
boxplot(exprs(eset)[VCAM1, ] ~ pData(eset)[, "condition"],
        main = fData(eset)[VCAM1, "SYMBOL"],
        xlab="Stress Condition",
        ylab="Gene Expression")

###CCL2
#find row containing CCL2
CCL2 <- which(fData(eset)[, "SYMBOL"] == "CCL2")
CCL2

#plot CCL2 expression vs condition
boxplot(exprs(eset)[CCL2, ] ~ pData(eset)[, "condition"],
        main = fData(eset)[CCL2, "SYMBOL"],
        xlab="Stress Condition",
        ylab="Gene Expression")



#plot priniciple components labelled by condition
plotMDS(eset, labels = pData(eset)[, "condition"], gene.selection = "common")

colnames(pData(eset))

#make the design matrix
design <- model.matrix(~condition, data = pData(eset))
colSums(design)
table(pData(eset)[, "condition"])


################
###limma time
################
#fit the model
fit <- lmFit(eset, design)

#calculate the t-statistics
fit <- eBayes(fit)

colnames(fit)

#summarize results
results <- decideTests(fit[, "conditionStatic Culture"])

summary(results)
#286 downregulated genes in static
#22959 genes with no significant change between conditions
#179 upregulated genes in static condition



### Genes

###ADAM17
#find row containing VCAM1
ADAM17 <- which(fData(eset)[, "SYMBOL"] == "ADAM17")
ADAM17

boxplot(exprs(eset)[ADAM17, ] ~ pData(eset)[, "condition"],
        main = fData(eset)[ADAM17, "SYMBOL"],
        xlab="Stress Condition",
        ylab="Gene Expression")


###CXCL3
CXCL3 <- which(fData(eset)[, "SYMBOL"] == "CXCL3")
CXCL3

boxplot(exprs(eset)[CXCL3, ] ~ pData(eset)[, "condition"],
        main = fData(eset)[CXCL3, "SYMBOL"],
        xlab="Stress Condition",
        ylab="Gene Expression")


###ATF4
ATF4 <- which(fData(eset)[, "SYMBOL"] == "ATF4")
ATF4

boxplot(exprs(eset)[ATF4, ] ~ pData(eset)[, "condition"],
        main = fData(eset)[ATF4, "SYMBOL"],
        xlab="Stress Condition",
        ylab="Gene Expression")

###eIF2a
eIF2a <- which(fData(eset)[, "SYMBOL"] == "EIF2A")
eIF2a

boxplot(exprs(eset)[eIF2a, ] ~ pData(eset)[, "condition"],
        main = fData(eset)[eIF2a, "SYMBOL"],
        xlab="Stress Condition",
        ylab="Gene Expression")


###DDIT3
DDIT3 <- which(fData(eset)[, "SYMBOL"] == "DDIT3")
DDIT3

boxplot(exprs(eset)[DDIT3, ] ~ pData(eset)[, "condition"],
        main = fData(eset)[DDIT3, "SYMBOL"],
        xlab="Stress Condition",
        ylab="Gene Expression")

###ERN1
ERN1 <- which(fData(eset)[, "SYMBOL"] == "ERN1")
ERN1

boxplot(exprs(eset)[ERN1, ] ~ pData(eset)[, "condition"],
        main = fData(eset)[ERN1, "SYMBOL"],
        xlab="Stress Condition",
        ylab="Gene Expression")

###XBP1
XBP1 <- which(fData(eset)[, "SYMBOL"] == "XBP1")
XBP1

boxplot(exprs(eset)[XBP1, ] ~ pData(eset)[, "condition"],
        main = fData(eset)[XBP1, "SYMBOL"],
        xlab="Stress Condition",
        ylab="Gene Expression")

### Heat map
# convert data frame into a matrix

gene <- norm_df1[c(5490, 12303, 9803, 3311, 7884, 2100, 10075, 12821, 8578), c(1:11)]

# name the rows after the gene name
row.names(gene) <- make.names(gene$Gene.symbol, unique = TRUE)

# remove data columns that unnecessary 
gene <- gene[, 4:11]

#log transform
log_gene <-  log2(gene)

# make the dataframe a matrix so a heatmap could be produced
gene_matrix <- as.matrix(log_gene)

# heatmap

gene_heatmap <- heatmap.2(gene_matrix, notecol = "black", density.info = "none", 
                           trace = "none", margins = c(12,9), dendrogram = "row", 
                           Colv = NA, ylab = "Gene", xlab = "Sample")









