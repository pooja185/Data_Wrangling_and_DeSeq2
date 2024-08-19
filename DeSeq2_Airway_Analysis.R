#Script to perform DESeq2 analysis

#setwd("C:/Users/pooja/Documents/BIoinformatics_learning/RNA-seq/raw_data/DeSeq2_project/Scripts")

#load libraries
library(DESeq2)
library(tidyverse)
library(airway)

#-------------------------------------------------------------------------------
# Step 1: Load data
#-------------------------------------------------------------------------------
# Loading counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)


# loading sample_info as column information
col_data <- read.csv('sample_info.csv')

# To make sure DESeq2 package does not throw any error, the row names in col_data should match with column names in counts_data
# Check if its correct or not

all(colnames(counts_data) %in% rownames(col_data))

#if true proceed to check if they are in same order
all(colnames(counts_data) == rownames(col_data))


#-------------------------------------------------------------------------------
# Step 2 : Construct DESeqDataSet object
#-------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = col_data,
                       design = ~ dexamethasone)
dds

#-------------------------------------------------------------------------------
# prefiltering (recommended step, not necessary but preferred) to remove rows with low gene counts
#-------------------------------------------------------------------------------

#Keeping rows that have atleast more than 10 reads
filtered <- rowSums(counts(dds)) >= 10
dds <- dds[filtered,]


#set factor level (treated vs untreated) by keeping one as baseline, usually control samples
dds$dexamethasone <- relevel(dds$dexamethasone, ref= "untreated")
dds$dexamethasone


#-------------------------------------------------------------------------------
# Step 3: Run DESeq
#-------------------------------------------------------------------------------
dds <- DESeq(dds)

#save results
res <- results(dds)

res

# Explore results

summary(res)

#results at alpha 0.01
res0.01 <- results(dds, alpha = 0.01) 
summary(res0.01)

resultsNames(dds)
#-----------------------------------------------------------------------------------
#for samples having multiple levels e.g. untreated, treated_4hrs, treated 20hrs etc use contrast

#contrasts

#results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

#-------------------------------------------------------------------------------
# Visualize results
#-------------------------------------------------------------------------------
#Creating (MA plot)

plotMA(res)

#Volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#PCA plot
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="dexamethasone") #using the DESEQ2 plotPCA fxn we can
