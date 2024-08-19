#Project : To explore gene expression in Breast cancer
#
# Script to manipulate gene expression data
#load libraries

library(dplyr)
library(tidyverse)
library(GEOquery)


# reading the data

raw_data <- read.csv(file="../raw_data/GSE183947_fpkm.csv")

# get metadata

gse <- getGEO(GEO='GSE183947', GSEMatrix = TRUE)
gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

#to avoid generating multiple variable for performing manipulation on dataset use"pipe function (%>%) of tidyverse"
modified_metadata <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%  #change column name of characteristics_ch1 (old) to tissue 
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%  #can create new column orchange values present in one column
  mutate(metastasis = gsub("metastasis: ", "", metastasis)) 


#reshape data
data_long <- raw_data %>%
  rename(gene = X) %>%
  gather(key='samples', value = 'FPKM', -gene)   #coverts y format to long format

# Join dataframes data_long and metadta modified to have information about all the samples under one roof


data_long <- data_long %>%
  left_join(., modified_metadata, by = c("samples" = "description"))




#Exploring the dataframe
#Since we are interested to know the gene expression of BRCA1 and BRCA2 we can just filter the information

data_long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM), 
            median_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM)

##----------------------------------------------------------------------------------------------------------------  
# Data visualization
##----------------------------------------------------------------------------------------------------------------
library(ggplot2)

#Creating barplot for gene of interest BRCA1 across all the samples

data_long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue)) + 
  geom_col() #do not provide column names in quotes


# Density plot

data_long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue)) + 
  geom_density(alpha = 0.3)


# Boxplot --- to confirm between different metastasis
data_long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM)) + 
  geom_boxplot()

#Scatterplot --comapre expression between 2 genes
data_long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  spread(key = gene, value= FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissue)) +
  geom_point() + 
  geom_smooth(method = 'lm', se = FALSE)


#heatmap
# For all multiple genes of interest, create a character vector

genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

p <- data_long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) + 
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')

ggsave(p, filename = 'heatmap_save1.pdf', width = 10, height = 8) #saves the file

#saving format no.2


data_long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')
