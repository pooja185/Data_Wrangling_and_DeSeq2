#Data to get from airway package R

#This package provides a RangedSummarizedExperiment object of read counts 
#in genes for an RNA-Seq experiment on four human airway smooth muscle cell lines treated with dexamethasone. 
#Details on the gene model and read counting procedure are provided in the package vignette. 
#The citation for the experiment is: Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q, 
#Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr, Tantisira KG, Weiss ST, Lu Q. 
#'RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates 
#'#Cytokine Function in Airway Smooth Muscle Cells.' PLoS One. 2014 Jun 13;9(6):e99625. PMID: 24926665. GEO: GSE52778.



library(airway)

data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)