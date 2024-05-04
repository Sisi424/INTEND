<<<<<<< HEAD
version https://git-lfs.github.com/spec/v1
oid sha256:9e88b63d33dad1d391e1d852282d385b5400097d489292667e2d80a48a8a6f97
size 2838
=======
#1] Import libraries
library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(SummarizedExperiment)
library(sesameData)
library(sesame)
############################################
#2] Set functions to download and label the TCGA data
clean_and_rename_columns <- function(df, common_samples) {
  # Extract the simplified column names (up to -01, inclusive)
  simplified_colnames <- sapply(colnames(df), function(name) {
    sub("^(TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-01).*", "\\1", name)
  })
  
  # Create a logical vector indicating whether each simplified name is in common_samples
  in_common_samples <- simplified_colnames %in% common_samples
  
  # Filter the dataframe columns based on the match to common_samples
  df <- df[, in_common_samples]
  
  # Rename columns to their simplified form, ensuring they exactly match common_samples
  colnames(df) <- simplified_colnames[in_common_samples]
  
  # Remove duplicate columns by keeping only the first occurrence of each unique identifier
  df <- df[, !duplicated(colnames(df))]
  
  return(df)
}
############################################
#3] Define a vector of cancer types to loop over
cancer_types <- c("BLCA","COAD","KIRP","LGG","LIHC","LUAD","MESO","PAAD","READ","SARC","STAD","UCEC","PRAD", "THCA", "SKCM")
############################################
#4] Loop through each cancer type and download its RNA and Methylation data common samples
for(cancer_type in cancer_types) {
  common_samples <- matchedMetExp(paste0("TCGA-", cancer_type))
  
  Methylation <- GDCquery(project = paste0("TCGA-", cancer_type),
                          data.category = "DNA Methylation",
                          platform = "Illumina Human Methylation 450",
                          data.type = "Methylation Beta Value",
                          barcode = common_samples)
  
  mRNA <- GDCquery(project = paste0("TCGA-", cancer_type),
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   barcode = common_samples)
  
  GDCdownload(Methylation)
  GDCdownload(mRNA)
  
  data_rnaseq <- GDCprepare(mRNA)
  data_methylation <- GDCprepare(Methylation)
  
  geneExp <- SummarizedExperiment::assay(data_rnaseq)
  geneExp = as.data.frame(geneExp)
  geneExp <- clean_and_rename_columns(geneExp, common_samples)
  
  MethExp <- SummarizedExperiment::assay(data_methylation)
  MethExp = as.data.frame(MethExp)
  MethExp <- clean_and_rename_columns(MethExp, common_samples)
  
  # Use getwd() to dynamically generate the path for saving files in the current working directory
  write.csv(geneExp, file = paste0(getwd(), "/", cancer_type, "_RNA.csv"))
  write.csv(MethExp, file = paste0(getwd(), "/", cancer_type, "_Methylation.csv"))
}
>>>>>>> 6439888a81728d3683f31c96e4077c49b53cecda
