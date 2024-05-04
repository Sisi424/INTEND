<<<<<<< HEAD
version https://git-lfs.github.com/spec/v1
oid sha256:5019ea12bcee232f91573be678c1846cd6c500665c62404d227a21f64e084254
size 9342
=======
#1] Run libraries
library(rtracklayer)
library(dplyr)
library(parallel)
library(pracma)
library(glmnet)
##################################
#2] Load data
HG19.REF.GENE.PATH = "hg19.refGene.gtf"
GENE.COORDINATES.PATH = "gene.coordinates.rda"
HUMAN.METHYLATION.450.MANIFEST.PATH = "HumanMethylation450_15017482_v1-2.csv"
GENE.TO.CG.SITE.MAPPING.PATH.FORMAT = "gene_to_cpg_site_mapping.rda"
upstream_margin = 10000
downstream_margin = 10000
expression=readRDS("expression.rda")
methylation=readRDS("methylation.rda")
expression_features = rownames(expression)
methylation_features = rownames(methylation)
BRCA_Methylation = readRDS("BRCA_Methylation.rda")
###################################
#3] Functions
get.full.ref.gene.df <- function() 
{
  hg19.ref.gene <- rtracklayer::import(HG19.REF.GENE.PATH)
  hg19.ref.gene.df <- as.data.frame(hg19.ref.gene)
  return(hg19.ref.gene.df)
}
get.gene.coordinates.df <- function(force.update = FALSE) 
{
  if (file.exists(GENE.COORDINATES.PATH) & !force.update) {
    return(readRDS(file = GENE.COORDINATES.PATH))
  }
  
  hg19.ref.gene.df <- get.full.ref.gene.df()
  
  first.chr.index <- which(levels(hg19.ref.gene.df$seqnames) == "chr1")
  last.chr.index <- which(levels(hg19.ref.gene.df$seqnames) == "chrY")
  
  gene.coordinates.df <- hg19.ref.gene.df %>% 
    filter(type == "transcript" &
             as.integer(seqnames) >= first.chr.index & as.integer(seqnames) <= last.chr.index) %>%
    select(gene_id, chr = seqnames, strand, start, end) %>% 
    distinct()
  
  saveRDS(gene.coordinates.df, file = GENE.COORDINATES.PATH)
  return(gene.coordinates.df)
}
is.valid.value <- function(value)
{
  return(!is.na(value) & value != "")
}
utils.log <- function(log.str) {
  
  print(paste(format(Sys.time(), "[%X]"), log.str))
}
create.gene.to.cg.site.mapping <- function(
    gene.coordinates.df, upstream.margin.in.bases, downstream.margin.in.bases, force.update = FALSE)
{
  mapping.file.path <- 
    sprintf(GENE.TO.CG.SITE.MAPPING.PATH.FORMAT, upstream.margin.in.bases, downstream.margin.in.bases)
  
  if (file.exists(mapping.file.path) & !force.update) {
    return(readRDS(file = mapping.file.path))
  }
  
  methylation.manifest <- read.csv(file = HUMAN.METHYLATION.450.MANIFEST.PATH)
  
  methylation.manifest.filtered <- methylation.manifest %>%
    filter(is.valid.value(CHR) & is.valid.value(MAPINFO) & Genome_Build == 37)
  
  gene.to.cg.site.mapping <- data.frame(gene = unique(gene.coordinates.df$gene_id))
  gene.to.cg.site.mapping$cg_sites <- lapply(1:nrow(gene.to.cg.site.mapping), function(x) c())
  
  cg.sites.per.gene <- mclapply(1:nrow(gene.to.cg.site.mapping), function(gene_idx) {
    gene <- gene.to.cg.site.mapping$gene[gene_idx]
    gene.coordinates.filtered <- gene.coordinates.df %>% filter(gene_id == gene)
    
    #convert chromosome to string comparable with the one in methylation.manifest.filtered$CHR
    gene.coordinates.filtered$chr <- sapply(gene.coordinates.filtered$chr, function(chr.str) 
      substr(as.character(chr.str), start = nchar("chr") + 1, stop = nchar(as.character(chr.str))))
    
    gene.cg.sites <- c()
    
    for (i in 1:nrow(gene.coordinates.filtered)) {
      chr <- gene.coordinates.filtered$chr[i]
      
      if (gene.coordinates.filtered$strand[i] == "+") {
        
        start <- gene.coordinates.filtered$start[i] - upstream.margin.in.bases
        end <- gene.coordinates.filtered$end[i] + downstream.margin.in.bases
        
      } else {
        
        start <- gene.coordinates.filtered$start[i] - downstream.margin.in.bases
        end <- gene.coordinates.filtered$end[i] + upstream.margin.in.bases
      }
      
      cg.sites.in.range <- (methylation.manifest.filtered %>% 
                              filter(tolower(CHR) == tolower(chr) & MAPINFO >= start & MAPINFO <= end) %>%
                              select(Name))$Name
      
      gene.cg.sites <- c(gene.cg.sites, cg.sites.in.range)
    }
    
    if (gene_idx %% 1000 == 1) {
      utils.log(sprintf("Done mapping gene number %d", gene_idx - 1))
    }
    
    return(gene.cg.sites)
  })
  
  gene.to.cg.site.mapping$cg_sites <- lapply(cg.sites.per.gene, unique)
  
  saveRDS(gene.to.cg.site.mapping, file = mapping.file.path)
  return(gene.to.cg.site.mapping)
}
train_predict_lasso <- function(gene, cg_sites, expression_data, methylation_data) {
  # Subset methylation data for relevant CpG sites and transpose for glmnet
  X_train <- t(methylation_data[cg_sites, , drop = FALSE])
  
  # Get corresponding expression data for the gene
  Y_train <- expression_data[gene, ]
  
  # Ensure Y_train is a numeric vector
  Y_train <- as.numeric(Y_train)
  
  # Train Lasso model with cross-validation
  model <- cv.glmnet(x = as.matrix(X_train), y = Y_train, alpha = 1, nfolds = 10)
  
  # Predict gene expression
  predictions <- predict(model, newx = as.matrix(X_train), s = "lambda.min")
  
  # Calculate R-squared
  R_squared <- cor(Y_train, predictions)^2
  
  # Return a named list including the gene, model, predictions, and R-squared
  return(list(Gene = gene,
              Fit_Parameters = model, 
              Predicted_Expression = predictions, 
              R_squared = R_squared))
}

# Wrapper function to apply the Lasso model to each gene and return a data frame
apply_lasso_to_genes <- function(expression_data, methylation_data, gene_to_cpg_site_mapping) {
  results_list <- lapply(1:nrow(gene_to_cpg_site_mapping), function(i) {
    gene <- gene_to_cpg_site_mapping$gene[i]
    cg_sites <- gene_to_cpg_site_mapping$cg_sites[[i]]
    result <- train_predict_lasso(gene, cg_sites, expression_data, methylation_data)
    result
  })
  
  # Convert the list of results to a data frame
  results_df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(Gene = x$Gene,
               Fit_Parameters = I(list(x$Fit_Parameters)), 
               Predicted_Expression = I(list(x$Predicted_Expression)), 
               R_squared = x$R_squared)
  }))
  colnames(results_df) <- c("Gene", "Fit_Parameters", "Predicted_Expression", "R_squared")
  
  colnames(results_df) <- c("Gene", "Fit_Parameters", "Predicted_Expression", "R_squared")
  return(results_df)
}
predict_expression <- function(methylation, gene_to_cpg_site_mapping, regression.output) {
  predicted_expression_list <- list()
  
  for (gene in rownames(regression.output)) {
    cg_sites <- gene_to_cpg_site_mapping$cg_sites[gene_to_cpg_site_mapping$gene == gene][[1]]
    
    if (length(cg_sites) == 0) {
      cat("No CpG sites found for gene:", gene, "\n")
      predicted_expression_list[[gene]] <- rep(NA, ncol(methylation))
      next
    }
    
    fit <- regression.output$Fit_Parameters[[which(rownames(regression.output) == gene)]]
    
    if (is.null(fit)) {
      cat("Invalid or NULL model for gene:", gene, "\n")
      predicted_expression_list[[gene]] <- rep(NA, ncol(methylation))
      next
    }
    
    if (!all(cg_sites %in% rownames(methylation))) {
      cat("Not all CpG sites found in methylation data for gene:", gene, "\n")
      predicted_expression_list[[gene]] <- rep(NA, ncol(methylation))
      next
    }
    
    X <- t(methylation[cg_sites, , drop = FALSE])
    gene_prediction <- predict(fit, newx = X, s = "lambda.min")
    predicted_expression_list[[gene]] <- as.vector(gene_prediction)
  }
  
  predicted_expression <- do.call(cbind, predicted_expression_list)
  rownames(predicted_expression) <- colnames(methylation)
  colnames(predicted_expression) <- names(predicted_expression_list)
  
  return(predicted_expression)
}
#################################################
#4] Call functions

gene_coords_df <- get.gene.coordinates.df()

gene_to_cpg_site_mapping <- create.gene.to.cg.site.mapping(gene.coordinates.df =  gene_coords_df,upstream.margin.in.bases =  upstream_margin,downstream.margin.in.bases =  downstream_margin)

gene_to_cpg_site_mapping <- gene_to_cpg_site_mapping %>% 
  filter(gene %in% expression_features)

gene_to_cpg_site_mapping$cg_sites <- lapply(gene_to_cpg_site_mapping$cg_sites, function(cpg_sites_for_gene) 
  intersect(cpg_sites_for_gene, methylation_features))

cpg_sites_lengths <- sapply(gene_to_cpg_site_mapping$cg_sites, length)

gene_to_cpg_site_mapping <- gene_to_cpg_site_mapping[which(cpg_sites_lengths > 1), ]

saveRDS(gene_to_cpg_site_mapping, "gene_to_cpg_site_mapping.rda")

regression.output <- apply_lasso_to_genes(expression_data = expression, methylation_data = methylation, gene_to_cpg_site_mapping = gene_to_cpg_site_mapping)

row.names(regression.output) <- regression.output[, 1]; regression.output <- regression.output[, -1]

saveRDS(regression.output, "regression.output.rda")

Genes_R_Squared <- regression.output[, 3, drop = FALSE]  

rownames(Genes_R_Squared) <- rownames(regression.output)  

saveRDS(Genes_R_Squared, "Genes_R_Squared.rda")

predicted_expression <- predict_expression(methylation=BRCA_Methylation, gene_to_cpg_site_mapping, regression.output=regression.output)

predicted_expression=t(predicted_expression)

predicted_expression=as.data.frame(predicted_expression)

saveRDS(predicted_expression,"Predicted.rda")


>>>>>>> 6439888a81728d3683f31c96e4077c49b53cecda
