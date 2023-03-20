# load libraries
library(stringr)
library(hgnc)

map.to.hgnc <- function(gene_expression, hgnc) {
  
  # Mapping ensembl ids to HGNC symbol.
  ensemblID <- row.names(gene_expression)
  ensemblID <- t(as.data.frame(str_split(ensemblID, pattern = "[.]")))
  ensemblID <- ensemblID[,1]
  
  # check for matches in hgnc file
  k <- ensemblID %in% hgnc_full_list$ensembl_gene_id

  ## generate dataset filtered by ensemblID match
  gene_expression <- gene_expression[k,]
  
  ## renaming rows to HGNC symbol
  k <- na.omit(match(ensemblID, hgnc_full_list$ensembl_gene_id))
  row.names(gene_expression) <- hgnc_full_list$symbol[k]
  
  # Normalizing data to counts per million
  gene_expression <- edgeR::cpm(gene_expression)
  
  global_mean <- mean(gene_expression)
  
  ## removing transcripts with total expression lower than global mean
  k <- rowSums(input_norm) >= global_mean
  gene_expression <- gene_expression[k,]
  
  ## Best statistical practice. removing samples with total expression lower than 10M
  k <- colSums(gene_expression) >= 10^7
  gene_expression <- gene_expression[, k]

  return(gene_expression)
  
}