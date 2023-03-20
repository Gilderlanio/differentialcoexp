# load libraries
library(stringr)
library(hgnc)
library(psych)
library(corto)
#source("my.diffcoexp.R")

diff.coexp <- function(project, samples, transcriptome, outdir, hgnc, target.genes) {
  
  allowWGCNAThreads()
  
  diagnosis <- unique(samples$diagnosis)
  tissues <- unique(samples$tissue)
  case <- NULL
  control <- NULL
  for (tissue in tissues) {
    for (diag in diagnosis)  {
      if (diag %in% c("AD", "CN")) {
        target.samples <- samples[samples$tissue == tissue & samples$diagnosis == diag, ]$specimenID
        target.samples <- intersect(colnames(transcriptome), target.samples)
        exp.data <- transcriptome[, target.samples]
        exp.data <- map.to.hgnc(exp.data, hgnc)
        exp.data <- exp.data[target.genes,]
        if (diag == "AD") {
          print(paste("Getting AD samples... ", dim(exp.data)[1], dim(exp.data)[2]))
          case <- exp.data
        } else {
          print(paste("Getting control samples... ", dim(exp.data)[1], dim(exp.data)[2]))
          control <- exp.data
        }
      }
    } # diag
    
    print(paste("Performing differential co-expression analysis. Tissue: ", tissue))
    target.genes <- intersect(rownames(control), rownames(case))
    control <- control[rownames(control) %in% target.genes, ]
    case <- case[rownames(case) %in% target.genes, ]
    
    control <- as.matrix(control[order(row.names(control)), ])
    case <- as.matrix(case[order(row.names(case)), ])

    if (dim(control)[2] >= 3 & dim(case)[2] >= 3) {
      res = diffcoexp(exprs.1 = control, exprs.2 = case, r.method = "spearman", q.method = "fdr", rth = 0.8, qth = 0.05)
      dcgs <- res$DCGs
      dcls <- res$DCLs
      
      write.table(dcgs, paste(paste(outdir, project, sep="/"), tissue, "DCGs.txt", sep = "_"),
                  row.names = T, sep = "\t", quote = F, append = F)
      write.table(dcls, paste(paste(outdir, project, sep="/"), tissue, "DCLs.txt", sep = "_"),
                  row.names = T, sep = "\t", quote = F, append = F)
      print("Done!")  
    }
    
  } # tissues
}

split.data <- function(project, samples, transcriptome, outdir, hgnc) {
  diagnosis <- unique(samples$diagnosis)
  tissues <- unique(samples$tissue)
  for (diag in diagnosis) {
    if (diag %in% c("AD", "CN")) {
      for (tissue in tissues) {
        target.samples <- samples[samples$tissue == tissue & samples$diagnosis == diag, ]$specimenID
        target.samples <- intersect(colnames(transcriptome), target.samples)
        exp.data <- transcriptome[,target.samples]
        exp.data <- map.to.hgnc(exp.data, hgnc)
        write.table(exp.data, paste(paste(outdir, project, sep="/"), diag, tissue, ".count", sep = "_"),
                    row.names = T, sep = "\t", quote = F, append = F)
      }
    }
  }
}

map.to.hgnc <- function(gene_expression, hgnc) {
  
  ensemblID <- row.names(gene_expression)
  
  # check for matches in hgnc data
  k <- ensemblID %in% hgnc$ensembl_gene_id
  
  ## generate dataset filtered by ensemblID match
  gene_expression <- gene_expression[k,]
  
  ## renaming rows to HGNC symbol
  k <- na.omit(match(ensemblID, hgnc$ensembl_gene_id))
  row.names(gene_expression) <- hgnc$symbol[k]
  
  ## Best statistical practice. Removing samples with total expression lower than 10M
  k <- colSums(gene_expression) >= 10^7
  gene_expression <- as.data.frame(gene_expression)
  gene_expression <- gene_expression[, k]
  
  # Normalizing data to counts per million
  gene_expression <- as.data.frame(edgeR::cpm(gene_expression))
  global_mean <- mean(as.matrix(gene_expression))
  ## removing transcripts with total expression lower than global mean
  k <- rowSums(gene_expression) >= global_mean
  gene_expression <- gene_expression[k,]
}