# install.packages("devtools")
# source("http://bioconductor.org/biocLite.R")
# biocLite("pcaMethods")
# 
# library(devtools)
# install_github("afukushima/DiffCorr")

comparecor <- function (exprs.1, exprs.2, r.method = c("pearson", "spearman")[1], 
          q.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", 
                       "BY", "fdr", "none")[1]) 
{
  require("psych")
  require("corto")
  require("DiffCorr")
  if (is(exprs.1, "SummarizedExperiment")) {
    exprs.1 <- assays(exprs.1)[[1]]
  }
  if (is(exprs.2, "SummarizedExperiment")) {
    exprs.2 <- assays(exprs.2)[[1]]
  }
  exprs.1 <- exprs.1[!is.na(rownames(exprs.1)), ]
  exprs.1 <- exprs.1[rownames(exprs.1) != "", ]
  exprs.2 <- exprs.2[!is.na(rownames(exprs.2)), ]
  exprs.2 <- exprs.2[rownames(exprs.2) != "", ]
  if (!all(rownames(exprs.1) == rownames(exprs.2))) {
    stop("Rownames of two expression matrices must be the same!")
  }
  genes <- rownames(exprs.1)
  exprs.1 <- as.matrix(exprs.1)
  exprs.2 <- as.matrix(exprs.2)
  if (sum(is.na(exprs.1)) == 0) {
    cor.1 <- cor(t(exprs.1), method = r.method, use = "all.obs")
    n.1 <- ncol(exprs.1)
  }
  else {
    cor.1 <- cor(t(exprs.1), method = r.method, use = "pairwise.complete.obs")
    n.1 <- pairwiseCount(t(exprs.1))
    n.1 <- n.1[lower.tri(n.1, diag = FALSE)]
  }
  if (sum(is.na(exprs.2)) == 0) {
    cor.2 <- cor(t(exprs.2), method = r.method, use = "all.obs")
    n.2 <- ncol(exprs.2)
  }
  else {
    cor.2 <- cor(t(exprs.2), method = r.method, use = "pairwise.complete.obs")
    n.2 <- pairwiseCount(t(exprs.2))
    n.2 <- n.2[lower.tri(n.2, diag = FALSE)]
  }
  cor.1 <- cor.1[lower.tri(cor.1, diag = FALSE)]
  cor.2 <- cor.2[lower.tri(cor.2, diag = FALSE)]
  rm(exprs.1)
  rm(exprs.2)
  name.row <- matrix(rep(genes, length(genes)), length(genes), 
                     length(genes))
  name.col <- matrix(rep(genes, length(genes)), length(genes), 
                     length(genes), byrow = TRUE)
  name.pairs <- matrix(paste(name.row, name.col, sep = ","), 
                       length(genes), length(genes))
  name.pairs <- name.pairs[lower.tri(name.pairs, diag = FALSE)]
  Gene.1 <- name.row[lower.tri(name.row, diag = FALSE)]
  Gene.2 <- name.col[lower.tri(name.col, diag = FALSE)]
  names(Gene.1) <- names(Gene.2) <- name.pairs
  rm(list = c("name.row", "name.col"))
  p.1 <- r2p(cor.1, n.1)
  p.2 <- r2p(cor.2, n.2)
  dc <- compcorr(n.1, cor.1, n.2, cor.2)
  res <- data.frame(Gene.1 = Gene.1, Gene.2 = Gene.2, cor.1 = cor.1, 
                    cor.2 = cor.2, cor.diff = cor.2 - cor.1, p.1 = p.1, 
                    p.2 = p.2, p.diffcor = dc$pval, stringsAsFactors = FALSE)
  res$q.1 <- p.adjust(res$p.1, method = q.method)
  res$q.2 <- p.adjust(res$p.2, method = q.method)
  res$q.diffcor <- p.adjust(res$p.diffcor, method = q.method)
  return(res)
}

coexpr <- function (exprs.1, exprs.2, r.method = c("pearson", "spearman")[1], 
          q.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", 
                       "BY", "fdr", "none")[1], rth = 0.5, qth = 0.1) 
{
  if (is(exprs.1, "SummarizedExperiment")) {
    exprs.1 <- assays(exprs.1)[[1]]
  }
  if (is(exprs.2, "SummarizedExperiment")) {
    exprs.2 <- assays(exprs.2)[[1]]
  }
  exprs.1 <- exprs.1[!is.na(rownames(exprs.1)), ]
  exprs.1 <- exprs.1[rownames(exprs.1) != "", ]
  exprs.2 <- exprs.2[!is.na(rownames(exprs.2)), ]
  exprs.2 <- exprs.2[rownames(exprs.2) != "", ]
  if (!all(rownames(exprs.1) == rownames(exprs.2))) {
    stop("Rownames of two expression matrices must be the same!")
  }
  x <- comparecor(exprs.1, exprs.2, r.method = r.method)
  if (!is.null(x)) {
    message("Finished running comparecor.")
  }
  x <- subset(x, subset = ((abs(x$cor.1) > rth & x$q.1 < qth) | 
                             (abs(x$cor.2) > rth & x$q.2 < qth)))
  return(x)
}

diffcoexp <- function(exprs.1, exprs.2, r.method = c("pearson", "kendall", 
                                         "spearman")[1], q.method = c("BH", "holm", "hochberg", "hommel", 
                                                                      "bonferroni", "BY", "fdr", "none")[1], rth = 0.5, qth = 0.1, 
          r.diffth = 0.5, q.diffth = 0.1, q.dcgth = 0.1) 
{
  if (is(exprs.1, "SummarizedExperiment")) {
    exprs.1 <- assays(exprs.1)[[1]]
  }
  if (is(exprs.2, "SummarizedExperiment")) {
    exprs.2 <- assays(exprs.2)[[1]]
  }
  exprs.1 <- exprs.1[!is.na(rownames(exprs.1)), ]
  exprs.1 <- exprs.1[rownames(exprs.1) != "", ]
  exprs.2 <- exprs.2[!is.na(rownames(exprs.2)), ]
  exprs.2 <- exprs.2[rownames(exprs.2) != "", ]
  if (!all(rownames(exprs.1) == rownames(exprs.2))) {
    stop("Rownames of two expression matrices must be the same!")
  }
  if (length(rownames(exprs.1)) == 0 | length(rownames(exprs.2)) == 
      0) {
    stop("The expression matrices must have row names specifying the gene names.")
  }
  if (min(ncol(exprs.1), ncol(exprs.2)) < 3) {
    stop("Each expression matrix must have at least three or more columns.")
  }
  else if (min(ncol(exprs.1), ncol(exprs.2)) < 5) {
    warning("The minimum number of columns is less than five and the result\n              may not be reliable.")
  }
  m <- nrow(exprs.1)
  genes <- rownames(exprs.1)
  colinks <- coexpr(exprs.1, exprs.2, r.method = r.method, 
                    rth = rth, qth = qth)
  if (!is.null(colinks)) {
    message("Finished running coexpr.")
  }
  if (nrow(colinks) == 0) {
    Result <- emptyresult()
    return(Result)
  }
  idx.same <- (colinks$cor.1 * colinks$cor.2) > 0
  idx.same[is.na(idx.same)] <- TRUE
  idx.diff <- (colinks$cor.1 * colinks$cor.2) < 0
  idx.diff[is.na(idx.diff)] <- FALSE
  idx.switched <- (colinks$cor.1 * colinks$cor.2 < 0) & (abs(colinks$cor.1) >= 
                                                           rth & abs(colinks$cor.2) >= rth & colinks$q.1 < qth & 
                                                           colinks$q.2 < qth)
  idx.switched[is.na(idx.switched)] <- FALSE
  cor.same <- colinks[idx.same, ]
  cor.switched <- colinks[idx.switched, ]
  cor.diff <- colinks[idx.diff & (!idx.switched), ]
  name.same <- NULL
  name.switched <- NULL
  name.diff <- NULL
  n.sameDCL <- 0
  if (nrow(cor.same) > 1) {
    idx.DCL.same <- cor.same$q.diffcor < q.diffth & abs(cor.same$cor.diff) > 
      r.diffth
    DCL.same <- cor.same[idx.DCL.same, ]
    name.same <- DCL.same[, c("Gene.1", "Gene.2")]
    n.sameDCL <- nrow(DCL.same)
  }
  else {
    DCL.same <- NULL
  }
  n.diffDCL <- 0
  if (nrow(cor.diff) > 1) {
    idx.DCL.diff <- cor.diff$q.diffcor < q.diffth & abs(cor.diff$cor.diff) > 
      r.diffth
    DCL.diff <- cor.diff[idx.DCL.diff, ]
    name.diff <- DCL.diff[, c("Gene.1", "Gene.2")]
    n.diffDCL <- nrow(DCL.diff)
  }
  else {
    DCL.diff <- NULL
  }
  n.switchedDCL <- 0
  if (nrow(cor.switched) > 1) {
    idx.DCL.switched <- cor.switched$q.diffcor < q.diffth & 
      abs(cor.switched$cor.diff) > r.diffth
    DCL.switched <- cor.switched[idx.DCL.switched, ]
    name.switched <- DCL.switched[, c("Gene.1", "Gene.2")]
    n.switchedDCL <- nrow(DCL.switched)
  }
  else {
    DCL.switched <- NULL
  }
  n.DCL <- n.sameDCL + n.diffDCL + n.switchedDCL
  message(nrow(colinks), " gene pairs remain after half thresholding.")
  if (n.DCL == 0) {
    message("No DCL meets the thresholds!")
    Result <- emptyresult()
    return(Result)
  }
  else {
    message(n.DCL, " DCLs identified.")
  }
  name.DCL <- rbind(name.same, name.diff, name.switched)
  name.colinks <- colinks[, c("Gene.1", "Gene.2")]
  g.colinks <- igraph::graph.data.frame(name.colinks)
  g.colinks.name <- as.matrix(igraph::V(g.colinks)$name)
  degree.colinks <- igraph::degree(g.colinks)
  g.DCL <- igraph::graph.data.frame(name.DCL)
  g.DCL.name <- as.matrix(igraph::V(g.DCL)$name)
  degree.DCL <- igraph::degree(g.DCL)
  if (n.sameDCL > 0) {
    g.same <- igraph::graph.data.frame(name.same)
    g.same.name <- as.matrix(igraph::V(g.same)$name)
    degree.same <- as.matrix(igraph::degree(g.same))
  }
  else {
    degree.same <- matrix(0, 1, 1)
  }
  if (n.diffDCL > 0) {
    g.diff <- igraph::graph.data.frame(name.diff)
    g.diff.name <- as.matrix(igraph::V(g.diff)$name)
    degree.diff <- as.matrix(igraph::degree(g.diff))
  }
  else {
    degree.diff <- matrix(0, 1, 1)
  }
  if (n.switchedDCL > 0) {
    g.switch <- igraph::graph.data.frame(name.switched)
    g.switch.name <- as.matrix(igraph::V(g.switch)$name)
    degree.switch <- as.matrix(igraph::degree(g.switch))
  }
  else {
    degree.switch <- matrix(0, 1, 1)
  }
  degree.bind <- data.frame(matrix(0, m, 5), stringsAsFactors = FALSE)
  row.names(degree.bind) <- genes
  colnames(degree.bind) <- c("CLs", "DCLs", "DCL.same", "DCL.diff", 
                             "DCL.switched")
  degree.bind[g.colinks.name, 1] <- degree.colinks
  degree.bind[g.DCL.name, 2] <- degree.DCL
  if (n.sameDCL > 0) {
    degree.bind[g.same.name, 3] <- degree.same
  }
  if (n.diffDCL > 0) {
    degree.bind[g.diff.name, 4] <- degree.diff
  }
  if (n.switchedDCL > 0) {
    degree.bind[g.switch.name, 5] <- degree.switch
  }
  prob <- nrow(name.DCL)/nrow(name.colinks)
  p.value <- pbinom(degree.bind[, "DCLs"] - 1, degree.bind[, 
                                                           "CLs"], prob, lower.tail = FALSE, log.p = FALSE)
  q.value <- p.adjust(p.value, method = q.method)
  degree.bind <- cbind(degree.bind, p.value, q.value)
  colnames(degree.bind) <- c("CLs", "DCLs", "DCL.same", "DCL.diff", 
                             "DCL.switch", "p", "q")
  DCGs <- degree.bind
  DCGs <- as.data.frame(DCGs)
  DCGs <- subset(DCGs, subset = q < q.dcgth)
  DCGs <- cbind(Gene = as.character(rownames(DCGs)), DCGs)
  DCGs$Gene <- as.character(DCGs$Gene)
  o <- order(DCGs$p)
  DCGs <- DCGs[o, ]
  message(length(DCGs$Gene), " DCGs identified.")
  DCLs <- data.frame()
  if (n.sameDCL > 0) {
    DCLs <- rbind(DCLs, data.frame(DCL.same, type = "same signed"))
  }
  if (n.diffDCL > 0) {
    DCLs <- rbind(DCLs, data.frame(DCL.diff, type = "diff signed"))
  }
  if (n.switchedDCL > 0) {
    DCLs <- rbind(DCLs, data.frame(DCL.switched, type = "switched opposites"))
  }
  DCLs$Gene.1 <- as.character(DCLs$Gene.1)
  DCLs$Gene.2 <- as.character(DCLs$Gene.2)
  Result <- list(DCGs = DCGs, DCLs = DCLs)
  return(Result)
}