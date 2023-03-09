library("xgboost")

expression <- read.table("GSE125583_cpm_fix_new.txt", header = T, sep = "\t")
expression <- expression[c("FAM153B", "CYP2C8", "CKMT1B"),]
control <- as.data.frame(t(expression[,215:283]))
case <- as.data.frame(t(expression[,1:214]))
control$Group <- 0
case$Group <- 1
expression <- rbind(control, case)

expression = expression[sample(1:nrow(expression)), ]
expression = expression[sample(1:nrow(expression)), ]
expression = expression[sample(1:nrow(expression)), ]

lab <- expression$Group
expression$Group <- NULL

cv.dcg <- xgb.cv(data = as.matrix(expression), label = lab, nrounds = 10, nthread = 2, nfold = 5, metrics = list("rmse","auc", "error"),
             max_depth = 2, eta = 1, objective = "binary:logistic")
print(cv.dcg)


########

expression <- read.table("GSE125583_cpm_fix_new.txt", header = T, sep = "\t")
expression <- expression[c("PKN2", "FNDC3A", "NRIP1", "TMTC2"), ]
control <- as.data.frame(t(expression[,215:283]))
case <- as.data.frame(t(expression[,1:214]))
control$Group <- 0
case$Group <- 1
expression <- rbind(control, case)

expression = expression[sample(1:nrow(expression)), ]
expression = expression[sample(1:nrow(expression)), ]
expression = expression[sample(1:nrow(expression)), ]

lab <- expression$Group
expression$Group <- NULL

cv.hub <- xgb.cv(data = as.matrix(expression), label = lab, nrounds = 10, nthread = 2, nfold = 5, metrics = list("rmse","auc", "error"),
             max_depth = 2, eta = 1, objective = "binary:logistic")

######

expression <- read.table("GSE125583_cpm_fix_new.txt", header = T, sep = "\t")
expression <- expression[c("FAM153B", "CYP2C8", "CKMT1B", "PKN2", "FNDC3A", "NRIP1", "TMTC2"), ]
control <- as.data.frame(t(expression[,215:283]))
case <- as.data.frame(t(expression[,1:214]))
control$Group <- 0
case$Group <- 1
expression <- rbind(control, case)

expression = expression[sample(1:nrow(expression)), ]
expression = expression[sample(1:nrow(expression)), ]
expression = expression[sample(1:nrow(expression)), ]

lab <- expression$Group
expression$Group <- NULL

cv.all <- xgb.cv(data = as.matrix(expression), label = lab, nrounds = 10, nthread = 2, nfold = 5, metrics = list("rmse","auc", "error"),
             max_depth = 2, eta = 1, objective = "binary:logistic")

####

network <- read.csv("Network.txt", header = F)

expression <- read.table("GSE125583_cpm_fix_new.txt", header = T, sep = "\t")
expression <- expression[network$V1, ]
control <- as.data.frame(t(expression[,215:283]))
case <- as.data.frame(t(expression[,1:214]))
control$Group <- 0
case$Group <- 1
expression <- rbind(control, case)

expression = expression[sample(1:nrow(expression)), ]
expression = expression[sample(1:nrow(expression)), ]
expression = expression[sample(1:nrow(expression)), ]

lab <- expression$Group
expression$Group <- NULL

cv.net <- xgb.cv(data = as.matrix(expression), label = lab, nrounds = 10, nthread = 2, nfold = 5, metrics = list("rmse","auc", "error"),
                 max_depth = 2, eta = 1, objective = "binary:logistic")

cv.dcg.df <- as.data.frame(cv.dcg$evaluation_log)
cv.hub.df <- as.data.frame(cv.hub$evaluation_log)
cv.all.df <- as.data.frame(cv.all$evaluation_log)
cv.net.df <- as.data.frame(cv.net$evaluation_log)

sd(cv.net.df$train_auc_mean)

cv.dcg.df$Genes <- "Diff. Co-expressed genes"
cv.hub.df$Genes <- "Co-expressed hubs"
cv.all.df$Genes <- "DCGs + co-expressed hubs"
cv.net.df$Genes <- "Diff. co-expressed network"

cv.all <- Reduce(rbind, list(cv.dcg.df, cv.hub.df, cv.all.df, cv.net.df))

cv.all.target <- cv.all[, c("test_auc_mean", "test_error_mean", "Genes", "iter")]

library("reshape2")
library("ggpubr")
cv.all.target <- melt(cv.all.target, id.vars = c("Genes", "iter"))

sca <- ggscatter(cv.all, x = "test_error_mean", 
          y = "test_auc_mean", color = "Genes", 
          ggtheme = theme_bw(), shape = "Genes", size = 3,
          ylab = "Cross-validation - Test AUC mean", title = "Evaluation of XGBoost on prediction of AD", 
          xlab = "Cross-validation - Test error mean")
ggpar(sca, xlim =c(0.2, 0.35), ylim=c(0, 1))

summary(cv.net.df)







