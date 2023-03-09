library(ggplot2)
library(cowplot)
library(stringr)

path_ad <- read.table("12_sig_pathways_ad.txt", header= T, sep = '\t')
path_cn <- read.table("12_sig_pathways_cn.txt", header= T, sep = '\t')
path_ad$Group <- "AD"
path_cn$Group <- "CN"
paths <- rbind(path_ad, path_cn)
exc <- setdiff(path_ad$ID, path_cn$ID)
exc12 <- exc
paths <- paths[paths$ID %in% exc,]
paths <- paths[order(paths$qvalue, decreasing = F),][1:10,]
ad_exc_12 <- ggplot(paths, aes(y = -log10(padjust), x = reorder(ID, -log10(padjust)), color = Module,)) +
  geom_point(stat = "identity", aes(size = Count)) + xlab("Pathway") + ylab("-log10(adj. p-value)") +
  geom_segment(aes(x=ID, xend=ID, y=0, yend=-log10(padjust))) + ylim(0, 5) +
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(text = element_text(size=10, color ="black", face = "bold")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + coord_flip() + facet_wrap(~Group)
ad_exc_12

path_ad <- read.table("12_sig_pathways_ad.txt", header= T, sep = '\t')
path_cn <- read.table("12_sig_pathways_cn.txt", header= T, sep = '\t')
path_ad$Group <- "AD"
path_cn$Group <- "CN"
paths <- rbind(path_ad, path_cn)
exc <- setdiff(path_cn$ID, path_ad$ID)
paths <- paths[paths$ID %in% exc,]
paths <- paths[order(paths$qvalue, decreasing = F),][1:10,]
cn_exc_12 <- ggplot(paths, aes(y = -log10(padjust), x = reorder(ID, -log10(padjust)), color = Module,)) +
  geom_point(stat = "identity", aes(size = Count)) + xlab("Pathway") + ylab("-log10( p-value adjust)") +
  geom_segment(aes(x=ID, xend=ID, y=0, yend=-log10(padjust))) + ylim(0, 5) +
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(text = element_text(size=10, color ="black", face = "bold")) + scale_color_manual(values = c('blue')) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + coord_flip() + facet_wrap(~Group)
cn_exc_12

plot_grid(cn_exc_12, ad_exc_12, nrow = 2)
ggsave("Final_12.pdf", scale = 1, width = 4, height = 8, units = "in", dpi = 300)

intersect(exc12, exc95)

path_ad <- read.table("95_sig_pathways_ad.txt", header= T, sep = '\t')
path_cn <- read.table("95_sig_pathways_cn.txt", header= T, sep = '\t')
path_ad$Group <- "AD"
path_cn$Group <- "CN"
paths <- rbind(path_ad, path_cn)
exc <- setdiff(path_ad$ID, path_cn$ID)
exc95 <- exc
paths <- paths[paths$ID %in% exc,]
paths <- paths[order(paths$qvalue, decreasing = F),][1:10,]
ad_exc_12 <- ggplot(paths, aes(y = -log10(padjust), x = reorder(ID, -log10(padjust)), color = Module,)) +
  geom_point(stat = "identity", aes(size = Count)) + xlab("Pathway") + ylab("-log10(adj. p-value)") +
  geom_segment(aes(x=ID, xend=ID, y=0, yend=-log10(padjust))) + ylim(0, 5) +
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(text = element_text(size=10, color ="black", face = "bold")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + coord_flip() + facet_wrap(~Group)
ad_exc_12

path_ad <- read.table("95_sig_pathways_ad.txt", header= T, sep = '\t')
path_cn <- read.table("95_sig_pathways_cn.txt", header= T, sep = '\t')
path_ad$Group <- "AD"
path_cn$Group <- "CN"
paths <- rbind(path_ad, path_cn)
exc <- setdiff(path_cn$ID, path_ad$ID)
paths <- paths[paths$ID %in% exc,]
paths <- paths[order(paths$qvalue, decreasing = F),][1:10,]
cn_exc_12 <- ggplot(paths, aes(y = -log10(padjust), x = reorder(ID, -log10(padjust)), color = Module,)) +
  geom_point(stat = "identity", aes(size = Count)) + xlab("Pathway") + ylab("-log10( p-value adjust)") +
  geom_segment(aes(x=ID, xend=ID, y=0, yend=-log10(padjust))) + ylim(0, 5) +
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(text = element_text(size=10, color ="black", face = "bold")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + coord_flip() + facet_wrap(~Group)
cn_exc_12

plot_grid(cn_exc_12, ad_exc_12, nrow = 2)
ggsave("Final_95.pdf", scale = 1, width = 4, height = 8, units = "in", dpi = 300)

install.packages("magick")
install.packages("knitr")
library("knitr")
library("magick")
library("ggplot2")
library("grid")


fig <- image_graph(res = 10)
plot_grid(cn_exc_12, ad_exc_12, nrow = 2)

manual <- image_read_pdf('Final_12.pdf', density = 300)
img = c(cn_exc_12, ad_exc_12, manual)

image_composite(image_scale(img, "x200"))


grid.raster(manual, just = "right")
grid.raster(manual, just = "right")











