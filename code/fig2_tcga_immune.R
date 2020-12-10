## Figure 2a, 2b, 2c, 2d, 2e; figure s4a, b

library(Biobase)
library(BiocGenerics)
library(parallel)

library(ggplot2)
library(pheatmap)
library(ggpubr)
library(ggbeeswarm)

## Read data-----------------
tcga_fit <- read.csv("../data/20191125TCGA_deconvolutionResults.csv", row.names = 1)
tcga_eset <- readRDS("../data/20181029TCGA_eset.rds")

aocs_eset <- readRDS("../data/20181029_tothill_eset.rds")
aocs_dec <- readRDS("../data/20190502Deconvolution_Tothill.rds")

## LM6 results ------
immune2 <- read.delim("../data/CIBERSORT.Output_Job19.txt", row.names = 1)
immune2 <- immune2[match(rownames(tcga_fit), rownames(immune2)),]


## heatmap (figure s4a)
cor.rls2 <- t(cor(tcga_fit[,1:5], immune2[,1:6]))
pheatmap::pheatmap(cor.rls2, 
                   filename = "../results/figS4a.pdf",
                   display_numbers = round(cor.rls2, 2), cellheight = 12, cellwidth = 30)

## violin plot (figure 2a) -----------
df <- cbind(tcga_fit[,1:5], immune2)
df$EMT.group <- "EMT-middle"
df$EMT.group[df$EMT < quantile(df$EMT, 1/3)] <- "EMT-low"
df$EMT.group[df$EMT > quantile(df$EMT, 2/3)] <- "EMT-high"
df$EMT.group <- factor(df$EMT.group, levels = c("EMT-low", "EMT-middle", "EMT-high"))

ggplot(df, aes(x = EMT.group, y =  Monocytes)) + geom_violin(scale = "area") + 
    geom_boxplot(width = 0.2) + geom_beeswarm(alpha = 0.4) + theme_classic2()
ggsave("../results/fig2a.pdf", width = 4, height = 3)

res <- t.test(df$Monocytes[df$EMT.group == "EMT-high"], df$Monocytes[df$EMT.group == "EMT-low"], alternative = "greater")
res_df <- data.frame(p = res$p.value,
                     mean.x = res$estimate[1],
                     mean.y = res$estimate[2],
                     method = res$method,
                     alternative = res$alternative,
                     data = res$data.name)
write.csv(res_df,  "../results/figb2a_ttest_res.csv", row.names = FALSE)



# LM22 deconvolution results -----
immune <- read.delim("../data/CIBERSORT.Output_Job18.txt", row.names = 1)
immune <- immune[match(rownames(tcga_fit), rownames(immune)),]

# macrophage markers used in LM22
LM22 <- read.delim("../data/LM22.txt", as.is = T, row.names = 1)
LM22_macrophage_markers <- names(which(apply(LM22, 1, which.max) == 16))

## correlation heatmap (Figure S4B) -----------------
cor.rls <- t(cor(tcga_fit[,1:5], immune[,1:22]))
pheatmap::pheatmap(cor.rls, 
                   filename = "../results/figS4b.pdf",
                   display_numbers = round(cor.rls, 2), 
                   cellheight = 12, cellwidth = 30)

## Barplot (Figure 2C) -----------------

LM22_macrophage_markers <- unique(c(LM22_macrophage_markers, "CD163","TGM2"))

cor_res <- c()
for(itor_gene in LM22_macrophage_markers){
    if(itor_gene %in% rownames(tcga_eset)){
        
        tmp <- cor(tcga_fit[,"EMT"], log1p(tcga_eset@assayData$exprs[itor_gene,]))
        cor_res <- rbind(cor_res, c(itor_gene, tmp))
    }
}
cor_res <- data.frame(gene = cor_res[,1],
                      cor = as.numeric(cor_res[,2]))
cor_res <- cor_res[order(cor_res$cor, decreasing = T),]
cor_res$sign <- cor_res$cor > 0
cor_res$gene <- factor(cor_res$gene, levels = cor_res$gene)


ggplot(cor_res, aes(x = gene, y = cor, fill = sign)) + geom_bar(stat = "identity") + 
    coord_flip() +
    theme_classic2() + 
    theme(legend.position = "none") + 
    xlab("") + ylab("Correlation with EMT scores")
ggsave("../results/fig2c.pdf", width = 3, height = 4)


## Violin plot (Figure 2B) ------
df <- cbind(tcga_fit[,1:5], immune)
df$EMT.group <- "EMT-middle"
df$EMT.group[df$EMT < quantile(df$EMT, 1/3)] <- "EMT-low"
df$EMT.group[df$EMT > quantile(df$EMT, 2/3)] <- "EMT-high"
df$EMT.group <- factor(df$EMT.group, levels = c("EMT-low", "EMT-middle", "EMT-high"))

ggplot(df, aes(x = EMT.group, y =  Macrophages.M2)) + geom_violin(scale = "area") + 
    geom_boxplot(width = 0.2) + geom_beeswarm(alpha = 0.4) + theme_classic2()
ggsave("../results/fig2b.pdf", width = 4, height = 3)

res <- t.test(df$Macrophages.M2[df$EMT.group == "EMT-high"], df$Macrophages.M2[df$EMT.group == "EMT-low"], alternative = "greater")
res_df <- data.frame(p = res$p.value,
                     mean.x = res$estimate[1],
                     mean.y = res$estimate[2],
                     method = res$method,
                     alternative = res$alternative,
                     data = res$data.name)
write.csv(res_df,  "../results/fig2b_ttest_res.csv", row.names = FALSE)

## AOCS analysis ----------------
aocs_lm22 <- read.csv("../data/CIBERSORT.Output_Job20_AOCS_LM22.csv")
aocs_eset$macrophageM2 <- aocs_lm22$Macrophages.M2

aocs_lm6 <- read.csv("../data/CIBERSORT.Output_Job21_AOCS_LM6.csv")
aocs_eset$monocyte <- aocs_lm6$Monocytes

## violin plot (Figure s2d)-------------
df <- cbind(aocs_dec[,1:5],aocs_lm6[,1:6])
df$EMT.group <- "EMT-middle"
df$EMT.group[df$EMT <= quantile(df$EMT, 1/3)] <- "EMT-low"
df$EMT.group[df$EMT > quantile(df$EMT, 2/3)] <- "EMT-high"
df$EMT.group <- factor(df$EMT.group, levels = c("EMT-low", "EMT-middle", "EMT-high"))

ggplot(df, aes(x = EMT.group, y =  Monocytes)) + geom_violin(scale = "width") + 
    geom_boxplot(width = 0.1) + geom_jitter(alpha = 0.4) + theme_classic2()
ggsave("../results/fig2d.pdf", width = 4, height = 3)

res <- t.test(df$Monocytes[df$EMT.group == "EMT-high"], df$Monocytes[df$EMT.group == "EMT-low"], alternative = "greater")

res_df <- data.frame(p = res$p.value,
                     mean.x = res$estimate[1],
                     mean.y = res$estimate[2],
                     method = res$method,
                     alternative = res$alternative,
                     data = res$data.name)
write.csv(res_df,  "../results/fig2d_ttest_res.csv", row.names = FALSE)


## violin plot (Figure s2e)-------------
df <- cbind(aocs_dec[,1:5],aocs_lm22[,1:23])
df$EMT.group <- "EMT-middle"
df$EMT.group[df$EMT <= quantile(df$EMT, 1/3)] <- "EMT-low"
df$EMT.group[df$EMT > quantile(df$EMT, 2/3)] <- "EMT-high"
df$EMT.group <- factor(df$EMT.group, levels = c("EMT-low", "EMT-middle", "EMT-high"))

ggplot(df, aes(x = EMT.group, y =  Macrophages.M2)) + geom_violin(scale = "width") + 
    geom_boxplot(width = 0.1) + geom_jitter(alpha = 0.4) + theme_classic2()
ggsave("../results/fig2e.pdf", width = 4, height = 3)

res <- t.test(df$Macrophages.M2[df$EMT.group == "EMT-high"], df$Macrophages.M2[df$EMT.group == "EMT-low"], alternative = "greater")

res_df <- data.frame(p = res$p.value,
                     mean.x = res$estimate[1],
                     mean.y = res$estimate[2],
                     method = res$method,
                     alternative = res$alternative,
                     data = res$data.name)
write.csv(res_df,  "../results/fig2e_ttest_res.csv", row.names = FALSE)
