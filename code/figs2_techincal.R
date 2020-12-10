## Figure S2
# Correlation analysis between RNA-seq data and NanoString readout
# Zhiyuan Hu
# 9 Oct 2020; last modify 10 Dec 2020

library(ggplot2)
library(cowplot)

# read nanostring data
ns_data <- read.csv("../data/OXOPCR_normalised_data.csv", as.is = TRUE, row.names = 1)
colnames(ns_data) <- gsub("X", "", colnames(ns_data))
colnames(ns_data) <- gsub("[.]", "", colnames(ns_data))
dim(ns_data)

# read RNA-Seq data
counts <- read.delim("../data/oxopcr_pre_sections_counts.txt")
counts <- counts[match(rownames(ns_data), rownames(counts)), ]
colnames(counts) <- sapply(colnames(counts), function(x) return(unlist(strsplit(x, "_"))[[2]]))
counts2 <- counts[match(rownames(ns_data), rownames(counts)),]

## correlation by samples
plist <- list()
for (i in 1:10) {
    df_plot <- data.frame(nanostring_readout = log2(ns_data[,i] + 1),
                          RNASeq_readout = log2(counts2[,i] + 1))
    pearsonr <- round(cor(df_plot[,1], df_plot[,2]), 2)
    pvalue <- signif(cor.test(df_plot[,1], df_plot[,2])$p.value, digits = 3)
    plist[[i]] <- ggplot(df_plot, aes(x = nanostring_readout, y = RNASeq_readout)) + geom_point() + 
        theme_linedraw()  + annotate("text", x = 5, y = 15, label = paste0("r = ", pearsonr, "\np = ", pvalue)) + xlim(0, 17) + ylim(0, 17) +
        xlab("log2(NanoString normalised count + 1)") + ylab("log2(RNA-Seq read count + 1)") + labs(title = colnames(ns_data)[i])
}
plot_grid(plotlist = plist, ncol = 5)
ggsave("../results/figS2.pdf", width = 15, height = 6, useDingbats = FALSE)
