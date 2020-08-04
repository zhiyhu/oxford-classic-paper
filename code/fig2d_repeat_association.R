library(ggplot2)
library(ggpubr)
# Figure 2D
# This analysis is aimed to verify that the association between the macrophage M2 signature and the mesenchymal signature is reproducible in multiple independent datasets. See Figure 2D 

## ----read data--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
expr_list <- readRDS("../data/Ovarian_expression_data_SPARC_CD163_list.rds")


## ----cor-test---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df <- data.frame(dataset = c("TCGA","AOCS","E.MTAB.386", "GSE13876", "GSE26193",
              "GSE26712", "GSE32062.GPL6480",
              "GSE49997", "GSE51088"),
              p = NA,
              r = NA)
for(i in 1:9){
    rls <- cor.test(expr_list[[i]]["CD163",], expr_list[[i]]["SPARC",])
    df$p[i] <- rls$p.value
    df$r[i] <- rls$estimate
}

write.csv(df, "../results/fig2d_pearson_correlation.csv", row.names = FALSE)


## ----prepare-data-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tmp <- log1p(expr_list[[1]][c("CD163","SPARC"),])
for(i in 2:9){
    tmp <- cbind(tmp, expr_list[[i]][c("CD163","SPARC"),])
}
tmp <- t(tmp)

id <- c( rep(1, ncol(expr_list[[1]])),
         rep(2, ncol(expr_list[[2]])),
         rep(3, ncol(expr_list[[3]])),
         rep(4, ncol(expr_list[[4]])),
         rep(5, ncol(expr_list[[5]])),
         rep(6, ncol(expr_list[[6]])),
         rep(7, ncol(expr_list[[7]])),
         rep(8, ncol(expr_list[[8]])),
         rep(9, ncol(expr_list[[9]]))
         )
id <- id[!(is.na(tmp[,1]) | is.na(tmp[,2]))]
tmp <- tmp[!(is.na(tmp[,1]) | is.na(tmp[,2])),]


## ----scatter plot-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_plot <- data.frame(tmp, id)
df_id <- data.frame(id = 1:9,
                    name = c("TCGA","AOCS","E.MTAB.386","GSE13876","GSE26193","GSE26712","GSE34062.GPL6480","GSE49997","GSE51088")
)
df_plot$name <- df_id$name[match(df_plot$id, df_id$id)]
df_plot$name <- factor(df_plot$name, levels = df_id$name)
ggplot(df_plot, aes(x = SPARC, y = CD163)) +
  geom_smooth(method="lm", alpha = 0.5) + 
    geom_point(size=1, pch = 21, col = "grey40", fill = alpha("grey60", 0.5)) + 
    facet_wrap(~name, scales = "free")  + theme_classic2() + xlab("Expression levels of SPARC") + 
    ylab("Expression levels of CD163") 
    
ggsave("../results/fig2d.pdf", width = 6, height = 5)


