## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggpubr)

## ----Low-grade-ciliated, fig.height=4, fig.width=3--------------------------------------------------------------------------------------------------------------------------------------------
clic_filtered <- read.csv("../data/clinical_data_filtered.csv", as.is = TRUE)
ggplot(clic_filtered[!is.na(clic_filtered$Summary_grade),], aes(x = Summary_grade, y = Ciliated)) + 
    geom_violin(scale = "width") + geom_boxplot(width = 0.1) + geom_jitter(alpha = 0.3) + 
    ylab("Ciliated scores") + xlab("Grade") + theme_classic2()

ggsave("../results/fig3a.pdf", width = 3, height = 3)


rls <- t.test(clic_filtered$Ciliated[clic_filtered$Summary_grade == "low"], clic_filtered$Ciliated[clic_filtered$Summary_grade == "high"], alternative = "greater")

res_df <- data.frame(type = c("p", "mean.x", "mean.y"),
                     value = c(rls$p.value, rls$estimate))

write.csv(res_df, "../results/fig3a_ttest_res.csv", row.names = FALSE)