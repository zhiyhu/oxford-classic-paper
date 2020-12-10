# fig S3A, B and C; tables S3, 4

library(survminer)
library(survival)
library(ggplot2)
library(ggbeeswarm)

# Read TCGA data ----
tcga_fit <- read.csv("../data/20191125TCGA_deconvolutionResults.csv", row.names = 1)
tcga_fit <- data.frame(tcga_fit)
tcga_eset <- readRDS("../data/20181029TCGA_eset.rds")

# BRCA status ---- 
## BRCA1
brca_data <- read.delim("../data/PATIENT_DATA_oncoprint_BRCA1.tsv",as.is = TRUE)
colnames(brca_data) <- gsub("[.]", "-", colnames(brca_data))
df_brca <- t(brca_data[,-1:-2])
colnames(df_brca) <- paste0(brca_data$track_name, "_", brca_data$track_type)
df_brca <- data.frame(df_brca)
df_brca <- df_brca[match(substr(colnames(tcga_eset), 1, 12), rownames(df_brca)),]
tcga_eset$Profiled.in.Mutations_BRCA1 <- df_brca$Profiled.in.Mutations_CLINICAL
tcga_eset$BRCA1_CNA <- df_brca$BRCA1_CNA
tcga_eset$BRCA1_MUTATIONS <- df_brca$BRCA1_MUTATIONS

## BRCA2
brca_data <- read.delim("../data/PATIENT_DATA_oncoprint_BRCA2.tsv",as.is = TRUE)
colnames(brca_data) <- gsub("[.]", "-", colnames(brca_data))

df_brca <- t(brca_data[,-1:-2])
colnames(df_brca) <- paste0(brca_data$track_name, "_", brca_data$track_type)
df_brca <- df_brca[rownames(df_brca) %in% substr(colnames(tcga_eset), 1, 12),]
df_brca <- data.frame(df_brca)

df_brca <- df_brca[match(substr(colnames(tcga_eset), 1, 12), rownames(df_brca)),]
tcga_eset$Profiled.in.Mutations_BRCA2 <- df_brca$Profiled.in.Mutations_CLINICAL
tcga_eset$BRCA2_CNA <- df_brca$BRCA2_CNA
tcga_eset$BRCA2_MUTATIONS <- df_brca$BRCA2_MUTATIONS

tcga_eset$BRCA_status <- NA

df_brca <- tcga_eset@phenoData@data[,-1:-110]
idx <- is.na(df_brca$BRCA1_CNA) & is.na(df_brca$BRCA2_CNA)
df_brca$BRCA_status[idx] <- "data_not_available"
idx <- which(df_brca$Profiled.in.Mutations_BRCA1 == "Yes" & 
                 df_brca$Profiled.in.Mutations_BRCA2 == "Yes" &
                 df_brca$BRCA1_CNA == "" &
                 df_brca$BRCA2_CNA == "" &
                 df_brca$BRCA1_MUTATIONS == "" &
                 df_brca$BRCA2_MUTATIONS == "" 
)
df_brca$BRCA_status[idx] <- "wildtype"

idx <- which(grepl("putative driver",df_brca$BRCA1_MUTATIONS)|
                 grepl("putative driver",df_brca$BRCA2_MUTATIONS))
df_brca$BRCA_status[idx] <- "mutated"

idx <- which(grepl("Deep Deletion",df_brca$BRCA1_CNA)|
                 grepl("Deep Deletion",df_brca$BRCA2_CNA))
df_brca$BRCA_status[idx] <- "mutated"
df_brca$BRCA_status[is.na(df_brca$BRCA_status)] <- "wildtype"

tcga_eset$brca_mutated <- df_brca$BRCA_status
metadata <- tcga_eset@phenoData@data
metadata <- cbind(metadata, tcga_fit)
metadata <- metadata[metadata$brca_mutated != "data_not_available",]

# Correlation test (fig S3A) ----
idx_mu <- metadata$brca_mutated == "mutated"
idx_wt <- metadata$brca_mutated == "wildtype"

df_ttest <- data.frame(sig = c("Differentiated", "EMT", "KRT17", "Cell_cycle", "Ciliated"), p_val = NA)
df_ttest$p_val[df_ttest$sig == "Differentiated"] <- t.test(metadata$Differentiated[idx_mu],
                                                           metadata$Differentiated[idx_wt])$p.value
df_ttest$p_val[df_ttest$sig == "EMT"] <- t.test(metadata$EMT[idx_mu],
                                                metadata$EMT[idx_wt])$p.value
df_ttest$p_val[df_ttest$sig == "KRT17"] <- t.test(metadata$KRT17[idx_mu],
                                                  metadata$KRT17[idx_wt])$p.value
df_ttest$p_val[df_ttest$sig == "Cell_cycle"] <- t.test(metadata$Cell_cycle[idx_mu],
                                                       metadata$Cell_cycle[idx_wt])$p.value
df_ttest$p_val[df_ttest$sig == "Ciliated"] <- t.test(metadata$Ciliated[idx_mu],
                                                     metadata$Ciliated[idx_wt])$p.value
write.csv(df_ttest, "../results/figS3a_pvalues.csv", row.names = FALSE)

# Violin plots (fig S3A) ----
p1 <- ggplot(metadata, aes(x = brca_mutated, y = EMT)) + geom_violin() + geom_boxplot(width = 0.1) + geom_jitter(alpha  = 0.5, width = 0.3) + theme_classic()
p2 <- ggplot(metadata, aes(x = brca_mutated, y = Differentiated)) + geom_violin() + geom_boxplot(width = 0.1) + geom_jitter(alpha  = 0.5, width = 0.3) + theme_classic()
p3 <- ggplot(metadata, aes(x = brca_mutated, y = KRT17)) + geom_violin() + geom_boxplot(width = 0.1) + geom_jitter(alpha  = 0.5, width = 0.3) + theme_classic()
p4 <- ggplot(metadata, aes(x = brca_mutated, y = Cell_cycle)) + geom_violin() + geom_boxplot(width = 0.1) + geom_jitter(alpha  = 0.5, width = 0.3) + theme_classic()
p5 <- ggplot(metadata, aes(x = brca_mutated, y = Ciliated)) + geom_violin() + geom_boxplot(width = 0.1) + geom_jitter(alpha  = 0.5, width = 0.3) + theme_classic()
cowplot::plot_grid(p1,p2,p3,p4,p5, ncol = 5)
ggsave("../results/figS3a.pdf", useDingbats = FALSE, width = 15, height = 3)

# KM-curves (fig S3B) -----
tcga_eset$EMT <- tcga_fit$EMT
tcga_eset$stage_summary <- NA
early_stages <- c("Stage IA","Stage IB","Stage IC","Stage IIA","Stage IIB","Stage IIC")
late_stages  <- c(c("Stage IIIA","Stage IIIB","Stage IIIC","Stage IV"))
tcga_eset$stage_summary[tcga_eset$clinical_stage %in% early_stages] <- "early"
tcga_eset$stage_summary[tcga_eset$clinical_stage %in% late_stages] <- "late"
tcga_eset$EMT.group[tcga_eset$EMT < quantile(tcga_eset$EMT, 1/2)] <- "EMT-low"
tcga_eset$EMT.group[tcga_eset$EMT > quantile(tcga_eset$EMT, 1/2)] <- "EMT-high"

df2 <- pData(tcga_eset)[tcga_eset$brca_mutated != "data_not_available",]
df2$brca_mutated <- as.character(df2$brca_mutated)
df2$EMT.group <- as.character(df2$EMT.group)
df2$group <- paste0(df2$EMT.group, "_", df2$brca_mutated)
fit1 <- survfit(Surv(OS.time, OS) ~ group, df2)
ggsurvplot(fit1, df2, pval = FALSE, risk.table = FALSE)
ggsave("../results/figS3b.pdf", width = 8, height = 6)

# Survival analysis (table s3) ----
fit <- coxph(Surv(OS.time, OS) ~ EMT + stage_summary + brca_mutated, data=pData(tcga_eset)[tcga_eset$brca_mutated != "data_not_available",])
df <- data.frame(sample.size = summary(fit)$n,
                 hr = summary(fit)$coefficients[,2],
                 p =  summary(fit)$coefficients[,5],
                 ci.lower = summary(fit)$conf.int[,3],
                 ci.upper = summary(fit)$conf.int[,4])
write.csv(df, "../results/tableS3.csv")

# Read TCGA microarray data ----
cibersort_fit133 <- readRDS("../data/TCGA_U133A_deconvolution.rds")
metadata <- readRDS("../data/TCGA_U133A_clinical_data.rds")
metadata <- metadata[which(metadata$brca_mutated != "data_not_available"),]

# Survival analysis (table s4) ----
fit <- coxph(Surv(OS.time, OS)~ EMT + summary_stage + brca_mutated, data = metadata)
df <- data.frame(sample.size = summary(fit)$n,
                 hr = summary(fit)$coefficients[,2],
                 p =  summary(fit)$coefficients[,5],
                 ci.lower = summary(fit)$conf.int[,3],
                 ci.upper = summary(fit)$conf.int[,4])
write.csv(df, "../results/tables4.csv")

# KM curve (fig s3c)----
metadata$EMT.group[metadata$EMT <= quantile(metadata$EMT, 1/2)] <- "EMT-low"
metadata$EMT.group[metadata$EMT > quantile(metadata$EMT, 1/2)] <- "EMT-high"
metadata$group <- paste0(metadata$EMT.group, "_", metadata$brca_mutated)
fit1 <- survfit(Surv(OS.time, OS) ~ group, metadata)
ggsurvplot(fit1, metadata, pval = FALSE, risk.table = FALSE)
ggsave("../results/figS3c.pdf", width = 8, height = 6)
