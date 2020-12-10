## ----setup, include=FALSE, cache=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------
library(survival)
library(viridis)
library(ggplot2)
library(pheatmap)


## ----read data--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## signature matrix
sig_matrix_new <-  read.csv('../data/TableS7_selected_signature_matrix_52genes20190204.csv', as.is = T, row.names = 1)
sig_matrix_new <- as.matrix(sig_matrix_new)

## clincial data
clic_dt <- read.csv('../data/clinical_data.csv', as.is = T) 
colnames(clic_dt) <- gsub(pattern = '[.]', replacement = '_', x = colnames(clic_dt))
nanostring_fit <- read.csv(file = '../data/cibersort_HHcohort150_output.csv', row.names = 1)
colnames(nanostring_fit)[c(1,2,4,5)] <- c('Differentiated','KRT17','Cell_cycle','Ciliated')
clic_dt <- cbind(clic_dt, nanostring_fit[match(clic_dt$Anonymised_code, rownames(nanostring_fit)),])

## Nanostring data
counts <- read.csv('../data/HH001-150_NormalizedData_cleaned.csv', as.is = T, row.names = 1)
rownames(counts)[7] <- 'MLF1IP'
colnames(counts) <- gsub(pattern = 'HH', replacement = 'HH-', x = colnames(counts))
colnames(counts) <- gsub(pattern = '[.]', replacement = '', x = colnames(counts))


## ----read-clinical dat------------------------------------------------------------------------------------------------------------------------------------------------------------------------
clic_filtered <- read.csv('../data/clinical_data_filtered.csv', as.is = TRUE)


## ----stacked-barplot, fig.height=3, fig.width=10----------------------------------------------------------------------------------------------------------------------------------------------
svr_fit_oxo <- clic_dt[,c('Differentiated','KRT17','EMT','Cell_cycle','Ciliated')]
svr_fit_new_diff <- svr_fit_oxo[svr_fit_oxo[,1] > 0.5,]
svr_fit_new_diff <- svr_fit_new_diff[order(svr_fit_new_diff[,1], decreasing = F),]
svr_fit_new_krt <- svr_fit_oxo[svr_fit_oxo[,2] > 0.5,]
svr_fit_new_krt <- svr_fit_new_krt[order(svr_fit_new_krt[,2], decreasing = F),]
svr_fit_new_emt <- svr_fit_oxo[svr_fit_oxo[,3] > 0.5,]
svr_fit_new_emt <- svr_fit_new_emt[order(svr_fit_new_emt[,3], decreasing = F),]
svr_fit_new_cc <- svr_fit_oxo[svr_fit_oxo[,4] > 0.5,]
svr_fit_new_cc <- svr_fit_new_cc[order(svr_fit_new_cc[,4], decreasing = F),]

svr_fit_new_others <- svr_fit_oxo[!rownames(svr_fit_oxo) %in% c(rownames(svr_fit_new_diff),rownames(svr_fit_new_krt),
                                                                rownames(svr_fit_new_emt),rownames(svr_fit_new_cc)),]

plot.data <- data.frame(sample = rep(rownames(svr_fit_oxo), 5),
                        state = rep(c('Differentiated','KRT17','EMT','Cell cycle','Ciliated'), 
                                    each = nrow(svr_fit_oxo)),
                        score = c(svr_fit_oxo[,1],
                                  svr_fit_oxo[,2],
                                  svr_fit_oxo[,3],
                                  svr_fit_oxo[,4],
                                  svr_fit_oxo[,5]))
plot.data$sample <- factor(plot.data$sample, levels = c(rownames(svr_fit_new_diff),rownames(svr_fit_new_krt),
                                                        rownames(svr_fit_new_emt),rownames(svr_fit_new_cc),
                                                        rownames(svr_fit_new_others)))
plot.data$state <- factor(plot.data$state, levels = c('Differentiated','KRT17','EMT','Cell cycle','Ciliated'))
# differentiated #A16BA3
# F6A000 KRT7
# FFDE00 EMT
# CC #4DD9FF
#F3A1B6 Ciliated


ggplot(plot.data) + geom_bar(aes(y = score, x = sample, fill = state),
                             stat='identity') + theme_classic() + 
    scale_fill_manual(values = c('#A16BA3','#F6A000','#FFDE00','#4DD9FF','#F3A1B6'), 
                      breaks = c('Differentiated','KRT17','EMT','Cell cycle','Ciliated'), 
                      labels = c('Differentiated','KRT17','EMT','Cell cycle','Ciliated')) +
    xlab('Bulk tumour samples \n (each column is one tumour case)') + ylab('Proportions of cell states') +
    theme(axis.ticks.x  = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size= 15),
          legend.text = element_text(size = 12),
          legend.title  = element_blank(),
          axis.line = element_line(size = 0.8),
          legend.position = 'top') 
ggsave('../results/fig1d.pdf', height = 3, width = 10)


## ----heatmap----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_plot <- counts
df_plot <- df_plot[rownames(sig_matrix_new), levels(plot.data$sample)]
df_plot <- t ( scale(t(df_plot + .5), center = T, scale = T))
df_plot[df_plot >= 5 ] <- 5
df_plot[df_plot <  -2 ] <- -2

pheatmap(df_plot,color = inferno(50),
         filename = '../results/fig1c.pdf',
                   show_colnames = FALSE,
                   scale = 'none', cluster_cols = F, 
                   cluster_rows = F)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

res_tb <- data.frame(Variate = c('EMT','Stage','Grade','Residual_disease'),
                     n            = NA,
                     Hazard_ratio = NA,
                     p            = NA,
                     CI_low       = NA,
                     CI_high      = NA)

surv_fit <- coxph(Surv(OS_time_new, OS_event_new)~EMT, data = as.data.frame(clic_filtered))
i <- 1
res_tb$n[i] <- surv_fit$n
res_tb$Hazard_ratio[i] <- summary(surv_fit)$coefficients[,2]
res_tb$p[i] <- summary(surv_fit)$coefficients[,5]
res_tb$CI_low[i] <- summary(surv_fit)$conf.int[,3]
res_tb$CI_high[i] <- summary(surv_fit)$conf.int[,4]

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
surv_fit <- coxph(Surv(OS_time_new, OS_event_new)~Summary_stage, data=as.data.frame(clic_filtered))
i <- 2
res_tb$n[i] <- surv_fit$n
res_tb$Hazard_ratio[i] <- summary(surv_fit)$coefficients[,2]
res_tb$p[i] <- summary(surv_fit)$coefficients[,5]
res_tb$CI_low[i] <- summary(surv_fit)$conf.int[,3]
res_tb$CI_high[i] <- summary(surv_fit)$conf.int[,4]

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
surv_fit <- coxph(Surv(OS_time_new, OS_event_new)~Summary_grade, data=as.data.frame(clic_filtered))
i <- 3
res_tb$n[i] <- surv_fit$n
res_tb$Hazard_ratio[i] <- summary(surv_fit)$coefficients[,2]
res_tb$p[i] <- summary(surv_fit)$coefficients[,5]
res_tb$CI_low[i] <- summary(surv_fit)$conf.int[,3]
res_tb$CI_high[i] <- summary(surv_fit)$conf.int[,4]


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
surv_fit <- coxph(Surv(OS_time_new, OS_event_new)~Residual_disease, data=as.data.frame(clic_filtered))
i <- 4
res_tb$n[i] <- surv_fit$n
res_tb$Hazard_ratio[i] <- summary(surv_fit)$coefficients[,2]
res_tb$p[i] <- summary(surv_fit)$coefficients[,5]
res_tb$CI_low[i] <- summary(surv_fit)$conf.int[,3]
res_tb$CI_high[i] <- summary(surv_fit)$conf.int[,4]

write.csv(res_tb, '../results/table2_part1.csv', row.names = F)


## ----multivariate-OS--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
surv_fit <- coxph(Surv(OS_time_new, OS_event_new)~EMT + Summary_stage + Summary_grade + Residual_disease, data = as.data.frame(clic_filtered))


## ----multivariate-res-table-------------------------------------------------------------------------------------------------------------------------------------------------------------------
res_tb <- data.frame(Variate = c('EMT','Stage','Grade','Residual_disease'),
                     Hazard_ratio = round(summary(surv_fit)$coefficients[,2],digits = 3),
                     p            = round(summary(surv_fit)$coefficients[,5],digits = 3),
                     CI_low       = round(summary(surv_fit)$conf.int[,3]    ,digits = 3),
                     CI_high      = round(summary(surv_fit)$conf.int[,4]    ,digits = 3) )

write.csv(res_tb, '../results/table2_part2.csv', row.names = F)


## ----PFS - multivariate analysis--------------------------------------------------------------------------------------------------------------------------------------------------------------
surv_fit <- coxph(Surv(PFS_time, PFS_event)~EMT+Summary_stage+Summary_grade+Residual_disease,
                  data=as.data.frame(clic_filtered))


## ----multivariate-res-table-PFS---------------------------------------------------------------------------------------------------------------------------------------------------------------
res_tb <- data.frame(Variate = c('EMT','Stage','Grade','Residual_disease'),
                     Hazard_ratio = round(summary(surv_fit)$coefficients[,2],digits = 3),
                     p            = round(summary(surv_fit)$coefficients[,5],digits = 3),
                     CI_low       = round(summary(surv_fit)$conf.int[,3]    ,digits = 3),
                     CI_high      = round(summary(surv_fit)$conf.int[,4]    ,digits = 3) )

write.csv(res_tb, '../results/tableS1.csv', row.names = F)

