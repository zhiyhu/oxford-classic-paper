library(pROC)
library(ggplot2)
library(caret)
library(lattice)

df_caps <- read.csv("../data/CAPS TMA 19.7.2019.csv", as.is = T)
df_caps$TMA[df_caps$TMA == "TA 172"] <- "CAPS 172"
colnames(df_caps) <- gsub(pattern = "[..]", replacement = "_", x = colnames(df_caps))
colnames(df_caps) <- gsub(pattern = "__", replacement = "_", x = colnames(df_caps))

df_caps$ratio_tumour_1 <- df_caps$Num_Tumor_1_/df_caps$Num_Tumor
df_caps$ratio_tumour_2 <- df_caps$Num_Tumor_2_/df_caps$Num_Tumor
df_caps$ratio_tumour_3 <- df_caps$Num_Tumor_3_/df_caps$Num_Tumor

# filter out some TMAs
df_caps <- df_caps[df_caps$TMA %in% c("CAPS 170", "CAPS 171", "CAPS 172", "CAPS 173", "CAPS 178", "CAPS 183"),]
table(df_caps$tumor_type)

set.seed(34521)
idx <- createDataPartition(df_caps$tumor_type, times = 2, p = 0.6)

training_dt <- df_caps[idx$Resample1,]
testing_dt <- df_caps[-idx$Resample1,]

## Logistic regression
model_glm = glm(tumor_type ~ ratio_tumour_1 + ratio_tumour_2 + ratio_tumour_3, data = training_dt, family = "binomial")
model_glm_pred = ifelse(predict(model_glm, type = "link") > 0, 1, 0)

train_tab = table(predicted = model_glm_pred, actual = training_dt$tumor_type)

train_con_mat = caret::confusionMatrix(train_tab)
c(train_con_mat$overall["Accuracy"],
  train_con_mat$byClass["Sensitivity"],
  train_con_mat$byClass["Specificity"])

# ROC (figure 3c) ------
test_prob = predict(model_glm, newdata = testing_dt, type = "response")
pdf("../results/fig3c.pdf")
roc(testing_dt$tumor_type, test_prob, plot = TRUE, print.auc = TRUE, legacy.axes=TRUE)
dev.off()

## distribution of AUROC (figure 3d)

set.seed(34521)
idx <- createFolds(df_caps$tumor_type, k = 8)

test_roc_list <- list()
model_list <- list()

for(itor in 1:length(idx)) {
    
    training_dt <- df_caps[-idx[[itor]],]
    testing_dt <- df_caps[idx[[itor]],]
    
    ## Logistic regression
    model_glm = glm(tumor_type ~ ratio_tumour_1 + ratio_tumour_2 + ratio_tumour_3, data = training_dt, family = "binomial")
    model_list[[itor]] <- model_glm
    
    # ROC
    test_prob <- predict(model_glm, newdata = testing_dt, type = "response")
    test_roc_list[[itor]] <- roc(testing_dt$tumor_type ~ test_prob, plot = F, print.auc = F)
}


auc_list <- data.frame(value = sapply(test_roc_list, function(x) return(x$auc)))
p <- ggplot(auc_list, aes(x = value))
# add command to produce a "faded" histogram and select the number of bins
# then overlay density curve (converted to common axis of count)
line_col <- RColorBrewer::brewer.pal(3, "Paired")
p  + 
    geom_histogram(aes(y=..density..),colour = "white", fill = line_col[1], binwidth = 0.05, alpha = 0.7) + geom_density(col = line_col[2]) + xlim(0,1)  + theme_classic() + xlab("AUC") + ylab("Frequency")

ggsave("../results/fig3d.pdf", width = 5, height = 3)

