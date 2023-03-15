library(plotly)
library(ggrepel)
library(dplyr)
library(tibble)
library(broom)
library(reshape2)
library(tidyr)
library(data.table)
library(here)


#Load input files
RBM14 <- read.csv(file = 'CCLE_RBM14expCN.csv')

PRISM <- read.csv(file = 'PRISMsecondary.csv')


#Extract just relevant info from PRISM dataset and remove duplicates based on screen name
myPRISM <- PRISM[PRISM$passed_str_profiling == TRUE, c("depmap_id", "ccle_name", "screen_id", "auc", "name", "moa", "target")]
head(myPRISM)

# For drugs with duplicate drug/cell line pairs, only use from screen = "MTS005" "MTS006"

duplicates <- duplicated(myPRISM[,c(1,5)])

myPRISMnodup <- myPRISM[!duplicates, ]


#Format RBM14 df
myRBM14 <- RBM14 %>% 
  select(c("depmap_id","cell_line_display_name", "lineage_1", "Expression.Public.22Q4.RBM14", "Copy.Number.Public.22Q4.RBM14")) %>% rename(expression = "Expression.Public.22Q4.RBM14") %>% rename(copy_number = "Copy.Number.Public.22Q4.RBM14")


##Only keep expression for cell lines with auc info

rbm14 <- myRBM14[myRBM14$depmap_id %in% myPRISM$depmap_id ,] 

rownames(rbm14) <-  rbm14$depmap_id

rbm14 <- rbm14[ ,-1]

#Format auc df. Reorder rbm14 df based on cell line order for auc

head(myPRISM)

auc <- myPRISM %>% select(c("depmap_id","auc","name")) %>% pivot_wider(names_from = "name", values_from = "auc") %>% column_to_rownames(var="depmap_id")
auc2 <- apply(auc,2,function(y) sapply(y,function(x) ifelse(is.null(x),0,x)))
auc2 <- as.data.frame(auc2)

reorder_idx <- match(rownames(auc2), rownames(rbm14))
reordered_rbm14 <- rbm14[reorder_idx, ]

head(auc2)
head(reordered_rbm14)

#lm for RBM14 exp
 lm <- apply(auc2, 2, function(x) broom::tidy(lm(x ~ reordered_rbm14$lineage_1 + reordered_rbm14$expression))) %>%
    bind_rows() %>%
    filter(term == "reordered_rbm14$expression") %>%
    mutate(p.adj = p.adjust(p.value, method="fdr"))
 
 head(lm)

 lm <- as.data.frame(lm)
 
rownames(lm) <- colnames(auc2)

lm <- lm %>% rownames_to_column(var="name")

drugmoa <- myPRISM %>% select(c("name","moa"))
drugmoa <- distinct(drugmoa)

lm.moa <- lm %>% inner_join(drugmoa, by="name")

#lm for RBM14 copy number
 lm.cn <- apply(auc2, 2, function(x) broom::tidy(lm(x ~ reordered_rbm14$lineage_1 + reordered_rbm14$copy_number))) %>%
    bind_rows() %>%
    filter(term == "reordered_rbm14$copy_number") %>%
    mutate(p.adj = p.adjust(p.value, method="fdr"))
 

 lm.cn <- as.data.frame(lm.cn)
 
rownames(lm.cn) <- colnames(auc2)

lm.cn <- lm.cn %>% rownames_to_column(var="name")

lm.cn.moa <- lm.cn %>% inner_join(drugmoa, by="name")



# Plot the correlations

PARP <- "PARP inhibitor"

PARP_volcano <- lm.moa[lm.moa$moa %in% PARP, ]
PARP_volcano$PARP <- "PARP inhibitor"

lm.moa$Significant <- ifelse(lm.moa$p.value <= 0.05, "p.val < 0.05", "Not Sig")
lm.cn.moa$Significant <- ifelse(lm.cn.moa$p.value <= 0.05, "p.val < 0.05", "Not Sig")

ggplot(lm.moa, aes(x=estimate, y=-log(p.value))) + xlab("Estimate") + ylab("Significance (-log10 pvalue)") + ggtitle("Correlation between RBM14 expression and drug sensitivity (corrected by lineage)") + geom_point(aes(color = Significant)) + scale_color_manual(values = c("grey", "red","blue")) + theme_classic() + geom_point(data= PARP_volcano, aes(x=estimate,y=-log(p.value), color= PARP,show.legend=FALSE)) +  geom_text_repel(size= 3, color="#2b35e2", data=PARP_volcano, aes(estimate,-log(p.value),label=name))

ggplot(lm.cn.moa, aes(x=estimate, y=-log(p.value))) + xlab("Estimate") + ylab("Significance (-log10 pvalue)") + ggtitle("Correlation between RBM14 CN and drug sensitivity (corrected by lineage)") + geom_point(aes(color = Significant)) + scale_color_manual(values = c("grey", "red","blue")) + theme_classic() + geom_point(data= PARP_volcano, aes(x=estimate,y=-log(p.value), color= PARP,show.legend=FALSE)) +  geom_text_repel(size= 3, color="#2b35e2", data=PARP_volcano, aes(estimate,-log(p.value),label=name))
