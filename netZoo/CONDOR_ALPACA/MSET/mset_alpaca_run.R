## MSET FOR EVERYTHING ##

options(stringsAsFactors = FALSE)
library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)

# Create cell types list:
# cell_types <- readxl::read_xlsx("~/Documents/projects/ms-victor/sex_differences_fetal_brain/netZoo/alpaca/_m/MSET/databases/cell_type(polioudakis2019)/1-s2.0-S0896627319305616-mmc5.xlsx",
#                                 sheet = 3)
# names_list <- unique(cell_types$Cluster)
# 
# for (i in names_list){
#   a <- cell_types[cell_types$Cluster == i,2]
#   write.table(a, file = paste0("~/Documents/projects/ms-victor/sex_differences_fetal_brain/netZoo/alpaca/_m/MSET/databases/cell_type(polioudakis2019)/genes_",i,".txt"),
#               sep = "", col.names = F, row.names = F, quote = F)
# }


### Input and databases must be selected by hand and parameters fed to the console prompt.


################################# MALE SPECIFIC
# Build the lists with genes and background for each cell type:
setwd("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca")
bg = as.data.frame(fread("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/male_modules_just_genes.tsv")) # 11151 elements

# Build the input for male modules:
bg = bg[bg$gene_class == "Target",] #10421 genes

for (i in unique(bg$module)){
  input = as.data.frame(bg[bg$module == i,1])
  a = as.data.frame(bg[!bg$module == i,1])
  colnames(input) = "genes"
  colnames(a) = "genes"
  print(paste(i,length(input$genes)),sep="_")
  input = rbind(input,a)
  write.table(input, file = paste0("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/MSET/module_",i,"_gene_list.txt"),
              sep = "", col.names = F, row.names = F, quote = F)
  
}

# Reference values:
# [1] "1 2690"
# [1] "2 2466"
# [1] "3 1010"
# [1] "4 1172"
# [1] "5 1680"
# [1] "6 579"
# [1] "7 261"
# [1] "8 94"
# [1] "9 717"
# [1] "10 482"

################# MALE CELL TYPE:
# Run the MSET for module 1 cell type:
setwd("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/MSET")
sink("MSET_results_mod1_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 2690 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 2 cell type:
sink("MSET_results_mod2_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 2466 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 3 cell type:
sink("MSET_results_mod3_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1010 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 4 cell type:
sink("MSET_results_mod4_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1172 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 5 cell type:
sink("MSET_results_mod5_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1680 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 6 cell type:
sink("MSET_results_mod6_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 579 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 7 cell type:
sink("MSET_results_mod7_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 261 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 8 cell type:
sink("MSET_results_mod8_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 94 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 9 cell type:
sink("MSET_results_mod9_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 717 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 10 cell type:
sink("MSET_results_mod10_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 482 (number of genes on the list)
# 2: 10000 (number of permutations)


################# MALE DISORDERS:
# Run the MSET for module 1 disorders:
setwd("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/MSET")
sink("MSET_results_mod1_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 2690 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 2 disorders:
sink("MSET_results_mod2_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 2466 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 3 disorders:
sink("MSET_results_mod3_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1010 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 4 disorders:
sink("MSET_results_mod4_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1172 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 5 disorders:
sink("MSET_results_mod5_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1680 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 6 disorders:
sink("MSET_results_mod6_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 579 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 7 disorders:
sink("MSET_results_mod7_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 261 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 8 disorders:
sink("MSET_results_mod8_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 94 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 9 disorders:
sink("MSET_results_mod9_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 717 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 10 disorders:
sink("MSET_results_mod10_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 482 (number of genes on the list)
# 2: 10000 (number of permutations)

################################# FEMALE SPECIFIC
# Build the lists with genes and background for each cell type:
bg = as.data.frame(fread("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/female_modules_genes.tsv")) # 11147 genes

# Build the input for male modules:
bg = bg[bg$gene_class == "Target",] #10417

for (i in unique(bg$module)){
  input = as.data.frame(bg[bg$module == i,1])
  a = as.data.frame(bg[!bg$module == i,1])
  colnames(input) = "genes"
  colnames(a) = "genes"
  print(paste(i,length(input$genes)),sep="_")
  input = rbind(input,a)
  write.table(input, file = paste0("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/MSET/female_module_",i,"_gene_list.txt"),
              sep = "", col.names = F, row.names = F, quote = F)
  
}

# Reference values:
# [1] "1 2399"
# [1] "2 1900"
# [1] "3 1472"
# [1] "4 1742"
# [1] "5 2007"
# [1] "6 478"
# [1] "7 419"

################# FEMALE CELL TYPE:
# Run the MSET for module 1 cell type:
sink("MSET_results_fmod1_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 2399 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 2 cell type:
sink("MSET_results_fmod2_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1900 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 3 cell type:
sink("MSET_results_fmod3_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1472 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 4 cell type:
sink("MSET_results_fmod4_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1742 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 5 cell type:
sink("MSET_results_fmod5_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 2007 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 6 cell type:
sink("MSET_results_fmod6_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 478 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 7 cell type:
sink("MSET_results_fmod7_cell_type", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 419 (number of genes on the list)
# 2: 10000 (number of permutations)


################# FEMALE DISORDERS:
# Run the MSET for module 1 disorders:
sink("MSET_results_fmod1_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 2399 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 2 disorders:
sink("MSET_results_fmod2_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1900 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 3 disorders:
sink("MSET_results_fmod3_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1472 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 4 disorders:
sink("MSET_results_fmod4_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 1742 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 5 disorders:
sink("MSET_results_fmod5_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 2007 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 6 disorders:
sink("MSET_results_fmod6_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 478 (number of genes on the list)
# 2: 10000 (number of permutations)

# Run the MSET for module 7 disorders:
sink("MSET_results_fmod7_disorders", append=F, split=F)
source("mset_v2.R")
sink()

# 1: 419 (number of genes on the list)
# 2: 10000 (number of permutations)

# Data was manually compiled in Google Sheets and downloaded as .tsv to generate the balloon plots below


################################## Cell type:
results_cell_type = read.csv("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/MSET_results_modules - cell_type.csv", header = TRUE)
results_female <- results_cell_type[results_cell_type$sex == "Female",]
results_male <- results_cell_type[results_cell_type$sex == "Male",]

fresults <- data.frame()
for (i in unique(results_female$Tissue)){
  a <- results_female[results_female$Tissue == i,]
  a$adj.p.valor <- p.adjust(a$p.valor, method = "BH")
  fresults <- rbind(fresults,a)
}

mresults <- data.frame()
for (i in unique(results_male$Tissue)){
  a <- results_male[results_male$Tissue == i,]
  a$adj.p.valor <- p.adjust(a$p.valor, method = "BH")
  mresults <- rbind(mresults,a)
}

### Female
# picking just the number of matches and renaming the columns
fresults$module <- paste0("M",fresults$module)
fresults_match = fresults[,c(1,5,11,8)]
colnames(fresults_match) = c("Tissues","Fold enrichment", "Adjusted p-value", "Female baseline modules")
fresults_match %>% 
  mutate(zscore = (`Fold enrichment` - mean(`Fold enrichment`))/sd(`Fold enrichment`)) -> fresults_match
colnames(fresults_match)[5] <- "Z-score"

# Plot:
library(ggplot2)
par(mar=(c(1,1,1,1)))
pdf("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/female_mset_tissue.pdf",
    width = 6, height = 5)
p = ggplot(fresults_match[which(fresults_match$`Fold enrichment`>0),], aes(x=`Female baseline modules`,y=`Tissues`)) +
  geom_point(aes(size = `Fold enrichment`, colour = `Adjusted p-value`), shape = 16) +
  theme_bw() +
  ggtitle("Midgestation brain cell types enrichment") +
  theme(plot.title = element_text(hjust = 0.65), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_continuous(limits=c(min(fresults_match$`Adjusted p-value`),0.05), high = "#F9A499", low = "#C13624", na.value = "grey80") +
  scale_y_discrete(limits=rev)
p
dev.off()


### Male
# picking just the number of matches and renaming the columns
mresults$module <- paste0("M",mresults$module)
mresults_match = mresults[,c(1,5,11,8)]
colnames(mresults_match) = c("Tissues","Fold enrichment", "Adjusted p-value", "Male baseline modules")
mresults_match %>% 
  mutate(zscore = (`Fold enrichment` - mean(`Fold enrichment`))/sd(`Fold enrichment`)) -> mresults_match
colnames(mresults_match)[5] <- "Z-score"
mresults_match$`Male baseline modules` <- factor(mresults_match$`Male baseline modules`,
                                                 levels=c("M1","M2","M3","M4",
                                                          "M5","M6","M7","M8",
                                                          "M9","M10"))

# Plot:
library(ggplot2)
par(mar=(c(1,1,1,1)))
pdf("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/male_mset_tissue.pdf",
    width = 6.5, height = 5)
p = ggplot(mresults_match[which(mresults_match$`Fold enrichment`>0),], aes(x=`Male baseline modules`,y=`Tissues`)) +
  geom_point(aes(size = `Fold enrichment`, colour = `Adjusted p-value`), shape = 16) +
  theme_bw() +
  ggtitle("Midgestation brain cell types enrichment") +
  theme(plot.title = element_text(hjust = 0.65), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_continuous(limits=c(min(mresults_match$`Adjusted p-value`),0.05), high = "#79ABF9", low = "#034DA2", na.value = "grey70") +
  scale_y_discrete(limits=rev)
p
dev.off()

################################## Disorders:
results_disorders = read.csv("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/MSET_results_modules - disorders.csv", header = TRUE)
results_disorders_female <- results_disorders[results_disorders$sex == "Female",]
results_disorders_male <- results_disorders[results_disorders$sex == "Male",]

fresults <- data.frame()
for (i in unique(results_disorders_female$disorder_database)){
  a <- results_disorders_female[results_disorders_female$disorder_database == i,]
  a$adj.p.valor <- p.adjust(a$p.valor, method = "BH")
  fresults <- rbind(fresults,a)
}

mresults <- data.frame()
for (i in unique(results_disorders_male$disorder_database)){
  a <- results_disorders_male[results_disorders_male$disorder_database == i,]
  a$adj.p.valor <- p.adjust(a$p.valor, method = "BH")
  mresults <- rbind(mresults,a)
}

### Female
# picking just the number of matches and renaming the columns
fresults$module <- paste0("M",fresults$module)
fresults_match = fresults[,c(1,5,11,8)]
colnames(fresults_match) = c("Disorders dataset","Fold enrichment", "Adjusted p-value", "Female baseline modules")
fresults_match %>% 
  mutate(zscore = (`Fold enrichment` - mean(`Fold enrichment`))/sd(`Fold enrichment`)) -> fresults_match
colnames(fresults_match)[5] <- "Z-score"

# Plot:
library(ggplot2)
pdf("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/female_mset_disorders.pdf",
    width = 6.5, height = 5)
p = ggplot(fresults_match[which(fresults_match$`Fold enrichment`>0),], aes(x=`Female baseline modules`,y=`Disorders dataset`)) +
  geom_point(aes(size = `Fold enrichment`, colour = `Adjusted p-value`), shape = 16) +
  theme_bw() +
  ggtitle("Disorders database enrichment") +
  theme(plot.title = element_text(hjust = 0.65), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_continuous(limits=c(min(fresults_match$`Adjusted p-value`),0.05), high = "#F9A499", low = "#C13624", na.value = "grey80") +
  scale_y_discrete(limits=rev)
p
dev.off()


### Male
# picking just the number of matches and renaming the columns
mresults$module <- paste0("M",mresults$module)
mresults_match = mresults[,c(1,5,11,8)]
colnames(mresults_match) = c("Disorders dataset","Fold enrichment", "Adjusted p-value", "Male baseline modules")
mresults_match %>% 
  mutate(zscore = (`Fold enrichment` - mean(`Fold enrichment`))/sd(`Fold enrichment`)) -> mresults_match
colnames(mresults_match)[5] <- "Z-escore"
mresults_match$`Male baseline modules` <- factor(mresults_match$`Male baseline modules`,
                                                 levels=c("M1","M2","M3","M4",
                                                          "M5","M6","M7","M8",
                                                          "M9","M10"))

# Plot:
library(ggplot2)
par(mar=(c(1,1,1,1)))
pdf("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/male_mset_disorders.pdf",
    width = 6.5, height = 5)
p = ggplot(mresults_match[which(mresults_match$`Fold enrichment`>0),], aes(x=`Male baseline modules`,y=`Disorders dataset`)) +
  geom_point(aes(size = `Fold enrichment`, colour = `Adjusted p-value`), shape = 16) +
  theme_bw() +
  ggtitle("Disorders database enrichment") +
  theme(plot.title = element_text(hjust = 0.65), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_continuous(limits=c(min(mresults_match$`Adjusted p-value`),0.05), high = "#79ABF9", low = "#034DA2", na.value = "grey70") +
  scale_y_discrete(limits=rev)
p
dev.off()

# Final disorders MSET results with adjusted p-value
final_disorders_results <- rbind(fresults,mresults)
write.table(final_disorders_results, file = paste0("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/MSET/final_disorders_MSET_results.txt"),
            sep = "\t", row.names = F, quote = F)
################################## hormones:
results_hormones = read.csv("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/MSET_results_modules - hormones.csv", header = TRUE)
results_hormones_female <- results_hormones[results_hormones$Baseline.Sex == "Female",]
results_hormones_male <- results_hormones[results_hormones$Baseline.Sex == "Male",]

fresults <- data.frame()
for (i in unique(results_hormones_female$Steroid.hormone.effector)){
  a <- results_hormones_female[results_hormones_female$Steroid.hormone.effector == i,]
  a$Adjusted.p.value <- p.adjust(a$p.value, method = "BH")
  fresults <- rbind(fresults,a)
}

mresults <- data.frame()
for (i in unique(results_hormones_male$Steroid.hormone.effector)){
  a <- results_hormones_male[results_hormones_male$Steroid.hormone.effector == i,]
  a$Adjusted.p.value <- p.adjust(a$p.value, method = "BH")
  mresults <- rbind(mresults,a)
}

### Female
# picking just the number of matches and renaming the columns
#fresults$module <- paste0("M",fresults$module)
fresults_match = fresults[,c(1,5,11,8)]
colnames(fresults_match) = c("hormones dataset","Fold enrichment", "Adjusted p-value", "Female baseline modules")
fresults_match %>% 
  mutate(zscore = (`Fold enrichment` - mean(`Fold enrichment`))/sd(`Fold enrichment`)) -> fresults_match
colnames(fresults_match)[5] <- "Z-score"

# Plot:
library(ggplot2)
pdf("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/female_mset_hormones.pdf",
    width = 6.5, height = 5)
p = ggplot(fresults_match[which(fresults_match$`Fold enrichment`>0),], aes(x=`Female baseline modules`,y=`hormones dataset`)) +
  geom_point(aes(size = `Fold enrichment`, colour = `Adjusted p-value`), shape = 16) +
  theme_bw() +
  ggtitle("hormones database enrichment") +
  theme(plot.title = element_text(hjust = 0.65), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_continuous(limits=c(min(fresults_match$`Adjusted p-value`),0.05), high = "#F9A499", low = "#C13624", na.value = "grey80") +
  scale_y_discrete(limits=rev)
p
dev.off()


### Male
# picking just the number of matches and renaming the columns
mresults$module <- paste0("M",mresults$module)
mresults_match = mresults[,c(1,5,11,8)]
colnames(mresults_match) = c("hormones dataset","Fold enrichment", "Adjusted p-value", "Male baseline modules")
mresults_match %>% 
  mutate(zscore = (`Fold enrichment` - mean(`Fold enrichment`))/sd(`Fold enrichment`)) -> mresults_match
colnames(mresults_match)[5] <- "Z-escore"
mresults_match$`Male baseline modules` <- factor(mresults_match$`Male baseline modules`,
                                                 levels=c("M1","M2","M3","M4",
                                                          "M5","M6","M7","M8",
                                                          "M9","M10"))

# Plot:
library(ggplot2)
par(mar=(c(1,1,1,1)))
pdf("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/male_mset_hormones.pdf",
    width = 6.5, height = 5)
p = ggplot(mresults_match[which(mresults_match$`Fold enrichment`>0),], aes(x=`Male baseline modules`,y=`hormones dataset`)) +
  geom_point(aes(size = `Fold enrichment`, colour = `Adjusted p-value`), shape = 16) +
  theme_bw() +
  ggtitle("hormones database enrichment") +
  theme(plot.title = element_text(hjust = 0.65), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_continuous(limits=c(min(mresults_match$`Adjusted p-value`),0.05), high = "#79ABF9", low = "#034DA2", na.value = "grey70") +
  scale_y_discrete(limits=rev)
p
dev.off()

# Final hormone MSET results with adjusted p-value
final_hormone_results <- rbind(fresults,mresults)
write.table(final_hormone_results, file = paste0("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/MSET/final_hormone_MSET_results.txt"),
            sep = ",", row.names = F, quote = F)
