###########################################
#title: "Analyzing ALPACA results by differential modularity matrix"
#author: "Victor Hugo Calegari de Toledo"
#date: '2021-09-10'
##########################################

#Load all the packages and set the folder:
setwd("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca")
options(stringsAsFactors = F)
options(bitmapType = 'cairo')
library(ggplot2)
library(ggrepel)
library(reshape2)

########################################## Male networks
###### Filter TFs and genes:
mcontrol <- read.delim("malecontrol_ALPACA_scores.txt", header = F, check.names = F)
mcontrol <- mcontrol[!mcontrol$V1 == "",]
mcontrol <- mcontrol[order(mcontrol$V2,decreasing = TRUE),]
colnames(mcontrol) <- c("Biological entity","Score")
#mcontrol$col <- "gray"
#mcontrol$col[mcontrol$Escore > 0.1] <- "midnightblue"
mcontrol$IDs <- gsub("_A","",mcontrol$`Biological entity`)
mcontrol$IDs <- gsub("_B","",mcontrol$IDs)
mcontrol$Type <- "Gene"
mcontrol$Type[grep("_A$",mcontrol$`Biological entity`)] <- "TF"

# Where are those elements?
membership <- read.delim("malecontrol_ALPACA_final_memb.txt", header = F, check.names = F)
mcontrol$Membership <- membership$V2[match(mcontrol$`Biological entity`, membership$V1)]

################### Top TFs by IQR
mcontrol_tfs <- mcontrol[mcontrol$Type == "TF",]
mcontrol_tfs$logScore <- log(mcontrol_tfs$Score)
mcontrol_tfs <- na.omit(mcontrol_tfs)
mcontrol_tfs$Membership <- as.factor(mcontrol_tfs$Membership)
mcontrol_tfs$`Z-score` <- NA
mcontrol_tfs$`Z-score` <- (mcontrol_tfs$logScore - mean(mcontrol_tfs$logScore)/sd(mcontrol_tfs$logScore))

p <- ggplot(mcontrol_tfs, aes(x=Membership, y=logScore)) + 
  geom_violin()
p

# # Select by total:
cut <- median(mcontrol_tfs$logScore) + 1.5*(IQR(mcontrol_tfs$logScore))
male_top_tfs <- mcontrol_tfs[mcontrol_tfs$logScore > cut,]

write.table(male_top_tfs, "~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/top_tfs_male.txt",
            sep = "\t", row.names = F, quote = F)

################### Top genes by IQR
mcontrol_genes <- mcontrol[mcontrol$Type == "Gene",]
mcontrol_genes$logScore <- log(mcontrol_genes$Score)
mcontrol_genes <- na.omit(mcontrol_genes)
mcontrol_genes$Membership <- as.factor(mcontrol_genes$Membership)

p <- ggplot(mcontrol_genes, aes(x=Membership, y=logScore)) + 
  geom_violin()
p

# Select by total:
cut <- median(mcontrol_genes$logScore) + 3*(IQR(mcontrol_genes$logScore))
male_top_genes <- mcontrol_genes[mcontrol_genes$logScore > cut,]

write.table(male_top_genes, "~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/top_genes_male.txt",
            sep = "\t", row.names = F, quote = F)

########################################## Female networks
setwd("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca")
fcontrol <- read.delim("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/femalecontrol_ALPACA_scores.txt", header = F, check.names = F)
fcontrol <- fcontrol[!fcontrol$V1 == "",]
fcontrol <- fcontrol[order(fcontrol$V2,decreasing = TRUE),]
colnames(fcontrol) <- c("Biological entity","Score")
#fcontrol$col <- "gray"
#fcontrol$col[fcontrol$Escore > 0.1] <- "darkred"
fcontrol$IDs <- gsub("_A","",fcontrol$`Biological entity`)
fcontrol$IDs <- gsub("_B","",fcontrol$IDs)

fcontrol$Type <- "Gene"
fcontrol$Type[grep("_A$",fcontrol$`Biological entity`)] <- "TF"

# Where are those elements?
fmembership <- read.delim("~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/femalecontrol_ALPACA_final_memb.txt", header = F, check.names = F)
fcontrol$Membership <- fmembership$V2[match(fcontrol$`Biological entity`, fmembership$V1)]

################### Top TFs by IQR
fcontrol_tfs <- fcontrol[fcontrol$Type == "TF",]
fcontrol_tfs$logScore <- log(fcontrol_tfs$Score)
fcontrol_tfs <- na.omit(fcontrol_tfs)
fcontrol_tfs$Membership <- as.factor(fcontrol_tfs$Membership)
fcontrol_tfs$`Z-score` <- NA
fcontrol_tfs$`Z-score` <- (fcontrol_tfs$logScore - mean(fcontrol_tfs$logScore)/sd(fcontrol_tfs$logScore))

# Select by total:
cut <- median(fcontrol_tfs$logScore) + 1.5*(IQR(fcontrol_tfs$logScore))
female_top_tfs <- fcontrol_tfs[fcontrol_tfs$logScore > cut,]

write.table(female_top_tfs, "~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/top_tfs_female.txt",
            sep = "\t", row.names = F, quote = F)

################### Top genes by IQR
fcontrol_genes <- fcontrol[fcontrol$Type == "Gene",]
fcontrol_genes$logScore <- log(fcontrol_genes$Score)
fcontrol_genes <- na.omit(fcontrol_genes)
fcontrol_genes$Membership <- as.factor(fcontrol_genes$Membership)

cut_genes <- median(fcontrol_genes$logScore) + 3*IQR(fcontrol_genes$logScore)
female_top_genes <- fcontrol_genes[fcontrol_genes$logScore > cut_genes,]
write.table(female_top_genes, "~/Documents/projects/ms-victor/sex_differences_fetal_brain/alpaca/top_genes_female.txt",
            sep = "\t", row.names = F, quote = F)


#### Shared top TFs:
shared_tfs <- male_top_tfs[male_top_tfs$IDs %in% female_top_tfs$IDs,]
table(shared_tfs$Membership)
shared_tfs$`Female Membership` <- female_top_tfs$IDs[match(shared_tfs$IDs, female_top_tfs$IDs),]