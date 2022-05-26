###########################################
#title: "Analyzing ALPACA results by differential modularity matrix"
#author: "Victor Hugo Calegari de Toledo"
#date: '2021-09-10'
##########################################

#Load all the packages and set the folder:
setwd("~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical_notf")
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

#Enrichment per se:
setwd("~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical_notf/enrichment_notf")

universe <- unique(mcontrol$IDs)
for (i in unique(mcontrol$Membership)) {
  
  x <- subset(mcontrol,Membership == i)
  x <- subset(x, Type == "Gene")
  type <- paste0(i,"_",nrow(x))
  x <- x$IDs
  
  ###########
  ##GeneOntology
  ###########
  
  ego_BP <- enrichGO(gene     = x, #Biological Process
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     keyType       = 'SYMBOL',
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.05,
                     minGSSize=5)
  #pool = TRUE)
  if (nrow(ego_BP)>=1) {
    write.table(ego_BP,paste(type,"ego_BP.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
  
  
  ego_MF <- enrichGO(gene          = x, #Molecular Function
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "MF",
                     keyType       = 'SYMBOL',
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.05)
  if (nrow(ego_MF)>=1) {
    write.table(ego_MF,paste(type,"ego_MF.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
  
  
  ego_CC <- enrichGO(gene          = x, #Cellular Component
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     keyType       = 'SYMBOL',
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.05)
  
  if (nrow(ego_CC)>=1) {
    write.table(ego_CC,paste(type,"ego_CC.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
  
  
  ##########
  #K E G G 
  ##########
  
  gene.df <- bitr(universe, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  y <- bitr(x, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
  
  kk <- enrichKEGG(gene         = y$ENTREZID,
                   organism     = 'hsa',
                   universe = gene.df$ENTREZID,
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05)
  #  qvalueCutoff  = 0.05)
  if (nrow(kk)>=1) {
    write.table(kk,paste(type,"KEGG.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
  
  
  ###########
  #REACTOME
  ###########
  reac <- enrichPathway(gene=y$ENTREZID,
                        organism = "human",
                        universe = gene.df$ENTREZID,
                        pAdjustMethod = "fdr",
                        pvalueCutoff=0.05)
  if (nrow(reac)>=1) {
    write.table(reac,paste(type,"REACTOME.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
}

################### Top TFs by IQR
mcontrol_tfs <- mcontrol[mcontrol$Type == "TF",]
mcontrol_tfs$logScore <- log(mcontrol_tfs$Score)
mcontrol_tfs <- na.omit(mcontrol_tfs)
mcontrol_tfs$Membership <- as.factor(mcontrol_tfs$Membership)
mcontrol_tfs$`Z-score` <- NA
mcontrol_tfs$`Z-score` <- (mcontrol_tfs$logScore - mean(mcontrol_tfs$logScore)/sd(mcontrol_tfs$logScore))

# tfs_male <- data.frame()
# 
# for (i in 1:max(mcontrol$Membership)){
#   x <- mcontrol_tfs[mcontrol_tfs$Membership == i,]
#   x$`Z-score` <- (x$logScore - mean(x$logScore)/sd(x$logScore))
#   tfs_male <- rbind(tfs_male,x)
# }
  
p <- ggplot(mcontrol_tfs, aes(x=Membership, y=logScore)) + 
  geom_violin()
p

# Loop to select top TFs:
# male_top_tfs <- data.frame()
# 
# for (i in 1:max(mcontrol$Membership)){
#   a <- tfs_male[tfs_male$Membership == i,]
#   cut <- median(a$logScore) + 1*(IQR(a$logScore))
#   b <- a[a$logScore > cut,]
#   male_top_tfs <- rbind(male_top_tfs,b)
#   print(i)
# }

# # Select by total:
cut <- median(mcontrol_tfs$logScore) + 1.5*(IQR(mcontrol_tfs$logScore))
male_top_tfs <- mcontrol_tfs[mcontrol_tfs$logScore > cut,]

write.table(male_top_tfs, "~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical_notf/top_tfs_male.txt",
            sep = "\t", row.names = F, quote = F)

################### Top genes by IQR
mcontrol_genes <- mcontrol[mcontrol$Type == "Gene",]
mcontrol_genes$logScore <- log(mcontrol_genes$Score)
mcontrol_genes <- na.omit(mcontrol_genes)
mcontrol_genes$Membership <- as.factor(mcontrol_genes$Membership)

# genes_male <- data.frame()
# 
# for (i in 1:max(mcontrol$Membership)){
#   x <- mcontrol_genes[mcontrol_genes$Membership == i,]
#   x$`Z-score` <- (x$logScore - mean(x$logScore)/sd(x$logScore))
#   genes_male <- rbind(genes_male,x)  
# }

p <- ggplot(mcontrol_genes, aes(x=Membership, y=logScore)) + 
  geom_violin()
p

# Loop to select top genes:
# male_top_genes <- data.frame()
# 
# for (i in 1:max(mcontrol$Membership)){
#   a <- genes_male[genes_male$Membership == i,]
#   cut <- median(a$logScore) + 1.5*(IQR(a$logScore))
#   b <- a[a$logScore > cut,]
#   male_top_genes <- rbind(male_top_genes,b)
#   print(i)
# }

# Select by total:
cut <- median(mcontrol_genes$logScore) + 3*(IQR(mcontrol_genes$logScore))
male_top_genes <- mcontrol_genes[mcontrol_genes$logScore > cut,]

write.table(male_top_genes, "~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical_notf/top_genes_male.txt",
            sep = "\t", row.names = F, quote = F)

# #Plotting distribution
# options(ggrepel.max.overlaps = Inf)
# pdf("top_tf_male.pdf", width = 8, height = 6)
# male_plot <- ggplot(a, aes(x = seq_along(Score),
#                                   y = Score, col = col), show.legend = F) +
#   geom_point() + theme_bw() +
#   scale_color_manual(values=c("gray" = "gray",
#                               "midnightblue" = "midnightblue")) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.title = element_blank(), legend.position = "none") +
#   labs (x = "Elementos", 
#         y = "Escores") +
#   # geom_text_repel(data = subset(a, Score > 0.01),
#   #                 aes(label = IDs),
#   #                 fontface = 'bold', color = 'grey30',
#   #                 segment.color = 'grey80',
#   #                 size = 3,
#   #                 box.padding = unit(0.5, "lines"),
#   #                 point.padding = unit(0.6, "lines")) +
#   ggtitle("Escore de contribuição para diferença de modularidade (controle masculino)") +
#   geom_vline(xintercept=1000) +
#   geom_vline(xintercept=2500)
# male_plot
# dev.off()

mnetwork <- read.delim("malecontrol_DWBM.txt", check.names = F)
mnetwork <- mnetwork[rownames(mnetwork) %in% m_top$Elemento,]
mod1 <- read.delim("module_1_top_alpaca_genes_male.txt", header = F, check.names = F)
mod2 <- read.delim("module_2_top_alpaca_genes_male.txt", header = F, check.names = F)
top_genes <- c(mod1$V1,mod2$V1)
top_genes <- paste0(top_genes,"_B")
mnetwork <- mnetwork[,colnames(mnetwork) %in% top_genes]
mnetwork$TFs <- rownames(mnetwork)
male_cytoscape <- melt(mnetwork, id = "TFs")

#Plot
male_cytoscape <- male_cytoscape[order(male_cytoscape$value,decreasing = TRUE),]
ggplot(subnet, aes(x = seq_along(value),
                     y = value), show.legend = F) +
  geom_point() + theme_bw() +
  #scale_color_manual(values=c("gray" = "gray",
  #                            "midnightblue" = "midnightblue")) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(), legend.position = "none") +
  labs (x = "Elementos", 
        y = "Escores")

subnet <- male_cytoscape[male_cytoscape$value > 0,]

ggplot(subnet, aes(x = variable,
                           y = value), show.legend = F) +
  geom_bar(stat='identity') + theme_bw() +
  #scale_color_manual(values=c("gray" = "gray",
  #                            "midnightblue" = "midnightblue")) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(), legend.position = "none") +
  labs (x = "Elementos", 
        y = "Escores")


male_cytoscape <- male_cytoscape[male_cytoscape$value > 0,]
male_cytoscape$TFs <- gsub("_A","",male_cytoscape$TFs)
male_cytoscape$variable <- gsub("_B","",male_cytoscape$variable)
male_cytoscape$module <- 1
male_cytoscape$module[male_cytoscape$variable %in% mod2$V1] <- 2
colnames(male_cytoscape) <- c("SOURCE","TARGET","WEIGTH","MODULE")
write.table(male_cytoscape, "~/Documents/projects/ms-victor/sex_differences_fetal_brain/netZoo/alpaca/_m/cytoscape_input/male_network_modules_cytoscape.txt",
            sep = "\t", row.names = F, quote = F)


########################################## Female networks
setwd("~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical_notf")
fcontrol <- read.delim("~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical/femalecontrol_ALPACA_scores.txt", header = F, check.names = F)
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
fmembership <- read.delim("~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical/femalecontrol_ALPACA_final_memb.txt", header = F, check.names = F)
fcontrol$Membership <- fmembership$V2[match(fcontrol$`Biological entity`, fmembership$V1)]

#Enrichment per se:
setwd("~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical_notf/enrichment_notf")

universe <- unique(fcontrol$IDs)
for (i in unique(fcontrol$Membership)) {
  
  x <- subset(fcontrol,Membership == i)
  x <- subset(x, Type == "Gene")
  type <- paste0(i,"_",nrow(x))
  x <- x$IDs
  
  ###########
  ##GeneOntology
  ###########
  
  ego_BP <- enrichGO(gene     = x, #Biological Process
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     keyType       = 'SYMBOL',
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.05,
                     minGSSize=5)
  #pool = TRUE)
  if (nrow(ego_BP)>=1) {
    write.table(ego_BP,paste(type,"ego_BP_fem.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
  
  
  ego_MF <- enrichGO(gene          = x, #Molecular Function
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "MF",
                     keyType       = 'SYMBOL',
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.05)
  if (nrow(ego_MF)>=1) {
    write.table(ego_MF,paste(type,"ego_MF_fem.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
  
  
  ego_CC <- enrichGO(gene          = x, #Cellular Component
                     universe      = universe,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     keyType       = 'SYMBOL',
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.05)
  
  if (nrow(ego_CC)>=1) {
    write.table(ego_CC,paste(type,"ego_CC_fem.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
  
  
  ##########
  #K E G G 
  ##########
  
  gene.df <- bitr(universe, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  y <- bitr(x, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
  
  kk <- enrichKEGG(gene         = y$ENTREZID,
                   organism     = 'hsa',
                   universe = gene.df$ENTREZID,
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05)
  #  qvalueCutoff  = 0.05)
  if (nrow(kk)>=1) {
    write.table(kk,paste(type,"KEGG_fem.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
  
  
  ###########
  #REACTOME
  ###########
  reac <- enrichPathway(gene=y$ENTREZID,
                        organism = "human",
                        universe = gene.df$ENTREZID,
                        pAdjustMethod = "fdr",
                        pvalueCutoff=0.05)
  if (nrow(reac)>=1) {
    write.table(reac,paste(type,"REACTOME_fem.txt",sep = "-"),sep="\t",quote=F) } else {
      paste0("NULL_",type)
    }
}



################### Top TFs by IQR
fcontrol_tfs <- fcontrol[fcontrol$Type == "TF",]
fcontrol_tfs$logScore <- log(fcontrol_tfs$Score)
fcontrol_tfs <- na.omit(fcontrol_tfs)
fcontrol_tfs$Membership <- as.factor(fcontrol_tfs$Membership)
fcontrol_tfs$`Z-score` <- NA
fcontrol_tfs$`Z-score` <- (fcontrol_tfs$logScore - mean(fcontrol_tfs$logScore)/sd(fcontrol_tfs$logScore))


# tfs_female <- data.frame()
# 
# for (i in 1:max(fcontrol$Membership)){
#   x <- fcontrol_tfs[fcontrol_tfs$Membership == i,]
#   x$`Z-score` <- (x$logScore - mean(x$logScore)/sd(x$logScore))
#   tfs_female <- rbind(tfs_female,x)  
# }
# 
# p <- ggplot(fcontrol_tfs, aes(x=Membership, y=logScore)) + 
#   geom_violin()
# p
# 
# # Loop to select top TFs:
# female_top_tfs <- data.frame()
# 
# for (i in 1:max(fcontrol$Membership)){
#   a <- tfs_female[tfs_female$Membership == i,]
#   cut <- median(a$logScore) + 1*(IQR(a$logScore))
#   b <- a[a$logScore > cut,]
#   female_top_tfs <- rbind(female_top_tfs,b)
#   print(i)
# }

# Select by total:
cut <- median(fcontrol_tfs$logScore) + 1.5*(IQR(fcontrol_tfs$logScore))
female_top_tfs <- fcontrol_tfs[fcontrol_tfs$logScore > cut,]

write.table(female_top_tfs, "~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical_notf/top_tfs_female.txt",
            sep = "\t", row.names = F, quote = F)

################### Top genes by IQR
fcontrol_genes <- fcontrol[fcontrol$Type == "Gene",]
fcontrol_genes$logScore <- log(fcontrol_genes$Score)
fcontrol_genes <- na.omit(fcontrol_genes)
fcontrol_genes$Membership <- as.factor(fcontrol_genes$Membership)
# fcontrol_genes$`Z-score` <- NA
# 
# genes_female <- data.frame()
# 
# for (i in 1:max(fcontrol$Membership)){
#   x <- fcontrol_genes[fcontrol_genes$Membership == i,]
#   x$`Z-score` <- (x$logScore - mean(x$logScore)/sd(x$logScore))
#   genes_female <- rbind(genes_female,x)  
# }
# 
# 
# p <- ggplot(genes_female, aes(x=Membership, y=logScore)) + 
#   geom_violin()
# p
# 
# # Loop to select top genes:
# female_top_genes <- data.frame()
# 
# for (i in 1:max(fcontrol$Membership)){
#   a <- genes_female[genes_female$Membership == i,]
#   cut <- median(a$logScore) + 1.5*(IQR(a$logScore))
#   b <- a[a$logScore > cut,]
#   female_top_genes <- rbind(female_top_genes,b)
#   print(i)
# }

cut_genes <- median(fcontrol_genes$logScore) + 3*IQR(fcontrol_genes$logScore)
female_top_genes <- fcontrol_genes[fcontrol_genes$logScore > cut_genes,]
write.table(female_top_genes, "~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca_canonical_notf/top_genes_female.txt",
            sep = "\t", row.names = F, quote = F)

# #Selecting top TFs:
# options(ggrepel.max.overlaps = Inf)
# pdf("top_tf_female.pdf", width = 8, height = 6)
# female_plot <- ggplot(fcontrol, aes(x = seq_along(Escore),
#                                   y = Escore, col = col), show.legend = F) +
#   geom_point() + theme_bw() +
#   scale_color_manual(values=c("gray" = "gray",
#                               "darkred" = "darkred")) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.title = element_blank(), legend.position = "none") +
#   labs (x = "Elementos", 
#         y = "Escores") +
#   geom_text_repel(data = subset(fcontrol, Escore > 0.1),
#                   aes(label = IDs),
#                   fontface = 'bold', color = 'grey30',
#                   segment.color = 'grey80',
#                   size = 3,
#                   box.padding = unit(0.5, "lines"),
#                   point.padding = unit(0.6, "lines")) +
#   ggtitle("Escore de contribuição para diferença de modularidade (controle feminino)")
# female_plot
# dev.off()


#### Shared top TFs:
shared_tfs <- male_top_tfs[male_top_tfs$IDs %in% female_top_tfs$IDs,]
table(shared_tfs$Membership)
shared_tfs$`Female Membership` <- female_top_tfs$IDs[match(shared_tfs$IDs, female_top_tfs$IDs),]












mnetwork <- read.delim("malecontrol_DWBM.txt", check.names = F)
mnetwork <- mnetwork[rownames(mnetwork) %in% m_top$Elemento,]
mod1 <- read.delim("module_1_top_alpaca_genes_male.txt", header = F, check.names = F)
mod2 <- read.delim("module_2_top_alpaca_genes_male.txt", header = F, check.names = F)
top_genes <- c(mod1$V1,mod2$V1)
top_genes <- paste0(top_genes,"_B")
mnetwork <- mnetwork[,colnames(mnetwork) %in% top_genes]
mnetwork$TFs <- rownames(mnetwork)
male_cytoscape <- melt(mnetwork, id = "TFs")
male_cytoscape <- male_cytoscape[male_cytoscape$value > 0,]
male_cytoscape$TFs <- gsub("_A","",male_cytoscape$TFs)
male_cytoscape$variable <- gsub("_B","",male_cytoscape$variable)
male_cytoscape$module <- 1
male_cytoscape$module[male_cytoscape$variable %in% mod2$V1] <- 2
colnames(male_cytoscape) <- c("SOURCE","TARGET","WEIGTH","MODULE")
write.table(male_cytoscape, "~/Documents/projects/ms-victor/sex_differences_fetal_brain/netZoo/alpaca/_m/cytoscape_input/male_network_modules_cytoscape.txt",
            sep = "\t", row.names = F, quote = F)