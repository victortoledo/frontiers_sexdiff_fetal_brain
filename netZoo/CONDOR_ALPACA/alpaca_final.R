###########################################
#title: "Running ALPACA on PANDA networks"
#author: "Victor Hugo Calegari de Toledo"
#date: '2021-09-10'
##########################################

#Load all the packages and set the folder:
setwd("~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca")
options(stringsAsFactors = F)
options(bitmapType = 'cairo')
library(data.table)
library(netZooR)
library(dplyr)
library(org.Hs.eg.db)
library(gtools)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(stringr)
library(biomaRt)
library(fgsea)
library(ReactomePA)

############################################### Get each regulatory network:
f <- fread("~/Documents/output_panda_female.txt")
m <- fread("~/Documents/output_panda_male.txt")
colnames(m) = c("tf","gene","canonical","score")
colnames(f) = c("tf","gene","canonical","score")

# Only canonical interactions
f <- f[f$canonical == 1,]
m <- m[m$canonical == 1,]

############################################### Using male samples as control:
malecontrol = cbind(m[,c(1,2,4)],f[,4])
#rm(f,m)
#gc()
colnames(malecontrol) = c("tf","gene","control_score","perturbed_score")
malecontrol = as.data.frame(malecontrol)
malealp = alpaca(net.table=malecontrol, file.stem = "~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca/malecontrol",
                 verbose = T)
save(malealp,file = "male_alpaca.rda")

###### Running enrichment on all elements:
mod.memb <- malealp[[1]]
mod.scores <- malealp[[2]]
mod.ord <- names(mod.scores)[order(mod.scores, decreasing = TRUE)]
mod.top <- NULL
mod.top.names <- NULL
count <- 0
lengths <- table(mod.memb)

  for (i in 1:max(mod.memb)) {
    set.lengths <- as.numeric(lengths[i])
    mod.top.names <- c(mod.top.names, paste(i, set.lengths, 
                                            sep = "_"))
    this.comm <- names(mod.memb)[mod.memb == i]
    this.comm.ord <- mod.ord[mod.ord %in% this.comm]
    for (j in 1:length(set.lengths)) {
      count <- count + 1
      if (length(this.comm.ord) < set.lengths[j]) 
        mod.top[[count]] <- sapply(this.comm.ord, alpaca.node.to.gene)
      else mod.top[[count]] <- sapply(this.comm.ord[1:set.lengths[j]], 
                                      alpaca.node.to.gene)
    }
  }

malealp_topgene <- list(mod.top, mod.top.names)

#Removing the NAs and nulls manually
malealp_topgene[[1]] <- malealp_topgene[[1]][-c(rep(10:15))]
malealp_topgene[[2]] <- malealp_topgene[[2]][-c(rep(10:15))]
#malealp_topgene[sapply(malealp_topgene, is.null)] <- NA

for(i in 1:length(malealp_topgene[[1]])){
  dat <- unlist(malealp_topgene[[1]][i])
  dat <- unname(dat)
  write.table(dat, paste0("module_",as.numeric(i),"_all_alpaca_genes_male.txt"), 
            col.names = F, quote = F, row.names = F)
}

universe <- unique(malecontrol$gene)
write.table(universe, paste0("universe_alpaca_genes_male.txt"), 
            col.names = F, quote = F, row.names = F)

#Enrichment per se:
setwd("~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca")

for (i in files) {
  
  type <- i
  type <- gsub(".txt","",type)
  
  #x = as.data.frame(malealp_topgene[[1]][[i]])
  #x = x[grep("_B$", row.names(x)),]
  
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

############################################### Using female samples as control:
setwd("~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca")
# f <- fread("~/Documents/output_panda_female.txt")
# m <- fread("~/Documents/output_panda_male.txt")
# colnames(m) = c("tf","gene","canonical","score")
# colnames(f) = c("tf","gene","canonical","score")

femalecontrol = cbind(f[,c(1,2,4)],m[,4])
rm(m)
femalecontrol = as.data.frame(femalecontrol)
femalp = alpaca(net.table=femalecontrol, file.stem = "~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca/femalecontrol",
                verbose = T)
save(femalp,file = "~/Documents/Documents - Victor’s MacBook Pro/projects/ms-victor/sex_differences_fetal_brain/alpaca/female_alpaca.rda")

###### Running enrichment on scores:
mod.memb <- femalp[[1]]
mod.scores <- femalp[[2]]
mod.ord <- names(mod.scores)[order(mod.scores, decreasing = TRUE)]
mod.top <- NULL
mod.top.names <- NULL
count <- 0
lengths <- table(mod.memb)

for (i in 1:max(mod.memb)) {
  set.lengths <- as.numeric(lengths[i])
  mod.top.names <- c(mod.top.names, paste(i, set.lengths, 
                                          sep = "_"))
  this.comm <- names(mod.memb)[mod.memb == i]
  this.comm.ord <- mod.ord[mod.ord %in% this.comm]
  for (j in 1:length(set.lengths)) {
    count <- count + 1
    if (length(this.comm.ord) < set.lengths[j]) 
      mod.top[[count]] <- sapply(this.comm.ord, alpaca.node.to.gene)
    else mod.top[[count]] <- sapply(this.comm.ord[1:set.lengths[j]], 
                                    alpaca.node.to.gene)
  }
}

femalealp_topgene <- list(mod.top, mod.top.names)
universe <- unique(f$gene)

#Removing the NAs and nulls manually
femalealp_topgene[[1]] <- femalealp_topgene[[1]][-c(rep(10:17))]
femalealp_topgene[[2]] <- femalealp_topgene[[2]][-c(rep(10:17))]
#malealp_topgene[sapply(malealp_topgene, is.null)] <- NA

for(i in 1:length(femalealp_topgene[[1]])){
  dat <- unlist(femalealp_topgene[[1]][i])
  dat <- unname(dat)
  write.table(dat, paste0("module_",as.numeric(i),"_all_alpaca_genes_female.txt"),
              col.names = F, quote = F, row.names = F)
}

#Enrichment per se:
setwd("~/Documents/projects/ms-victor/sex_differences_fetal_brain/netZoo/alpaca/_m/enrichment")
for (i in 1:length(femalealp_topgene[[1]])) {
  
  type <- femalealp_topgene[[2]][i]
  
  x = as.data.frame(femalealp_topgene[[1]][[i]])
  colnames(x) = "GENE_SYMBOL"
  
  ###########
  ##GeneOntology
  ###########
  
  ego_BP <- enrichGO(gene     = x$GENE_SYMBOL, #Biological Process
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
  
  
  ego_MF <- enrichGO(gene          = x$GENE_SYMBOL, #Molecular Function
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
  
  
  ego_CC <- enrichGO(gene          = x$GENE_SYMBOL, #Cellular Component
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
  y <- bitr(x$GENE_SYMBOL, fromType = "SYMBOL",
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

####################################################### Creating readable files:
setwd("~/Documents/projects/ms-victor/sex_differences_fetal_brain/netZoo/alpaca/_m/")
files = list.files(pattern = "final_memb|scores")

## Elements with _A are transcription factors and with _B are genes
alpaca = lapply(files, function(x){
  tab = read.delim(x, sep='\t', header=FALSE)
  if(is.integer(tab[,2])){
    colnames(tab) = c('element', 'membership')
  }else{colnames(tab) = c('element', 'score')}
  return(tab)
})
names(alpaca) = files

###################################### Merge gene membership with scores - for each sex
# MALE
alpaca[[grep("^male.*memb", files)]] = alpaca[[grep("^male.*memb", files)]]  %>%
  left_join(alpaca[[grep("^male.*scores", files)]], by="element")

# FEMALE
alpaca[[grep("^female.*memb", files)]] = alpaca[[grep("^female.*memb", files)]]  %>%
  left_join(alpaca[[grep("^female.*scores", files)]], by="element")

#Get merged tables 
alpaca = alpaca[grep('memb',names(alpaca))]

names(alpaca) <- c('female', 'male')
# Identify sex of dataset
alpaca[['male']]$sex <- 'male'
alpaca[['female']]$sex <- 'female'

# Remove modules and genes w/ no assigned score
#alpaca <- lapply(alpaca, function(x){
#keep <- x %>% group_by(membership) %>% summarize(sum=sum(score)) %>% filter(!is.na(sum))
#return (x[x$membership %in% keep$membership,])
#})

alpaca$female <- na.omit(alpaca$female) #18833
alpaca$male <- na.omit(alpaca$male) #18858

# Function to split modules
modul <- function(membership, module_tab){
  
  module<-lapply(1:length(membership), function(x) {
    alp <- module_tab[[x]]
    m <- memb[[x]]
    mods<-lapply(m, function(z){
      mod <- data.frame(alp[alp$membership %in% z,])
      mod[grep("_A$",mod$element),'gene_class'] <- 'TF'
      mod[grep("_A$",mod$element, invert=TRUE),'gene_class'] <- 'Target'
      mod <- mod[,c(1,3,2,5,4)]
      colnames(mod)[3] <- 'module'
      #mod <- mod %>% left_join(annot, by=c("gene"="ensembl_gene_id"))
      #mod[mod$gene_class == "TF", "hgnc_symbol"] <- mod[mod$gene_class == "TF", "gene"]
      return(mod)
    })
    return(mods)
  })
  return(module)
}

memb <- lapply(alpaca, function(x){unique(x$membership)})
mod <- modul(memb, alpaca)
dfmod <- lapply(mod, rbindlist)

#Saving the results
setwd("~/Documents/projects/ms-victor/sex_differences_fetal_brain/netZoo/alpaca")
dfmod[[1]]$element = gsub("_A","",dfmod[[1]]$element)
dfmod[[1]]$element = gsub("_B","",dfmod[[1]]$element)
dfmod[[2]]$element = gsub("_A","",dfmod[[2]]$element)
dfmod[[2]]$element = gsub("_B","",dfmod[[2]]$element)
write.table(dfmod[[1]], 'female_modules_genes.tsv',sep='\t', row.names=FALSE)
write.table(dfmod[[2]], 'male_modules_genes.tsv',sep='\t', row.names=FALSE)

################################################################# Network details:

# Size / Number of modules
library(hrbrthemes)
library(babynames)
library(viridis)
fem_scores <- read.delim('female_modules_genes.tsv', sep='\t')
men_scores <- read.delim('male_modules_genes.tsv', sep='\t')
fem <- data.frame(table(fem_scores$module,fem_scores$gene_class),"A. Female-baseline differential network")
colnames(fem) <- c("modules","element","value","sex")
men <- data.frame(table(men_scores$module,men_scores$gene_class),"B. Male-baseline differential network")
colnames(men) <- c("modules","element","value","sex")
input_plot <- rbind(fem,men)
levels(input_plot$element)[match("Target",levels(input_plot$element))] <- "Gene"

pdf("ALPACA_module_size_sex.pdf", width = 8, height = 6)
input_plot  %>%
  mutate(modules=as.factor(modules)) %>%
  ggplot( aes(x=modules, y=value, fill=element)) +
  geom_bar(position="stack", stat="identity", width = 0.5) +
  facet_wrap(~sex) +
  theme_bw() +
  scale_fill_manual(values=c("#4B0082","#008080")) +
  ylab("Tamanho dos modules")
dev.off()

# Size of each module:
table(fem_scores$module,fem_scores$gene_class)
# Target   TF
# 1   4196  154
# 2   5404  308
# 3   1289   54
# 4   1162   31
# 5    127    3
# 6   2649  176
# 7    115    4
# 8     85 3076
table(men_scores$module,men_scores$gene_class)
# Target   TF
# 1    4768  181
# 2    5361  279
# 3    2690  182
# 4    1604   74
# 5      54    0
# 6     177    6
# 7      77    2
# 8      42    0
# 9      96    0
# 10     83    3
# 11    103 3076
  