#------------------------------------------------------------------------
#title: "Cemitool"
#author: Victor Toledo victor.toledo.btos@gmail.com"
#date: "23/09/2021"
#------------------------------------------------------------------------

################### Load all the packages and set the folder:
setwd("~/projects/sex_differences_fetal_brain/cemitool")
options(stringsAsFactors = F)
options(bitmapType = 'cairo')
library(CEMiTool)
library(fgsea)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(WGCNA)
library(data.table)
##########################################################
# Reading input:
ppi <- read.delim("ppi_ready.txt",sep="\t",check.names = F)
expr <- read.delim("expression_ready.txt",sep="\t",check.names = F)
load("classes.Rda")

pathways <- gmtPathways("c2.cp.reactome.v7.4.symbols.gmt.txt")
pathways <- melt(pathways)
gmt_in <- pathways[,c(2,1)]
colnames(gmt_in) <- c("term","gene")

#Running CEMiTool:

sample_annot <- data.frame(SampleName = colnames(expr), Class = classes)

# run cemitool
cem <- cemitool(expr, annot = sample_annot, interactions=ppi, filter = F,
                apply_vst = F, plot_diagnostics = T,
                plot=T, verbose=TRUE, gsea_max_size = 6000)

hubs <- get_hubs(cem,10)
summary <- mod_summary(cem)

# write analysis results into files
write_files(cem, directory = '~/projects/sex_differences_fetal_brain/cemitool/Tables',force=TRUE)

# save all plots
cem <- plot_beta_r2(cem)
cem <- plot_hist(cem)
cem <- plot_interactions(cem)
cem <- plot_mean_k(cem)
cem <- plot_mean_var(cem)
cem <- plot_profile(cem)
cem <- plot_qq(cem)
cem <- plot_sample_tree(cem)
cem <- mod_ora(cem, gmt_in)
cem <- plot_ora(cem)
cem <- mod_gsea(cem, gsea_max_size = 6000)
cem <- plot_gsea(cem)

# create report as html document
generate_report(cem, directory = '~/projects/sex_differences_fetal_brain/cemitool/Report',force=TRUE)

save_plots(cem, 'all',
           directory = '~/projects/sex_differences_fetal_brain/cemitool/Plots',force=TRUE)

#save hubs
a <- data.frame()
for (i in names(hubs)){
  a <- rbind(a,names(hubs[[i]]))
  a <- rbind(a,hubs[[i]])
}
colnames(a) <- c("hub1","hub2","hub3","hub4","hub5","hub6","hub7","hub8","hub9","hub10")


write.csv(a, file = "~/projects/sex_differences_fetal_brain/cemitool/hubs.txt")

sessionInfo()

# R version 4.1.0 (2021-05-18)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.0.dev.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] gplots_3.1.1                Matrix_1.3-4                WGCNA_1.70-3               
# [4] fastcluster_1.2.3           dynamicTreeCut_1.63-1       reshape2_1.4.4             
# [7] forcats_0.5.1               purrr_0.3.4                 readr_2.0.0                
# [10] tidyr_1.1.3                 tibble_3.1.4                tidyverse_1.3.1            
# [13] fgsea_1.18.0                CEMiTool_1.16.0             ggupset_0.3.0              
# [16] enrichplot_1.12.2           UpSetR_1.4.0                cowplot_1.1.1              
# [19] GO.db_3.13.0                biomaRt_2.48.2              stringr_1.4.0              
# [22] DOSE_3.18.1                 clusterProfiler.dplyr_0.0.2 clusterProfiler_4.0.2      
# [25] ggplot2_3.3.5               gtools_3.9.2                org.Hs.eg.db_3.13.0        
# [28] AnnotationDbi_1.54.1        IRanges_2.26.0              S4Vectors_0.30.0           
# [31] dplyr_1.0.7                 netZooR_0.99.1              pandaR_1.19.4              
# [34] Biobase_2.52.0              BiocGenerics_0.38.0         reticulate_1.20            
# [37] igraph_1.2.6                data.table_1.14.0          
# 
# loaded via a namespace (and not attached):
#   [1] rappdirs_0.3.3         SparseM_1.81           pbdZMQ_0.3-5           AnnotationForge_1.34.0
# [5] ggthemes_4.2.4         R.methodsS3_1.8.1      coda_0.19-4            ggpmisc_0.4.0         
# [9] bit64_4.0.5            knitr_1.33             R.utils_2.10.1         rpart_4.1-15          
# [13] KEGGREST_1.32.0        RCurl_1.98-1.4         doParallel_1.0.16      generics_0.1.0        
# [17] preprocessCore_1.54.0  RSQLite_2.2.7          shadowtext_0.0.8       chron_2.3-56          
# [21] tzdb_0.1.2             bit_4.0.4              base64url_1.4          lubridate_1.7.10      
# [25] xml2_1.3.2             ggpp_0.4.1             assertthat_0.2.1       viridis_0.6.1         
# [29] STRINGdb_2.4.1         xfun_0.24              hms_1.1.0              evaluate_0.14         
# [33] fansi_0.5.0            progress_1.2.2         readxl_1.3.1           caTools_1.18.2        
# [37] dbplyr_2.1.1           Rgraphviz_2.36.0       DBI_1.1.1              htmlwidgets_1.5.3     
# [41] reshape_0.8.8          hash_2.2.6.1           ellipsis_0.3.2         crosstalk_1.1.1       
# [45] backports_1.2.1        signal_0.7-7           permute_0.9-5          annotate_1.70.0       
# [49] vctrs_0.3.8            quantreg_5.86          penalized_0.9-51       cachem_1.0.5          
# [53] withr_2.4.2            ggforce_0.3.3          checkmate_2.0.0        sna_2.6               
# [57] vegan_2.5-7            treeio_1.16.1          prettyunits_1.1.1      cluster_2.1.2         
# [61] ape_5.5                IRdisplay_1.0          lazyeval_0.2.2         crayon_1.4.1          
# [65] uchardet_1.1.0         genefilter_1.74.0      labeling_0.4.2         pkgconfig_2.0.3       
# [69] tweenr_1.0.2           GenomeInfoDb_1.28.1    nlme_3.1-152           nnet_7.3-16           
# [73] rlang_0.4.11           RJSONIO_1.3-1.5        lifecycle_1.0.0        MatrixModels_0.5-0    
# [77] downloader_0.4         filelock_1.0.2         BiocFileCache_2.0.0    modelr_0.1.8          
# [81] GOstats_2.58.0         cellranger_1.1.0       polyclip_1.10-0        matrixStats_0.60.0    
# [85] graph_1.70.0           aplot_0.0.6            IRkernel_1.2           reprex_2.0.1          
# [89] base64enc_0.1-3        png_0.1-7              viridisLite_0.4.0      bitops_1.0-7          
# [93] R.oo_1.24.0            KernSmooth_2.23-20     dplR_1.7.2             Biostrings_2.60.2     
# [97] blob_1.2.2             qvalue_2.24.0          jpeg_0.1-9             scales_1.1.1          
# [101] memoise_2.0.0          GSEABase_1.54.0        magrittr_2.0.1         plyr_1.8.6            
# [105] hexbin_1.28.2          zlibbioc_1.38.0        compiler_4.1.0         scatterpie_0.1.6      
# [109] RColorBrewer_1.1-2     plotrix_3.8-1          intergraph_2.0-2       cli_3.0.1             
# [113] XVector_0.32.0         Category_2.58.0        patchwork_1.1.1        htmlTable_2.2.1       
# [117] Formula_1.2-4          MASS_7.3-54            mgcv_1.8-36            tidyselect_1.1.1      
# [121] stringi_1.7.3          highr_0.9              yaml_2.2.1             GOSemSim_2.18.0       
# [125] latticeExtra_0.6-29    ggrepel_0.9.1          grid_4.1.0             polynom_1.4-0         
# [129] fastmatch_1.1-3        tools_4.1.0            rstudioapi_0.13        uuid_0.1-4            
# [133] foreach_1.5.1          foreign_0.8-81         gridExtra_2.3          farver_2.1.0          
# [137] ggraph_2.0.5           digest_0.6.27          rvcheck_0.1.8          BiocManager_1.30.16   
# [141] pracma_2.3.3           proto_1.0.0            Rcpp_1.0.7             broom_0.7.9           
# [145] httr_1.4.2             ggdendro_0.1.22        colorspace_2.0-2       rvest_1.0.1           
# [149] fs_1.5.0               XML_3.99-0.7           splines_4.1.0          RBGL_1.68.0           
# [153] tidytree_0.3.4         conquer_1.0.2          graphlayouts_0.7.1     xtable_1.8-4          
# [157] jsonlite_1.7.2         ggtree_3.0.2           tidygraph_1.2.0        R6_2.5.1              
# [161] RUnit_0.4.32           Hmisc_4.5-0            gsubfn_0.7             pillar_1.6.2          
# [165] htmltools_0.5.1.1      glue_1.4.2             fastmap_1.1.0          DT_0.18               
# [169] BiocParallel_1.26.1    codetools_0.2-18       utf8_1.2.2             lattice_0.20-44       
# [173] network_1.17.1         sqldf_0.4-11           curl_4.3.2             survival_3.2-11       
# [177] rmarkdown_2.9          statnet.common_4.5.0   repr_1.1.3             munsell_0.5.0         
# [181] DO.db_2.9              GenomeInfoDbData_1.2.6 iterators_1.0.13       RCy3_2.12.4           
# [185] impute_1.66.0          haven_2.4.3            gtable_0.3.0          
