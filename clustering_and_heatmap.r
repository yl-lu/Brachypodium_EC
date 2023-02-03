============
# WT_LD/WT_SD 193R/144R # cluster reference
# elf3_LD/WT_LD 193R/193R
# elf3_SD/WT_SD 144R/144R
# elf3_OX_LD/WT_LD 169R/169R
# phyc_LD/WT_LD 201R/201R

library(ClustVarLV)
library(ClusterR)
library(pheatmap)
library(grid)
library(amap)
library(RColorBrewer)

setwd("C:/project-Brachypodium/clustering")

df.mrg <- read.table("WT_LD_193R.WT_SD_144R.elf3_SD_144R.elf3OX_LD_169R.WT_LD_169R.phyC_LD.xls",
                      header = T,
                      row.names = 1,
                      sep = "\t")
colnames(df.mrg)

wt <- df.mrg[,c("WT_ZT00_LD_193R",
                 "WT_ZT04_LD_193R",
                 "WT_ZT08_LD_193R",
                 "WT_ZT12_LD_193R",
                 "WT_ZT16_LD_193R",
                 "WT_ZT20_LD_193R",
                 "WT_ZT00_SD_144R",
                 "WT_ZT04_SD_144R",
                 "WT_ZT08_SD_144R",
                 "WT_ZT12_SD_144R",
                 "WT_ZT16_SD_144R",
                 "WT_ZT20_SD_144R")]

# write.table(wt,
#             file = paste("WT_LD_193R.WT_SD_144R.xls",sep = ""),
#             row.names = T, col.names = T, quote = FALSE, sep = '\t')

wt <- log2(wt + 1)
wt <- wt[apply(wt,1,max)>2, ]
WT_LDvsSD <- data.frame(WT_ZT00_LDvsSD = wt$WT_ZT00_LD_193R - wt$WT_ZT00_SD_144R,
                        WT_ZT04_LDvsSD = wt$WT_ZT04_LD_193R - wt$WT_ZT04_SD_144R,
                        WT_ZT08_LDvsSD = wt$WT_ZT08_LD_193R - wt$WT_ZT08_SD_144R,
                        WT_ZT12_LDvsSD = wt$WT_ZT12_LD_193R - wt$WT_ZT12_SD_144R,
                        WT_ZT16_LDvsSD = wt$WT_ZT16_LD_193R - wt$WT_ZT16_SD_144R,
                        WT_ZT20_LDvsSD = wt$WT_ZT20_LD_193R - wt$WT_ZT20_SD_144R)
rownames(WT_LDvsSD) <- rownames(wt)

WT_LDvsSD[WT_LDvsSD > 2] <- 2
WT_LDvsSD[WT_LDvsSD < -2] <- -2
data <- WT_LDvsSD
opt_gmm <- Optimal_Clusters_GMM(data,
                                max_clusters = 30,
                                criterion = "BIC",
                                dist_mode = "eucl_dist",
                                seed_mode = "static_spread",
                                km_iter = 30,
                                em_iter = 30,
                                var_floor = 1e-10,
                                plot_data = T)
for (clusterN in c(12)){
  gmm <- GMM(data,
             gaussian_comps = clusterN,
             dist_mode = "eucl_dist",
             seed_mode = "static_spread",
             km_iter = 30,
             em_iter = 30,
             verbose = F)
  pr <- predict_GMM(data,
                    gmm$centroids,
                    gmm$covariance_matrices,
                    gmm$weights)
  names(pr$cluster_labels) <- rownames(data)
  clustersTree <- sort(pr$cluster_labels)
  gaps <- which((clustersTree[-1] - clustersTree[-length(clustersTree)]) != 0)
  #---hclust (optional use)--#
  # d <- Dist(data, method = "maximum")
  # tree <- hclust(d, method = "ward.D")
  # clustersTree <- cutree(tree,clusterN)[tree$order]
  # clustersTree <- sort(clustersTree)
  # gaps <- which((clustersTree[-1] - clustersTree[-length(clustersTree)]) !=0)
  #--------------------------#
  gene.cluster <- as.data.frame(clustersTree)
  df.cluster <- cbind(gene.cluster,data[rownames(gene.cluster),])
  data_clusters <- cbind(gene.cluster,data[rownames(gene.cluster),])
  write.table(data_clusters,
              file = paste("df.mrg.WT_LD_193RvsWT_SD_144R.cluster",clusterN,".gmm.plotOrder.xls",sep = ""),
              row.names = T, col.names = T, quote = FALSE, sep = '\t')
  
  my_gene_clusters <- data.frame(cluster = df.cluster[,1])
  rownames(my_gene_clusters) <- rownames(df.cluster)
  colorsVec <- c(brewer.pal(12, "Set3"),
                 brewer.pal(8, "Accent"))
  names(colorsVec) <- c(1:clusterN)
  my_colors <- list(cluster = colorsVec)
  
  # # WT clusters as ref, other samples follow the order of WT clusters
  # wk3 <- read.table("TPM_week3_169Rbackward.txt",
  #                   header = T,
  #                   row.names = 1,
  #                   sep = "\t")
  # wk3 <- read.table("WT_LD_193R.WT_SD_144R.elf3_SD_144R.elf3OX_LD_169R.WT_LD_169R.phyC_LD.xls",
  #                   header = T,
  #                   row.names = 1,
  #                   sep = "\t")
  # wk3_wk4$elf3_OX_ZT00_LD_169R <- wk3$ELF3_OX_ZT00_LD
  # wk3_wk4$elf3_OX_ZT04_LD_169R <- wk3$ELF3_OX_ZT04_LD
  # wk3_wk4$elf3_OX_ZT12_LD_169R <- wk3$ELF3_OX_ZT12_LD
  # wk3_wk4$elf3_OX_ZT20_LD_169R <- wk3$ELF3_OX_ZT20_LD
  # wk3_wk4$WT_ZT00_LD_169R <- wk3$WT_ZT00_LD
  # wk3_wk4$WT_ZT04_LD_169R <- wk3$WT_ZT04_LD
  # wk3_wk4$WT_ZT12_LD_169R <- wk3$WT_ZT12_LD
  # wk3_wk4$WT_ZT20_LD_169R <- wk3$WT_ZT20_LD
  
  df.show.df.mrg <- df.mrg[,c("WT_ZT00_LD_193R",
                                "WT_ZT04_LD_193R",
                                "WT_ZT08_LD_193R",
                                "WT_ZT12_LD_193R",
                                "WT_ZT16_LD_193R",
                                "WT_ZT20_LD_193R",
                                "elf3_ZT00_SD_144R",
                                "elf3_ZT04_SD_144R",
                                "elf3_ZT08_SD_144R",
                                "elf3_ZT12_SD_144R",
                                "elf3_ZT16_SD_144R",
                                "elf3_ZT20_SD_144R",
                                # "elf3_ZT00_LD",
                                # "elf3_ZT04_LD",
                                # "elf3_ZT08_LD",
                                # "elf3_ZT12_LD",
                                # "elf3_ZT16_LD",
                                # "elf3_ZT20_LD",
                                "WT_ZT00_SD_144R",
                                "WT_ZT04_SD_144R",
                                "WT_ZT08_SD_144R",
                                "WT_ZT12_SD_144R",
                                "WT_ZT16_SD_144R",
                                "WT_ZT20_SD_144R",
                                "elf3_OX_ZT00_LD_169R",
                                "elf3_OX_ZT04_LD_169R",
                                "elf3_OX_ZT12_LD_169R",
                                "elf3_OX_ZT20_LD_169R",
                                "WT_ZT00_LD_169R",
                                "WT_ZT04_LD_169R",
                                "WT_ZT12_LD_169R",
                                "WT_ZT20_LD_169R",
                                "phyC_ZT00_LD",
                                "phyC_ZT04_LD",
                                "phyC_ZT08_LD",
                                "phyC_ZT12_LD",
                                "phyC_ZT16_LD",
                                "phyC_ZT20_LD")]
  
  # write.table(df.show.wk3_wk4,
  #             file = paste("WT_LD_193R.WT_SD_144R.elf3_SD_144R.elf3OX_LD_169R.WT_LD_169R.phyC_LD.xls",sep = ""),
  #             row.names = T, col.names = T, quote = FALSE, sep = '\t')
  
  head(df.show.df.mrg)
  head(df.cluster)
  # length(df.cluster$clustersTree[df.cluster$clustersTree==2])
  df.show.df.mrg <- log2(df.show.df.mrg + 1)
  
  # wk3_wk4.elf3_LDvsWT_LD <- data.frame(
  #   wk3_wk4.elf3_LDvsWT_LD_ZT00 = df.show.wk3_wk4$elf3_ZT00_LD - df.show.wk3_wk4$WT_ZT00_LD_193R,
  #   wk3_wk4.elf3_LDvsWT_LD_ZT04 = df.show.wk3_wk4$elf3_ZT04_LD - df.show.wk3_wk4$WT_ZT04_LD_193R,
  #   wk3_wk4.elf3_LDvsWT_LD_ZT08 = df.show.wk3_wk4$elf3_ZT08_LD - df.show.wk3_wk4$WT_ZT08_LD_193R,
  #   wk3_wk4.elf3_LDvsWT_LD_ZT12 = df.show.wk3_wk4$elf3_ZT12_LD - df.show.wk3_wk4$WT_ZT12_LD_193R,
  #   wk3_wk4.elf3_LDvsWT_LD_ZT16 = df.show.wk3_wk4$elf3_ZT16_LD - df.show.wk3_wk4$WT_ZT16_LD_193R,
  #   wk3_wk4.elf3_LDvsWT_LD_ZT20 = df.show.wk3_wk4$elf3_ZT20_LD - df.show.wk3_wk4$WT_ZT20_LD_193R
  # )
  
  df.mrg.elf3_SDvsWT_SD <- data.frame(
    df.mrg.elf3_SDvsWT_SD_ZT00 = df.show.df.mrg$elf3_ZT00_SD_144R - df.show.df.mrg$WT_ZT00_SD_144R,
    df.mrg.elf3_SDvsWT_SD_ZT04 = df.show.df.mrg$elf3_ZT04_SD_144R - df.show.df.mrg$WT_ZT04_SD_144R,
    df.mrg.elf3_SDvsWT_SD_ZT08 = df.show.df.mrg$elf3_ZT08_SD_144R - df.show.df.mrg$WT_ZT08_SD_144R,
    df.mrg.elf3_SDvsWT_SD_ZT12 = df.show.df.mrg$elf3_ZT12_SD_144R - df.show.df.mrg$WT_ZT12_SD_144R,
    df.mrg.elf3_SDvsWT_SD_ZT16 = df.show.df.mrg$elf3_ZT16_SD_144R - df.show.df.mrg$WT_ZT16_SD_144R,
    df.mrg.elf3_SDvsWT_SD_ZT20 = df.show.df.mrg$elf3_ZT20_SD_144R - df.show.df.mrg$WT_ZT20_SD_144R
  )
  
  df.mrg.elf3_OX_LDvsWT_LD <- data.frame(
    df.mrg.elf3_OX_LDvsWT_LD_ZT00 = df.show.df.mrg$elf3_OX_ZT00_LD_169R - df.show.df.mrg$WT_ZT00_LD_169R,
    df.mrg.elf3_OX_LDvsWT_LD_ZT04 = df.show.df.mrg$elf3_OX_ZT04_LD_169R - df.show.df.mrg$WT_ZT04_LD_169R,
    df.mrg.elf3_OX_LDvsWT_LD_ZT12 = df.show.df.mrg$elf3_OX_ZT12_LD_169R - df.show.df.mrg$WT_ZT12_LD_169R,
    df.mrg.elf3_OX_LDvsWT_LD_ZT20 = df.show.df.mrg$elf3_OX_ZT20_LD_169R - df.show.df.mrg$WT_ZT20_LD_169R
  )
  
  df.mrg.phyc_LDvsWT_LD <- data.frame(
    df.mrg.phyc_LDvsWT_LD_ZT00 = df.show.df.mrg$phyC_ZT00_LD - df.show.df.mrg$WT_ZT00_LD_193R,
    df.mrg.phyc_LDvsWT_LD_ZT04 = df.show.df.mrg$phyC_ZT04_LD - df.show.df.mrg$WT_ZT04_LD_193R,
    df.mrg.phyc_LDvsWT_LD_ZT08 = df.show.df.mrg$phyC_ZT08_LD - df.show.df.mrg$WT_ZT08_LD_193R,
    df.mrg.phyc_LDvsWT_LD_ZT12 = df.show.df.mrg$phyC_ZT12_LD - df.show.df.mrg$WT_ZT12_LD_193R,
    df.mrg.phyc_LDvsWT_LD_ZT16 = df.show.df.mrg$phyC_ZT16_LD - df.show.df.mrg$WT_ZT16_LD_193R,
    df.mrg.phyc_LDvsWT_LD_ZT20 = df.show.df.mrg$phyC_ZT20_LD - df.show.df.mrg$WT_ZT20_LD_193R
  )
  
  # rownames(wk3_wk4.elf3_LDvsWT_LD) <- rownames(df.show.wk3_wk4)
  rownames(df.mrg.elf3_SDvsWT_SD) <- rownames(df.show.df.mrg)
  rownames(df.mrg.elf3_OX_LDvsWT_LD) <- rownames(df.show.df.mrg)
  rownames(df.mrg.phyc_LDvsWT_LD) <- rownames(df.show.df.mrg)
  
  
  # colnames(wk3_wk4.elf3_LDvsWT_LD)
  colnames(df.mrg.elf3_SDvsWT_SD)
  colnames(df.mrg.elf3_OX_LDvsWT_LD)
  colnames(df.mrg.phyc_LDvsWT_LD)
  
  
  # wk3_wk4.elf3_LDvsWT_LD <- wk3_wk4.elf3_LDvsWT_LD[rownames(df.cluster),]
  df.mrg.elf3_SDvsWT_SD <- df.mrg.elf3_SDvsWT_SD[rownames(df.cluster),]
  df.mrg.elf3_OX_LDvsWT_LD <- df.mrg.elf3_OX_LDvsWT_LD[rownames(df.cluster),]
  df.mrg.phyc_LDvsWT_LD <- df.mrg.phyc_LDvsWT_LD[rownames(df.cluster),]
  
  
  # df.cluster is the WT_LD - WT_SD and clustered
  df.cluster.allsamples <- cbind(df.cluster,
                                 # df.mrg.elf3_LDvsWT_LD,
                                 df.mrg.elf3_SDvsWT_SD,
                                 df.mrg.elf3_OX_LDvsWT_LD,
                                 df.mrg.phyc_LDvsWT_LD)
  
  my_gene_clusters <- data.frame(cluster = df.cluster[,1])
  rownames(my_gene_clusters) <- rownames(df.cluster)
  
  #apply(df.cluster.allsamples,1,)>2
  df.cluster.allsamples[df.cluster.allsamples > 2] <- 2
  df.cluster.allsamples[df.cluster.allsamples < -2] <- -2
  
  colnames(df.cluster.allsamples[,-1])
  
  my_gaps_col <- c(6,12,16)
  
  pdf(paste("df.mrg.pannel2.1.WT_LD_193R_vs_WT_SD_144R.elf3LD_deleted.cluster",clusterN,".gmm.pdf",sep = ""),
      height = 11,
      width = 6)
  phm <- pheatmap(#df.cluster[,-1],
    df.cluster.allsamples[,-1],
    color = colorRampPalette(c("#c1207e","snow","#669900"))(200),
    border_color = "snow",
    cluster_rows = F,
    cluster_cols = F,                
    treeheight_row = 30, 
    treeheight_col = 30,
    fontsize_col = 8,
    fontsize_row = 10,
    cutree_rows = 12,
    angle_col = 90,
    gaps_row = gaps,
    gaps_col = my_gaps_col,
    annotation_colors = my_colors,
    annotation_row = my_gene_clusters,
    annotation_legend = T,
    labels_col = c(rep(c("ZT0","ZT4","ZT8","ZT12",
                         "ZT16","ZT20"),2),
                   rep(c("ZT0","ZT4","ZT12","ZT20"),1),
                   rep(c("ZT0","ZT4","ZT8","ZT12",
                         "ZT16","ZT20"),1)),
    annotation_names_row = T,                  
    annotation_names_col = T,
    show_rownames = F,
    display_numbers = F,
    drop_levels = T
  )
  dev.off()
} 


#============================================================================
# 2022-09-01
# Clutering by using ELF3 dependent genes (elf3SD vs WTSD) and marked the ELF3 bound genes
# define up-regulated gene in elf3 mutant vs WT

mycolnames <- c("elf3_ZT00_LD_193R_week3","elf3_ZT00_SD_144R_week3",
                "elf3_ZT04_LD_193R_week3","elf3_ZT04_SD_144R_week3",
                "elf3_ZT08_LD_193R_week3","elf3_ZT08_SD_144R_week3",
                "elf3_ZT12_LD_193R_week3","elf3_ZT12_SD_144R_week3",
                "elf3_ZT16_LD_193R_week3","elf3_ZT16_SD_144R_week3",
                "elf3_ZT20_LD_193R_week3","elf3_ZT20_SD_144R_week3",
                "elf3_ZT22_LD_193R_week3","WT_ZT00_LD_193R_week3",
                "WT_ZT00_SD_144R_week3","WT_ZT04_LD_193R_week3",
                "WT_ZT04_SD_144R_week3","WT_ZT08_LD_193R_week3",
                "WT_ZT08_SD_144R_week3","WT_ZT12_LD_193R_week3",
                "WT_ZT12_SD_144R_week3","WT_ZT16_LD_193R_week3",
                "WT_ZT16_SD_144R_week3","WT_ZT20_LD_193R_week3",
                "WT_ZT20_SD_144R_week3","WT_ZT22_LD_193R_week3")

setwd("C:/project-Brachypodium/clustering")

tpmAll <- read.table("TPM_all.txt",
                     header = T,
                     row.names = 1,
                     sep = "\t")
colnames(tpmAll)
wk3 <- tpmAll[,mycolnames]

wk3_SD_wt_elf3 <- wk3[,mycolnames[grep("SD",mycolnames)]]

wk3_SD_wt_elf3 <- log2(wk3_SD_wt_elf3 + 1)
wk3_SD_wt_elf3 <- wk3_SD_wt_elf3[apply(wk3_SD_wt_elf3,1,max)>2, ]

wk3_elf3SD_vs_WTSD <- data.frame(
  wk3.elf3_SDvsWT_SD_ZT00 = wk3_SD_wt_elf3$elf3_ZT00_SD - wk3_SD_wt_elf3$WT_ZT00_SD,
  wk3.elf3_SDvsWT_SD_ZT04 = wk3_SD_wt_elf3$elf3_ZT04_SD - wk3_SD_wt_elf3$WT_ZT04_SD,
  wk3.elf3_SDvsWT_SD_ZT08 = wk3_SD_wt_elf3$elf3_ZT08_SD - wk3_SD_wt_elf3$WT_ZT08_SD,
  wk3.elf3_SDvsWT_SD_ZT12 = wk3_SD_wt_elf3$elf3_ZT12_SD - wk3_SD_wt_elf3$WT_ZT12_SD,
  wk3.elf3_SDvsWT_SD_ZT16 = wk3_SD_wt_elf3$elf3_ZT16_SD - wk3_SD_wt_elf3$WT_ZT16_SD,
  wk3.elf3_SDvsWT_SD_ZT20 = wk3_SD_wt_elf3$elf3_ZT20_SD - wk3_SD_wt_elf3$WT_ZT20_SD)
rownames(wk3_elf3SD_vs_WTSD) <- rownames(wk3_SD_wt_elf3)

# select ELF3 dependent genes by setting cutoff
library(ClustVarLV)
library(ClusterR)
library(pheatmap)
library(grid)
library(amap)
library(RColorBrewer)

ELF3_dependent_SD <- as.data.frame(
  matrix(
    nrow = length(rownames(wk3_elf3SD_vs_WTSD)),
    ncol = length(colnames(wk3_elf3SD_vs_WTSD))
  )
)

for(i in 1:length(rownames(wk3_elf3SD_vs_WTSD))){
  if(sort(as.matrix(wk3_elf3SD_vs_WTSD[i,]), decreasing = T)[2] > 1){
    ELF3_dependent_SD[i,] <- wk3_elf3SD_vs_WTSD[i,]
  }
}
rownames(ELF3_dependent_SD) <- rownames(wk3_elf3SD_vs_WTSD)
colnames(ELF3_dependent_SD) <- colnames(wk3_elf3SD_vs_WTSD)
ELF3_dependent_SD <- na.omit(ELF3_dependent_SD)

dim(ELF3_dependent_SD)

write.table(ELF3_dependent_SD,
            file = paste("ELF3_dependent_genes_week3_elf3SD_vs_WTSD.xls",sep = ""),
            row.names = T, col.names = T, quote = FALSE, sep = '\t')

ELF3_dependent_SD[ELF3_dependent_SD > 2] <- 2
ELF3_dependent_SD[ELF3_dependent_SD < -2] <- -2
data_SD <- ELF3_dependent_SD

data.list <- list(SD = data_SD)

for (photoperiod in c("SD")){
  data <- data.list[[photoperiod]]
  
  for (clusterN in c(6)){
    gmm <- GMM(data,
               gaussian_comps = clusterN,
               dist_mode = "eucl_dist",
               seed_mode = "static_spread",
               km_iter = 30,
               em_iter = 30,
               verbose = F)
    pr <- predict_GMM(data,
                      gmm$centroids,
                      gmm$covariance_matrices,
                      gmm$weights)
    names(pr$cluster_labels) <- rownames(data)
    clustersTree <- sort(pr$cluster_labels)
    gaps <- which((clustersTree[-1] - clustersTree[-length(clustersTree)]) != 0)

    gene.cluster <- as.data.frame(clustersTree)
    df.cluster <- cbind(gene.cluster,data[rownames(gene.cluster),])
    data_clusters <- cbind(gene.cluster,data[rownames(gene.cluster),])
    write.table(data_clusters,
                file = paste("ELF3_dependent_genes_week3.elf3",photoperiod,"_vs_WT",photoperiod,".cluster",clusterN,".gmm.plotOrder.xls",sep = ""),
                row.names = T, col.names = T, quote = FALSE, sep = '\t')
    my_gene_clusters <- data.frame(cluster = df.cluster[,1])
    rownames(my_gene_clusters) <- rownames(df.cluster)
    colorsVec <- c(brewer.pal(12, "Set3"),
                   brewer.pal(8, "Accent"))
    names(colorsVec) <- c(1:clusterN)
    my_colors <- list(cluster = colorsVec)
    
    if(stringr::str_detect(colnames(data_clusters[,-1]), "LD")[1]){
      my_gaps_col <- c(7)
      my_labels_col <- c(rep(c("ZT0","ZT4","ZT8","ZT12",
                               "ZT16","ZT20","ZT22"),1))
    } else if(stringr::str_detect(colnames(data_clusters[,-1]), "SD")[1]){
      my_gaps_col <- c(6)
      my_labels_col <- c(rep(c("ZT0","ZT4","ZT8","ZT12",
                               "ZT16","ZT20"),1))
    }
    pdf(paste("ELF3_dependent_genes_week3.elf3",photoperiod,"_vs_WT",photoperiod,".cluster",clusterN,".gmm.pdf",sep = ""), height = 11, width = 3)
    phm <- pheatmap(df.cluster[,-1],
                    # df.cluster.allsamples[,-1],
                    color = colorRampPalette(c("#c1207e","snow","#669900"))(200),
                    border_color = "snow",
                    cluster_rows = F,
                    cluster_cols = F,                
                    treeheight_row = 30, 
                    treeheight_col = 30,
                    fontsize_col = 8,
                    fontsize_row = 10,
                    cutree_rows = 12,
                    angle_col = 90,
                    gaps_row = gaps,
                    gaps_col = my_gaps_col,
                    annotation_colors = my_colors,
                    annotation_row = my_gene_clusters,
                    annotation_legend = T,
                    labels_col = my_labels_col,
                    annotation_names_row = T,                  
                    annotation_names_col = T,
                    show_rownames = F,
                    display_numbers = F,
                    drop_levels = T
    )
    dev.off()
  }
}

# mark out ELF3 bound genes
ELF3BoundGenes <- read.table("C:/project-Brachypodium/macs2Output/ELF3_OX_ZT20_SD_rep2.targetGenes.bed", header = F, sep = "\t")
rownames(ELF3BoundGenes) <- ELF3BoundGenes$V4

head(ELF3BoundGenes)
head(df.cluster)

colorsVec <- brewer.pal(8, "Set3")
names(colorsVec) <- c(1:8)
my_colors <- list(cluster = colorsVec)

df.cluster$ELF3_bound[rownames(df.cluster) %in% rownames(ELF3BoundGenes)] <- 1
df.cluster$ELF3_bound[!(rownames(df.cluster) %in% rownames(ELF3BoundGenes))] <- 0


write.table(df.cluster,
            file = paste("ELF3_dependent_genes.elf3SD_vs_WTSD.cluster6.gmm.plotOrder.ELF3_bound.xls",sep = ""),
            row.names = T, col.names = T, quote = FALSE, sep = '\t')

df.cluster.sort <- df.cluster[order(df.cluster$clustersTree,df.cluster$ELF3_bound,decreasing = c(F,F)),]

write.table(df.cluster.sort,
            file = paste("ELF3_dependent_genes.elf3SD_vs_WTSD.cluster8.gmm.plotOrder.sorted_by_ELF3_bound.xls",sep = ""),
            row.names = T, col.names = T, quote = FALSE, sep = '\t')

my_gaps_col <- c(1,4)

my_gene_clusters <- data.frame(cluster = df.cluster.sort[,1],
                               bound = df.cluster.sort[,8])
rownames(my_gene_clusters) <- rownames(df.cluster.sort)

clusterColors <- c(rep(c(brewer.pal(12, "Set3")[9]),8))
names(clusterColors) <- c(1:clusterN)
# boundColors <- c("grey60","indianred2")
boundColors <- c(brewer.pal(12, "Set3")[c(9,4)])
names(boundColors) <- c(0,1)
my_colors <- list(cluster = clusterColors,
                  bound = boundColors)

pdf(paste("ELF3_dependent_genes.elf3SD_vs_WTSD.cluster8.gmm.plotOrder.ELF3_bound",clusterN,".gmm.pdf",sep = ""),
    height = 8,
    width = 3)
phm <- pheatmap(df.cluster.sort[,c(2:7)],
                color = colorRampPalette(c("#c1207e","snow","#669900"))(200),
                border_color = "snow",
                # scale = "row",
                cluster_rows = F,
                cluster_cols = F,                
                treeheight_row = 30, 
                treeheight_col = 30,
                fontsize_col = 8,
                fontsize_row = 10,
                cutree_rows = 12,
                angle_col = 90,
                gaps_row = gaps,
                gaps_col = my_gaps_col,
                annotation_colors = my_colors,
                annotation_row = my_gene_clusters,
                annotation_legend = T,
                labels_col = gsub("wk3.elf3_SDvsWT_SD_","",colnames(df.cluster.sort[,c(2:7)])),
                main = paste("n = ",length(rownames(df.cluster.sort)), sep = ""),
                annotation_names_row = T,                  
                annotation_names_col = T,
                show_rownames = F,
                display_numbers = F,
                drop_levels = T
)
dev.off()
