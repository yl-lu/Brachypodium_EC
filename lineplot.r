#-----------------------------
# line plot WITH WITH WITH cluster info 
# in WT elf3 and phyC 
# in WT elf3 and phyC
# in WT elf3 and phyC
#-----------------------------

# cluster info come from WT_LD/WT_SD (193R/144R)

library(ggplot2)
# library(ggridges)
library(fields)
setwd("C:/project-Brachypodium/clustering")
df <- read.table("TPM_original.week3_week4.WT_LD_193RvsWT_SD_144R.cluster12.gmm.plotOrder.xls", header = T, row.names = 1)
# mycolnames <- colnames(df)

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
                "WT_ZT20_SD_144R_week3","WT_ZT22_LD_193R_week3",
                "phyC_ZT00_LD_201R_week4",
                "phyC_ZT01_LD_201R_week4",
                "phyC_ZT04_LD_201R_week4",
                "phyC_ZT08_LD_201R_week4",
                "phyC_ZT12_LD_201R_week4",
                "phyC_ZT16_LD_201R_week4",
                "phyC_ZT20_LD_201R_week4",
                "phyC_ZT22_LD_201R_week4")

# ,
# "phyC_ZT00_LD_149R_week2",
# "phyC_ZT01_LD_149R_week2",
# "phyC_ZT04_LD_149R_week2",
# "phyC_ZT08_LD_149R_week2",
# "phyC_ZT12_LD_149R_week2",
# "phyC_ZT16_LD_149R_week2",
# "phyC_ZT20_LD_149R_week2"

# df <- data.frame(log2(df[,1:18]+1),
#                  df$clustersTree)
df <- data.frame(log2(df[,mycolnames] + 1),
                 df$clustersTree)
colnames(df) <- c(mycolnames,"cluster")
head(df)

# df.select.all <- read.table("Phil_selected.geneID.geneSymbol.TPM.clusterInfo.xls",
                            # header = T, row.names = 1, sep = "\t")
# df.select.all <- read.table("Phil_selected_round2.geneID.geneSymbol.TPM.clusterInfo.xls",
#                             header = T, row.names = 1, sep = "\t")
df.select.all <- read.table("Phil_selected_round3.geneID.geneSymbol.TPM.clusterInfo.xls",
                            header = T, row.names = 1, sep = "\t")
colnames(df.select.all)

cluster <- NULL
for (k in 1:12){
  cluster[[k]] <- df[df$cluster==k,]
}

df.select.wk3 <- df.select.all[,c(mycolnames,"cluster")]

# optional
df.select.wk3 <- data.frame(log2(df.select.wk3[,mycolnames] + 1),
                            df.select.wk3$cluster)
colnames(df.select.wk3) <- c(mycolnames,"cluster")

gene <- NULL
p <- NULL

for (i in 1:length(df.select.wk3[,1])){
  gene[[i]] <- df.select.wk3[i,]
  reshapeList <- strsplit(colnames(gene[[i]]), split = '_')
  myline <- NULL
  for (j in 1:(length(reshapeList)-1)){
    myline[[j]] <- as.vector(unlist(c(reshapeList[[j]], gene[[i]][j])))
  }
  df.gene <- data.frame(matrix(unlist(myline),
                               nrow = length(reshapeList)-1,
                               byrow = TRUE),
                        stringsAsFactors = FALSE)
  colnames(df.gene) <- c("Genotype","ZT","Day_length","LibID","Week","TPM")
  df.gene$TPM <- as.numeric(df.gene$TPM)
  df.gene$cluster <- rep(gene[[i]]$cluster,length(df.gene$TPM))
  df.gene$ZT <- as.integer(gsub("ZT","",df.gene$ZT))
  df.gene$Genotype <- factor(df.gene$Genotype, levels = c("WT","elf3","phyC"))
  df.gene$Day_length <- as.factor(df.gene$Day_length)
  # df.gene$symbol <- rep(rownames(df.select.wk3[i,]),length(df.gene$TPM))
  df.gene$symbol <- rep(rownames(df.select.wk3[i,]),length(df.gene$TPM))
  
  mysymbol <- unique(df.gene$symbol)
  myclusterName <- unique(df.gene$cluster)
  
  if(myclusterName %in% df$cluster){
    mycluster <- cluster[[unique(df.gene$cluster)]]
    mysummary <- summary(mycluster)
    
    df.gene$sd <- apply(mycluster,2,sd)[1:length(colnames(mycluster))-1]
    df.gene$mean <- apply(mycluster,2,mean)[1:length(colnames(mycluster))-1]
    df.gene$median<- apply(mycluster,2,median)[1:length(colnames(mycluster))-1]
    df.gene$upperquartile <- as.numeric(gsub("1st Qu.:","",mysummary[2,][1:length(colnames(mycluster))-1]))
    df.gene$lowerquartile <- as.numeric(gsub("3rd Qu.:","",mysummary[5,][1:length(colnames(mycluster))-1]))
    
    p[[i]] <- ggplot(df.gene,
                     aes(x = ZT, y = TPM,
                         color = Day_length,
                         group = Day_length)) +
      geom_ribbon(aes(x = ZT, y = TPM,
                      ymin = median-sd, ymax = median+sd,
                      alpha = 0.02,
                      group = Day_length,
                      fill = Day_length),
                  inherit.aes = F) +
      # geom_smooth(aes(fill = Day_length),
      #             method = lm,
      #             formula = y ~ splines::bs(x, 3),
      #             se = T,
      #             size = 2,
      #             span = 0.1) +
      geom_line(size = 2) +
      scale_color_manual(values = c("red2","deepskyblue3")) +
      scale_fill_manual(values = c("lightcoral","lightskyblue")) +
      facet_wrap(~Genotype, scales = "fixed") +
      scale_x_continuous(breaks = unique(df.gene$ZT),
                         labels = as.character(unique(df.gene$ZT)),
                         position = 'bottom') +
      #ylim(0,3) + 
      theme(axis.text.x = element_text(size = 21,
                                       color = "black",
                                       face = "plain",
                                       vjust = 0,
                                       hjust = 0.5,
                                       angle = 0)) +
      theme(axis.text.y = element_text(size = 21,
                                       color = "black",
                                       face = "plain",
                                       vjust = 0.5,
                                       hjust = 1,
                                       angle = 0)) +
      theme(axis.title.x = element_text(size = 21,
                                        color = "black",
                                        face = "plain",
                                        vjust = -0.3,
                                        hjust = 0.5,
                                        angle = 0)) +
      theme(legend.title = element_text(size = 21,
                                        color = "black",
                                        face = "plain",
                                        vjust = 5,
                                        hjust = 0.5,
                                        angle = 0)) +
      theme(legend.text = element_text(size = 18,
                                       color = "black",
                                       face = "plain",
                                       angle = 0)) +
      theme_bw() +
      theme(panel.border = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")) +
      labs(x = paste(mysymbol," Cluster ",myclusterName, sep = "")) +
      # labs(x = paste(mysymbol, sep = "")) +
      labs(y = "log2(TPM+1)")
    # labs(y = "TPM")
  } else {
    p[[i]] <- ggplot(df.gene,
                     aes(x = ZT, y = TPM,
                         color = Day_length,
                         group = Day_length)) +
      # geom_ribbon(aes(x = ZT, y = TPM,
      #                 ymin = median-sd, ymax = median+sd,
      #                 alpha = 0.02,
      #                 group = Day_length,
      #                 fill = Day_length),
      #             inherit.aes = F) +
      # geom_smooth(aes(fill = Day_length),
      #             method = lm,
      #             formula = y ~ splines::bs(x, 3),
      #             se = T,
      #             size = 2,
    #             span = 0.1) +
    geom_line(size = 2) +
      scale_color_manual(values = c("red2","deepskyblue3")) +
      scale_fill_manual(values = c("lightcoral","lightskyblue")) +
      facet_wrap(~Genotype, scales = "fixed") +
      scale_x_continuous(breaks = unique(df.gene$ZT),
                         labels = as.character(unique(df.gene$ZT)),
                         position = 'bottom') +
      #ylim(0,3) + 
      theme(axis.text.x = element_text(size = 21,
                                       color = "black",
                                       face = "plain",
                                       vjust = 0,
                                       hjust = 0.5,
                                       angle = 0)) +
      theme(axis.text.y = element_text(size = 21,
                                       color = "black",
                                       face = "plain",
                                       vjust = 0.5,
                                       hjust = 1,
                                       angle = 0)) +
      theme(axis.title.x = element_text(size = 21,
                                        color = "black",
                                        face = "plain",
                                        vjust = -0.3,
                                        hjust = 0.5,
                                        angle = 0)) +
      theme(legend.title = element_text(size = 21,
                                        color = "black",
                                        face = "plain",
                                        vjust = 5,
                                        hjust = 0.5,
                                        angle = 0)) +
      theme(legend.text = element_text(size = 18,
                                       color = "black",
                                       face = "plain",
                                       angle = 0)) +
      theme_bw() +
      theme(panel.border = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")) +
      labs(x = paste(mysymbol," Cluster ",myclusterName, sep = "")) +
      # labs(x = paste(mysymbol, sep = "")) +
      labs(y = "log2(TPM+1)")
    # labs(y = "TPM")
  }
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
pdf("lineplots.week3.WT_elf3_phyC-wk4.Phil_selectedGenes.clusterBackground.pdf", width = 8, height = length(rownames(df.select.all))*(30/11))
multiplot(p[[1]],p[[2]],p[[3]],
          p[[4]],p[[5]],p[[6]],
          p[[7]],p[[8]],p[[9]],
          p[[10]],p[[11]],p[[12]],
          p[[13]],p[[14]],p[[15]],
          p[[16]],p[[17]],p[[18]],
          p[[19]],p[[20]],p[[21]],
          p[[22]],p[[23]],p[[24]],
          p[[25]],p[[26]],p[[27]],
          p[[28]],p[[29]],p[[30]],
          p[[31]],p[[32]],p[[33]],
          p[[34]],p[[35]],p[[36]],
          p[[37]],p[[38]],p[[39]],
          p[[40]],p[[41]],p[[42]],
          p[[43]],p[[44]],p[[45]],
          p[[46]],p[[47]],p[[48]],
          p[[49]],p[[50]],p[[51]],
          p[[52]],p[[53]],p[[54]],
          p[[55]],p[[56]],p[[57]],
          p[[58]],p[[59]],p[[60]],
          p[[61]],p[[62]],p[[63]],
          p[[64]],p[[65]],p[[66]],
          p[[67]],p[[68]],p[[69]],
          p[[70]],p[[71]],
          cols = 1)
dev.off()


#-----------------------------
# line plot WITH WITH WITH cluster info 
# WT_week3_LD + ppd1_week3_LD; WT_week6_SD + PPD1OX_week6_SD
# WT_week3_LD + ppd1_week3_LD; WT_week6_SD + PPD1OX_week6_SD
# WT_week3_LD + ppd1_week3_LD; WT_week6_SD + PPD1OX_week6_SD
#-----------------------------

# cluster info come from week3 WT_LD/WT_SD (193R/144R)

library(ggplot2)
# library(ggridges)
library(fields)
setwd("C:/project-Brachypodium/clustering")
df <- read.table("TPM_original.week3_week4.WT_LD_193RvsWT_SD_144R.cluster12.gmm.plotOrder.xls", header = T, row.names = 1)
df.all <- read.table("TPM_all.txt", header = T, row.names = 1)
# mycolnames <- colnames(df)
# colnames(df.all)

mycolnames <- c(grep("193R_week3|199R_week6", grep("WT", colnames(df.all), value = T), value = T),
                grep("ppd1", colnames(df.all), value = T),
                grep("PPD1", colnames(df.all), value = T))

# "WT_ZT00_SD_144R_week3",
# "WT_ZT04_SD_144R_week3",
# "WT_ZT08_SD_144R_week3",
# "WT_ZT12_SD_144R_week3",
# "WT_ZT16_SD_144R_week3",
# "WT_ZT20_SD_144R_week3",

# ,
# "phyC_ZT00_LD_149R_week2",
# "phyC_ZT01_LD_149R_week2",
# "phyC_ZT04_LD_149R_week2",
# "phyC_ZT08_LD_149R_week2",
# "phyC_ZT12_LD_149R_week2",
# "phyC_ZT16_LD_149R_week2",
# "phyC_ZT20_LD_149R_week2"

# df <- data.frame(log2(df[,1:18]+1),
#                  df$clustersTree)
df <- data.frame(log2(df[,mycolnames] + 1),
                 df$clustersTree)
colnames(df) <- c(mycolnames,"cluster")
head(df)

# df.select.all <- read.table("select_all.geneID.geneSymbol.TPM.clusterInfo.xls",
#                             header = T, row.names = 1, sep = "\t")
df.select.all <- read.table("select_all.geneID.geneSymbol.TPM.round2.clusterInfo.xls",
                            header = T, row.names = 1, sep = "\t")
head(df.select.all)

cluster <- NULL
for (k in 1:12){
  cluster[[k]] <- df[df$cluster==k,]
}

df.select.wk3 <- df.select.all[,c(mycolnames,"cluster")]

write.table(df.select.wk3,
            file = paste("TPM_table.week3_LD_WT_ppd1.week6_SD_WT_PPD1OX.selectedGenes_all.round2.clusterBackground.xls",sep = ""),
            row.names = T, col.names = T, quote = FALSE, sep = '\t')

# optional
df.select.wk3 <- data.frame(log2(df.select.wk3[,mycolnames] + 1),
                            df.select.wk3$cluster)
colnames(df.select.wk3) <- c(mycolnames,"cluster")

gene <- NULL
p <- NULL

for (i in 1:length(df.select.wk3[,1])){
  gene[[i]] <- df.select.wk3[i,]
  reshapeList <- strsplit(colnames(gene[[i]]), split = '_')
  myline <- NULL
  for (j in 1:(length(reshapeList)-1)){
    myline[[j]] <- as.vector(unlist(c(reshapeList[[j]], gene[[i]][j])))
  }
  df.gene <- data.frame(matrix(unlist(myline),
                               nrow = length(reshapeList)-1,
                               byrow = TRUE),
                        stringsAsFactors = FALSE)
  colnames(df.gene) <- c("Genotype","ZT","Day_length","LibID","Week","TPM")
  df.gene$TPM <- as.numeric(df.gene$TPM)
  df.gene$cluster <- rep(gene[[i]]$cluster,length(df.gene$TPM))
  df.gene$ZT <- as.integer(gsub("ZT","",df.gene$ZT))
  df.gene$Genotype <- factor(df.gene$Genotype, levels = c("WT","ppd1","PPD1OX"))
  df.gene$LibID <- factor(df.gene$LibID, levels = c("193R","196R","199R"))
  df.gene$Day_length <- as.factor(df.gene$Day_length)
  df.gene$symbol <- rep(rownames(df.select.wk3[i,]),length(df.gene$TPM))
  
  mysymbol <- unique(df.gene$symbol)
  myclusterName <- unique(df.gene$cluster)
  
  if(myclusterName %in% df$cluster){
    mycluster <- cluster[[unique(df.gene$cluster)]]
    mysummary <- summary(mycluster)
    
    df.gene$sd <- apply(mycluster,2,sd)[1:length(colnames(mycluster))-1]
    df.gene$mean <- apply(mycluster,2,mean)[1:length(colnames(mycluster))-1]
    df.gene$median <- apply(mycluster,2,median)[1:length(colnames(mycluster))-1]
    df.gene$upperquartile <- as.numeric(gsub("1st Qu.:","",mysummary[2,][1:length(colnames(mycluster))-1]))
    df.gene$lowerquartile <- as.numeric(gsub("3rd Qu.:","",mysummary[5,][1:length(colnames(mycluster))-1]))
    
    # df.gene <- subset(df.gene, df.gene$ZT != 1)
    
    p[[i]] <- ggplot(df.gene,
                     aes(x = ZT, y = TPM,
                         group = Genotype,
                         color = Genotype)) +
      geom_ribbon(aes(ymin = median-sd, ymax = median+sd,
                      group = Genotype,
                      fill = Genotype,
                      alpha = 0.02),
                  linetype = 0) +
      geom_line(size = 2) +
      scale_color_manual(values = c("black","orangered1","olivedrab3")) +
      scale_fill_manual(values = c("grey73","peachpuff","#d4fab2")) +
      facet_wrap(~Day_length, scales = "fixed") +
      scale_x_continuous(breaks = unique(df.gene$ZT),
                         labels = as.character(unique(df.gene$ZT)),
                         position = 'bottom') +
      #ylim(0,3) + 
      theme(axis.text.x = element_text(size = 21,
                                       color = "black",
                                       face = "plain",
                                       vjust = 0,
                                       hjust = 0.5,
                                       angle = 0)) +
      theme(axis.text.y = element_text(size = 21,
                                       color = "black",
                                       face = "plain",
                                       vjust = 0.5,
                                       hjust = 1,
                                       angle = 0)) +
      theme(axis.title.x = element_text(size = 21,
                                        color = "black",
                                        face = "plain",
                                        vjust = -0.3,
                                        hjust = 0.5,
                                        angle = 0)) +
      theme(legend.title = element_text(size = 21,
                                        color = "black",
                                        face = "plain",
                                        vjust = 5,
                                        hjust = 0.5,
                                        angle = 0)) +
      theme(legend.text = element_text(size = 18,
                                       color = "black",
                                       face = "plain",
                                       angle = 0)) +
      theme_bw() +
      theme(panel.border = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")) +
      labs(x = paste(mysymbol," Cluster ",myclusterName, sep = "")) +
      # labs(x = paste(mysymbol, sep = "")) +
      labs(y = "log2(TPM+1)")
    # labs(y = "TPM")
  } else {
    p[[i]] <- ggplot(df.gene,
                     aes(x = ZT, y = TPM,
                         # group = Day_length,
                         color = Genotype)) +
      # geom_ribbon(aes(x = ZT, y = TPM,
      #                 ymin = median-sd, ymax = median+sd,
      #                 alpha = 0.02,
      #                 group = Day_length,
      #                 fill = Day_length),
      #             inherit.aes = F) +
      # geom_smooth(aes(fill = Day_length),
      #             method = lm,
      #             formula = y ~ splines::bs(x, 3),
      #             se = T,
      #             size = 2,
    #             span = 0.1) +
    geom_line(size = 2) +
      scale_color_manual(values = c("black","orangered1","olivedrab3")) +
      scale_fill_manual(values = c("grey73","peachpuff","palegreen")) +
      facet_wrap(~Day_length, scales = "fixed") +
      scale_x_continuous(breaks = unique(df.gene$ZT),
                         labels = as.character(unique(df.gene$ZT)),
                         position = 'bottom') +
      #ylim(0,3) + 
      theme(axis.text.x = element_text(size = 21,
                                       color = "black",
                                       face = "plain",
                                       vjust = 0,
                                       hjust = 0.5,
                                       angle = 0)) +
      theme(axis.text.y = element_text(size = 21,
                                       color = "black",
                                       face = "plain",
                                       vjust = 0.5,
                                       hjust = 1,
                                       angle = 0)) +
      theme(axis.title.x = element_text(size = 21,
                                        color = "black",
                                        face = "plain",
                                        vjust = -0.3,
                                        hjust = 0.5,
                                        angle = 0)) +
      theme(legend.title = element_text(size = 21,
                                        color = "black",
                                        face = "plain",
                                        vjust = 5,
                                        hjust = 0.5,
                                        angle = 0)) +
      theme(legend.text = element_text(size = 18,
                                       color = "black",
                                       face = "plain",
                                       angle = 0)) +
      theme_bw() +
      theme(panel.border = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")) +
      labs(x = paste(mysymbol," Cluster ",myclusterName, sep = "")) +
      # labs(x = paste(mysymbol, sep = "")) +
      labs(y = "log2(TPM+1)")
    # labs(y = "TPM")
  }
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
pdf("lineplots.week3_LD_WT_ppd1.week6_SD_WT_PPD1OX.selectedGenes_all.round2.clusterBackground.pdf", width = 6, height = length(rownames(df.select.all))*(30/11))
multiplot(p[[1]],p[[2]],p[[3]],
          p[[4]],p[[5]],p[[6]],
          p[[7]],p[[8]],p[[9]],
          p[[10]],p[[11]],p[[12]],
          p[[13]],p[[14]],p[[15]],
          p[[16]],p[[17]],p[[18]],
          p[[19]],p[[20]],p[[21]],
          p[[22]],p[[23]],p[[24]],
          p[[25]],p[[26]],p[[27]],
          p[[28]],p[[29]],p[[30]],
          p[[31]],p[[32]],p[[33]],
          p[[34]],p[[35]],p[[36]],
          p[[37]],p[[38]],p[[39]],
          p[[40]],p[[41]],p[[42]],
          p[[43]],p[[44]],p[[45]],
          p[[46]],p[[47]],p[[48]],
          p[[49]],p[[50]],p[[51]],
          p[[52]],p[[53]],p[[54]],
          p[[55]],p[[56]],p[[57]],
          p[[58]],p[[59]],p[[60]],
          p[[61]],p[[62]],p[[63]],
          p[[64]],p[[65]],p[[66]],
          p[[67]],p[[68]],p[[69]],
          p[[70]],p[[71]],
          cols = 1)
dev.off()