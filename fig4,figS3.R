library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

#fig4a
#refer to fig1c

#fig4b
#refer to fig1f

#fig4d
library(ggplot2)
p <- DotPlot(data.combined_fib, features = rev(gene),
             group.by = "fincell", split.by = "group")+coord_flip()
exp <- p$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)
exp$id <- as.factor(exp$id)
exp$id <- fct_inorder(exp$id)
ggplot(exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_point(aes(size=`pct.exp`,color=`avg.exp.scaled`),
             shape=21,color="black",stroke=1)+
  theme(panel.background =element_blank(),
        axis.line=element_line(colour="black", size = 1),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=11,color="black",angle=90), 
        axis.text.y=element_text(size=11,color="black"))+
  scale_color_gradientn(colors = colorRampPalette(c("white", "#00C1D4", "#FFED99","#FF7600"))(10))+
  labs(x=NULL,y=NULL)

p1 <- ggplot(exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_point(aes(size=`pct.exp`,color=`avg.exp.scaled`),
             shape=21,color="black",stroke=1)+
  theme(panel.background =element_blank(),
        axis.line=element_line(colour="black"),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=11,color="black"),
        panel.margin=unit(c(3,0,0,1), "cm"),
        plot.margin = margin(t = 5,r = 1,b = 1,l = 1,unit = 'cm'),
        axis.text.x = element_text(angle=45,
                                   vjust=1, 
                                   size=11,
                                   hjust=1,
                                   
                                   color = 'black')
  )+
  coord_cartesian(clip = 'off') +
  scale_color_gradientn(colors = colorRampPalette(c("#1E3163","#00C1D4","#FFED99","#FF7600"))(10))+
  labs(x=NULL,y=NULL)+ 
  geom_vline(xintercept=c(2.5,4.5,6.5,8.5,10.5,12.5), linetype="dotted",size=1)

library(jjAnno)
library(aplot)
p2 <- annoSegment(object = p1,
                  annoPos = 'top',
                  annoManual = F,
                  xPosition = c(1:45.5),
                  yPosition = c(15),
                  segWidth = 0.7,
                  #pCol = c(useMyCol('stallion',20),'orange'),
                  pCol=rep(c(rep("#1CC5FE",1),rep("#FB7D80",1)),7)
)
p3 <- annoSegment(object = p2,
                  annoPos = 'top',
                  annoManual = T,
                  xPosition = list(c(1,3,5,7,9,11,13),
                                   c(2,4,6,8,10,12,14)),
                  yPosition = 15.6,
                  segWidth = 0.7,
                  pCol = celltype_col[1:7],
                  addBranch = T,
                  branDirection = -1,
                  lwd = 2.5)
p4 <- annoRect(object = p3,
               annoPos = 'left',
               annoManual = T,
               yPosition = list(c(0.5,8.5),
                                c(8.5,14.5)),
               xPosition = c(-2.5,0.4),
               pCol = rep('white',3),
               pFill = useMyCol('calm',3.5),
               alpha = 0.5)

#fig4f
sigScores <- as.matrix(t(getSignatureScores(vis)))
data=as.data.frame(t(sigScores))
id=data.combined_fib@meta.data$name
type=factor(data.combined_fib@meta.data$fincell,levels=celltype_sub)
group=as.vector(data.combined_fib@meta.data$group)
for (i in gene_use) {
  Signature_score=data[i]
  Signature_score=Signature_score[,1]
  pdata_melt=data.frame(id,group,type,Signature_score)
    c <- ggplot(pdata_melt,
              aes(x=type, y=Signature_score, 
                  fill = group )) +
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = sample_color) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")+stat_compare_means(label = "p.signif")  
  ggsave(paste0("cell&group1_",i,".pdf"),plot = c, width=14, height=5)
  
  c<-ggplot(pdata_melt, aes(x = type, y = Signature_score, fill = group)) +
    geom_violin( aes(color=group),width=.8)+
    scale_fill_manual(values = sample_color)+
    scale_color_manual(values = sample_color)+
    theme_classic() +
    geom_boxplot( aes(color=group),fill="white",width=.8,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+stat_compare_means(label = "p.signif")
  ggsave(paste0("cell&group2_",i,".pdf"),plot = c, width=20, height=5)
   c <- ggplot(pdata_melt, aes(x = group, y = Signature_score, fill = group)) +
    geom_violin( aes(color=group))+
    scale_fill_manual(values = sample_color)+
    scale_color_manual(values = sample_color)+
    theme_classic() +
    geom_boxplot( aes(color=group),fill="white",width=.8,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+stat_compare_means()
  p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
  ggsave(paste0("group2_",i,".pdf"),plot = p, width=3.5, height=5)
  c <- ggplot(pdata_melt,
              aes(x=type, y=Signature_score, 
                  fill = type, 
              )) + 
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = celltype_col) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")  
  ggsave(paste0("cell1_",i,".pdf"),plot = c, width=8, height=5)
  
  c <- ggplot(pdata_melt, aes(x = type, y = Signature_score, fill = type)) +
    geom_violin( aes(color=type))+
    scale_fill_manual(values = celltype_col)+
    scale_color_manual(values = celltype_col)+
    theme_classic() +
    geom_boxplot( aes(color=type),fill="white",width=.8,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+stat_compare_means()
  ggsave(paste0("cell2_",i,".pdf"),plot = c, width=9, height=4)
    
  
}

#fig4g
#refer to fig3i

#fig4h
#refer to figS2g

#figS3a
#refer to fig4f

#figS3b
#refer to fig4f

#figS3c

male <- data.combined_1
femlale <- data.combined_2
male_DGEs <- list()
cells1 <- unique(male$fincell)
table(mlale$fincell,mlale$group)
for (i in 1:length(unique(male$fincell))) {
  data = subset(male, fincell==cells1[i])
  df <- FindMarkers(data,
                    group.by="group",
                    ident.1="B_ca",
                    ident.2="A_ca",
                    logfc.threshold = 0,
                    min.cells.group = 1,
                    min.pct = 0)
  male_DGEs[[i]] <- df
}

names(male_DGEs) <- cells1
femlale_DGEs <- list()
cells2 <- unique(femlale$fincell)
table(femlale$fincell,femlale$group)
for (i in 1:length(unique(femlale$fincell))) {
  data = subset(femlale, fincell==cells2[i])
  df <- FindMarkers(data,
                    group.by="group",
                    ident.1="B_lymph",
                    ident.2="A_lymph",
                      min.cells.group = 1,
                    logfc.threshold = 0,
                    min.pct = 0)
  femlale_DGEs[[i]] <- df
}

names(femlale_DGEs) <- cells2
for(i in seq_along(male_DGEs)){
  male_DGEs[[i]]$sex <- "M"
  male_DGEs[[i]]$cluster <- names(male_DGEs)[i]
  male_DGEs[[i]]$gene <- rownames(male_DGEs[[i]])
}
for(i in seq_along(femlale_DGEs)){
  femlale_DGEs[[i]]$sex <- "F"
  femlale_DGEs[[i]]$cluster <- names(femlale_DGEs)[i]
  femlale_DGEs[[i]]$gene <- rownames(femlale_DGEs[[i]])
}

myorder <- celltype_sub
male_DGEs <- male_DGEs[myorder]
femlale_DGEs <- femlale_DGEs[myorder]
for(i in seq_along(male_DGEs)){
  male_DGEs[[i]] <- male_DGEs[[i]][which(male_DGEs[[i]]$p_val < 1),]
  male_DGEs[[i]] <- male_DGEs[[i]][order(male_DGEs[[i]]$avg_log2FC, decreasing = TRUE), ]
}


for(i in seq_along(femlale_DGEs)){
    femlale_DGEs[[i]] <- femlale_DGEs[[i]][which(femlale_DGEs[[i]]$p_val < 1),]
  femlale_DGEs[[i]] <- femlale_DGEs[[i]][order(femlale_DGEs[[i]]$avg_log2FC, decreasing = TRUE), ]
  
}

top_male <- list()
for(i in seq_along(male_DGEs)){
  top_male[[i]] <- rbind(head(male_DGEs[[i]], 5), tail(male_DGEs[[i]],15))
}
top_female <- list()
for(i in seq_along(femlale_DGEs)){
  top_female[[i]] <- rbind(head(femlale_DGEs[[i]], 5), tail(femlale_DGEs[[i]],15))  
}

 male_DEGs_file <- do.call(rbind, male_DGEs)
 femlale_DGEs_file <- do.call(rbind, femlale_DGEs)
 write.csv(male_DEGs_file, file = "male_DEGs_file.csv")
 write.csv(femlale_DGEs_file, file = "femlale_DGEs_file.csv")
top_female <- do.call(rbind, top_female)
top_female <- top_female[,c(2,6:8)]
top_female$group <- paste0(top_female$sex,"_",top_female$cluster)
top_female <- top_female[,-c(2:3)]
top_female<-spread(top_female, group, avg_log2FC)
rownames(top_female) <- top_female$gene
top_female <- top_female[,-1]
library(ComplexHeatmap)
library(circlize)

cell=gsub("^..","",colnames(A))

annotation_col = data.frame(
  group = c(rep("B_ca.vs.A_ca",length(unique(male$fincell))),rep("B_lymph.vs.A_lymph",length(unique(male$fincell)))),
  celltype = c(cell))

row.names(annotation_col) <- colnames(A)
groupcolor <- c("#2f5688","#CC0000")

celltypecolor <- celltype_col

ann_colors <- list(group=groupcolor, celltype=celltypecolor) 
col_fun=colorRamp2(c(-2,0,2),c("#2f5688", "white","#CC0000"))

pdf("heat_two_compare_1.pdf",5,4.5)
ComplexHeatmap::pheatmap(as.matrix(A), 
                         scale = "none",
                         show_rownames = F,
                         show_colnames = F,
                         #col = colorRampPalette(c("#2f5688", "white","#CC0000"))(10),
                         col =col_fun,
                         column_title_gp  = gpar(fontsize = 10),
                         fontsize_row   =  gpar(fontsize = 8),
                         row_title = NULL,
                         cluster_rows = T，
                         cluster_cols = F,
                         column_split = annotation_col$group,
                         heatmap_legend_param = list(
                           title='Log2FC'),
                         annotation_col  = annotation_col,
                         annotation_colors = ann_colors,
                         annotation_legend = T,
                         annotation_names_col = F,
                         border_color = NA
                         )
dev.off()

#figS3d
#refer to fig3h

#figS3e
#refer to figS2h









