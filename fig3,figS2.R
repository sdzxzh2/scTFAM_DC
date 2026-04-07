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


#fig3a
#refer to fig1c

#fig3b
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
infercnv_obj = readRDS("./run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$cCON
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- as.vector(unlist(test_loc))
anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)
gn <- rownames(expr)
geneFile <- read.table("./all_sam_geneFile.txt", header=FALSE, row.names=NULL)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]
expr=expr[,colnames(expr)[test_loc]]
geneFile_sub=geneFile[rownames(expr),]
kmeans_df_s=data.frame(colnames(expr))
rownames(kmeans_df_s)=kmeans_df_s$colnames.expr.
meta=data.combined@meta.data
kmeans_df_s$group=meta[rownames(kmeans_df_s),"group"]plot_lis=list()
for (i in c(1:length(expr_lis))) {
  dat1=geneFile_sub[,2:4]
  dat1$cnv=rowMeans(expr_lis[[i]])
  dat1_up=dat1[dat1$cnv>1,]
  dat1_no=dat1[dat1$cnv==1,]
  dat1_down=dat1[dat1$cnv<1,]
  dat1_down$cnv=-dat1_down$cnv
  dat1_plot=list(dat1_up,dat1_no,dat1_down)
  plot_lis[[i]]=dat1_plot
  
}
pdf(file = "cric.pdf",
    width = 7.5,height = 7.5)
circos.initializeWithIdeogram(chromosome.index = paste0("chr",c(1:22)))
circos.track(ylim = c(0, 1),bg.col = 'white', 
             bg.border = 'white', 
             track.height = 0.05)
circos.genomicTrackPlotRegion(plot_lis[[1]], 
                              ylim = c(-1,1.5),
                              
                              track.height = 0.1,
                              
                              
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "h",baseline = 0, border =  if (value>1) {"#F36F75"
                                  
                                }else if(value<1){"#73A2F5"}else{"#1CC5FE"},lwd = 0.1,area = F,col=if (value>1) {"#F36F75"
                                  
                                }else{"#73A2F5"})},
                              
                              bg.col = '#1CC5FE', 
                              bg.border = '#1CC5FE'
)
circos.genomicTrackPlotRegion(plot_lis[[2]], 
                              ylim = c(-1,1.5),
                              
                              track.height = 0.1,
                              
                              
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "h",baseline = 0, border =  if (value>1) {"#F36F75"
                                  
                                }else if(value<1){"#73A2F5"}else{"#6FC7CF"},lwd = 0.1,area = F,col=if (value>1) {"#F36F75"
                                  
                                }else{"#73A2F5"})},
                              
                              bg.col = '#6FC7CF', 
                              bg.border = '#6FC7CF'
)
circos.genomicTrackPlotRegion(plot_lis[[3]], 
                              ylim = c(-1,1.5),
                              
                              track.height = 0.1,
                              
                              
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "h",baseline = 0, border =  if (value>1) {"#F36F75"
                                  
                                }else if(value<1){"#73A2F5"}else{"#FBA27D"},lwd = 0.1,area = F,col=if (value>1) {"#F36F75"
                                  
                                }else{"#73A2F5"})},
                              
                              bg.col = '#FBA27D', 
                              bg.border = '#FBA27D'
)

circos.genomicTrackPlotRegion(plot_lis[[4]], 
                              ylim = c(-1,1.5),
                              
                              track.height = 0.1,
                              
                              
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "h",baseline = 0, border =  if (value>1) {"#F36F75"
                                  
                                }else if(value<1){"#73A2F5"}else{"#FB7D80"},lwd = 0.1,area = F,col=if (value>1) {"#F36F75"
                                  
                                }else{"#73A2F5"})},
                              
                              bg.col = '#FB7D80', 
                              bg.border = '#FB7D80'
)

circos.clear()
dev.off()

#fig3d
library(ggplot2)
library(dplyr)
library(reshape2)
bar_plot_input$N=-bar_plot_input$N
x=melt(bar_plot_input,id.vars = "cell")
x$ylab=x$cell
x$ylab=factor(x$ylab,levels = bar_plot_input[order(bar_plot_input$T),]$cell)
x$variable=factor(x$variable,levels = c("T","N"))
x$variable=factor(x$variable,levels = c("T","N"))
col <- c("#FB7D80","#1CC5FE")
ggplot(x,aes(x=ylab,y=value,fill=variable))+
  geom_col(position=position_dodge(0),width = 1)+
  coord_flip()+
  ylim(-1.65,1.65)+
  labs(x=NULL,y=NULL)+
  theme_bw()+
  scale_fill_manual(guide = guide_legend(title = NULL),values = col)+
  theme(legend.position ="top", 
       # axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(), 
        legend.key.size=unit(0.8,"line"),
        axis.text.y = element_text(color =celltype_col[levels(x$ylab)],size = 12))
  
#fig3e,f,figS2d,e,f
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(GeneNMF) 
library(Seurat)
library(msigdbr)
library(fgsea)
library(viridis)
seu.list <- SplitObject(data.combined_fib, split.by = "orig.ident")
geneNMF.programs <- multiNMF(seu.list, 
                             assay="RNA", slot="data", 
                             k=4:8, L1=c(0,0),
                             do_centering=TRUE, 
                             nfeatures = 2000)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nprograms=9,
                                        max.genes=50,
                                        hclust.method="ward.D2",
                                        min.confidence=0.1)

genes=geneNMF.metaprograms$metaprograms.genes
writeGmtPathways(genes,"genes.gmt")
ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0.1,0.8),
                       palette = viridis(100, option = "A", direction = -1),
                       show_rownames = F)

geneNMF.metaprograms$metaprograms.metrics
t(as.data.frame(lapply(geneNMF.metaprograms$metaprograms.genes, head)))
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(data.combined_fib), category = "H",species = "Homo sapiens", pval.thr = 0.2)
  })

head(top_p)

top_p_path=lapply(top_p,FUN = function(x){x[[1]]})
writeGmtPathways(top_p_path,"top_hall.gmt")
mp.genes <- geneNMF.metaprograms$metaprograms.genes
data.combined_fib <- AddModuleScore_UCell(data.combined_fib, features = mp.genes, ncores=14, name = "")
VlnPlot(data.combined_fib, features=names(mp.genes), group.by = "group",
        pt.size = 0, ncol=5)
library(viridis)
FeaturePlot(data.combined_fib, features = names(mp.genes), ncol=4) &
  scale_color_viridis(option="B") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())
celltype_sub=celltype
sigScores <- data.combined_fib@meta.data[,names(mp.genes)]
data=sigScores
id=data.combined_fib@meta.data$name
type=factor(data.combined_fib@meta.data$fincell,levels=celltype_sub)
group=as.vector(data.combined_fib@meta.data$group)
gene_use=colnames(data)
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
}


library(GSVA)
gmt_list <- split(gmt$gene, gmt$term) 
gsva_mat <- gsva(expr=as.matrix(expr), 
                 gset.idx.list=        , 
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores()-2)
risk=dat$sig
res.cut=surv_cutpoint(dat,time="OS.time",
                      event ="OS",variables="sig")
res.cut=res.cut$cutpoint$cutpoint

risk<-as.vector(ifelse(risk >res.cut,"high","low"))
dat$group<-risk
sur=dat
fitd <- survdiff(Surv(OS.time,OS) ~ group,
                 data      = dat,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = dat,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

ps <- pairwise_survdiff(Surv(OS.time, OS)~ group,
                        data            = dat,
                        p.adjust.method = "none") 
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit= fit,
                conf.int= F, 
                risk.table= T, 
                risk.table.col    = "strata",
                palette= mycol, 
                data= dat,
                  size= 1,
                break.time.by= 5, 
                legend.title= "",
                xlab= "Time (years)",
                ylab= "Overall survival",
                risk.table.y.text = FALSE,
                tables.height     = 0.3) 
p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.001, " < 0.001", 
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 0, y = 0.55, 
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)

#fig3h

library(DESeq2)
library(Seurat)
library(IHW)
library(monocle)
library(tidyverse)
library(magrittr)
library(Seurat) 
library(pheatmap)
library(RColorBrewer)
library(aplot)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

data.combined=data.combined_fib
RA_matrix<-as(as.matrix(data.combined@assays$RNA@counts), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(RA_matrix),gene_short_name=rownames(RA_matrix))
rownames(feature_ann)<-rownames(RA_matrix)
RA_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<- data.combined@meta.data
rownames(sample_ann)<-colnames(RA_matrix)
RA_pd<-new("AnnotatedDataFrame", data =sample_ann)
monocle_cds<-newCellDataSet(RA_matrix,phenoData =RA_pd,featureData =RA_fd,expressionFamily=negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- setOrderingFilter(monocle_cds, gene)


p1 <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 1.5)+
  scale_color_gradientn(colors  = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))
p2 <- plot_cell_trajectory(monocle_cds, color_by = "State", cell_size = 2)

expr <- infercnv_obj@expr.data
expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)
#kmeans_df_s$CB=rownames(kmeans_df_s)
meta=data.combined@meta.data
meta$CB=rownames(meta)
CNV_score=CNV_score%>%inner_join(meta,by="CB")
meta=CNV_score
df <- data.combined@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = data.combined@meta.data$fincell)
colnames(df)
df2=df
df2$CB=rownames(df2)
meta=meta%>%inner_join(df2,by="CB")
score=quantile(meta$CNV_score)
cnv_type=c()
for (i in 1:nrow(meta)) {
  tmp=meta$CNV_score[i]
  if (tmp>score[4]|tmp==score[4]) {
    cnv_type[i]="CNV_D"
  }
  if (tmp>score[3]&tmp<score[4]|tmp==score[3]) {
    cnv_type[i]="CNV_C"
  }
  if (tmp>score[2]&tmp<score[3]|tmp==score[2]) {
    cnv_type[i]="CNV_B"
  }
  if (tmp>score[1]&tmp<score[2]|tmp==score[1]) {
    cnv_type[i]="CNV_A"
  }
  
  
}
meta$group=cnv_type
p<-plot_cell_trajectory(monocle_cds, color_by = "State", cell_size = 1.5)+
  scale_colour_manual(values = col)
meta_monocle2=p$data
meta_monocle2$CB=meta_monocle2$sample_name
meta=meta%>%inner_join(meta_monocle2,by="CB")
meta=meta[order(meta$group2.x),]
ggplot(meta,aes(data_dim_1,data_dim_2,color=group.x))+ geom_point(size=0.6,alpha=0.7)+theme_classic()+
  scale_color_manual(values = c("#71BC7A","#7DC6C5","#FFFA94","#FF935F"))
ggsave("monocle2_CNV_state2.pdf", width=4, height=3)


#fig3i
library(VISION)
sigScores <- as.matrix(t(getSignatureScores(vis)))
path=rownames(as.matrix(t(getSignatureScores(vis))))

for (i in path) {
  ene=sigScores[i,]
  plot.dat$gene=gene
  p1<-ggplot(data=plot.dat,aes(x=Pseudotime,y=gene))+
    geom_smooth(aes(color=group), 
                span = 1, size = 0.7,
                method = "loess",
                formula = 'y~x')+
    scale_color_manual(values=sample_color)+
    theme_bw(base_size = 15)+
    theme(panel.grid = element_blank())+
    labs(y="Expression")
  
  ggsave(plot=p1,filename=paste0(i,".pdf"),width = 4.8, height = 3.5)
  
 
}

#fig3j,figS2i
library(CellChat)
library(ktplots)
library(SingleCellExperiment)
library(reshape2)
library(Seurat)
library(circlize)
library(igraph)
library(ggraph)
library(ggrepel)
library(RColorBrewer)
library(grid)
library(ggplot2)

names(col)=unique(plot_source$target)
vertices <- data.frame(
  name = c(plot.data$from, plot.data$to),
  celltype = c(plot.data$source, plot.data$target),
  type = rep(c("ligand", "receptor"), each = nrow(plot.data)),
  frac = c(plot.data$L_Fraction, plot.data$R_Fraction),
  label = c(plot.data$ligand, plot.data$receptor),
  logpadj = c(plot.data$logpadj, plot.data$logpadj),
  logFC = c(plot.data$logFC, plot.data$logFC)
)
vertices <- vertices[!duplicated(vertices$name), ]
vertices$cellident <- paste0(vertices$celltype, "_", vertices$type)
vertices.head <- data.frame(matrix(ncol = ncol(vertices), 
                                   nrow = length(unique(c(vertices$celltype, vertices$cellident)))+1))
colnames(vertices.head) <- colnames(vertices)
vertices.head$name <- c("root", unique(vertices$celltype), unique(vertices$cellident))
vertices <- rbind(vertices.head[, colnames(vertices)], vertices)

edgelist <- rbind(
  data.frame("from" = paste0(plot.data$source, "_ligand"), "to" = plot.data$from),
  data.frame("from" = paste0(plot.data$target, "_receptor"), "to" = plot.data$to)
)
edgelist.head <- data.frame("to" = unique(na.omit(vertices$cellident)))
edgelist.head$from <- gsub("_receptor|_ligand", "", edgelist.head$to)
edgelist.head <- rbind(edgelist.head, data.frame("from" = "root", "to" = unique(na.omit(vertices$celltype))))

mygraph <- graph_from_data_frame(d = edgelist, directed = T, vertices=vertices)
d <- igraph::as_data_frame(mygraph, "both")
from  <-  match(plot.data$from, vertices$name)
to  <-  match(plot.data$to, vertices$name)
e.fun <- get_con(from = from, to = to)
e.df <- e.fun(create_layout(mygraph, layout = "dendrogram"))
e.df$logpadj[is.na(e.df$logpadj)] <- 0

g <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to),
                   aes(colour = logFC), width = 2, tension = 1,
                   arrow = arrow(length = unit((e.df$logpadj)*0.2, "mm"))) +
   geom_text_repel(aes(x = x, y = y, label = label), 
                  segment.square = TRUE, segment.inflect = TRUE, 
                  segment.size = 0.2, force = 0.5, size = 3.5, 
                  force_pull = 0) +
  scale_edge_colour_gradientn(colours =colorRampPalette(c("#1E3163","#00C1D4","#FFED99","#FF7600"))(10), limits = c(0, 1.8)) +   
  scale_shape_manual(values = c("ligand" = 1, "receptor" = 19)) +   
  scale_color_manual(values = celltype_col) +                          
  theme_void()

g

#figS2a
kmeans.result <- kmeans(t(cnv), 7)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) 
head(kmeans_df_s)
meta=data.combined@meta.data
group=meta[rownames(kmeans_df_s),"group"]
sam=meta[rownames(kmeans_df_s),"stim"]
kmeans_df_s$group=group
kmeans_df_s$sam=sam
sam_col=all_col[1:length(unique(kmeans_df_s$sam))]
names(sam_col)=sort(unique(kmeans_df_s$sam))
kmeans_df_s2=kmeans_df_s[,c("kmeans_class","group","sam")]
kmeans_df_s2=kmeans_df_s2[order(kmeans_df_s2$kmeans_class),]
kmeans_df_s2=kmeans_df_s2[order(kmeans_df_s2$group),]
kmeans_df_s2=kmeans_df_s2[,c(2:3)]

pdf("p.pdf",width = 15,height = 8)
ht = Heatmap(t(expr)[rownames(kmeans_df_s2),], 
             col = colorRamp2(c(0.5,1,1.5), c("#377EB8","#F0F0F0","#E41A1C")), 
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
             column_gap = unit(2, "mm"),  
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.5,1,1.5),legend_height = unit(3, "cm")),
             top_annotation = top_anno,left_annotation = left_anno, 
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

#figS2b

library(tidyverse)
r="r30"
sam=unique(data.combined$group)
alltype=c()
for (t in sam) {
  
  group_cnvtype_raw=read.table(paste0("./",t,"/",r,"/CNVgroup_and_CNVtype_in_sampleA.txt"),header = T,sep = "\t",stringsAsFactors = F)
  alltype=c(alltype,group_cnvtype_raw$cnv_type)
}

some.sample.stat=data.frame(rep("tmp",length(alltype)),1)

for (tt in sam) {
  
  dir=paste0("./",tt,".cell_groupings")
  cell_group=read.csv(dir,header = T,sep = "\t",stringsAsFactors = F)
  cell_group=cell_group[!str_detect(cell_group$cell_group_name,"references"),]
  cell_group$cell_group_name=str_replace(cell_group$cell_group_name,"all.*observations\\.","")
  
  group_cellcount=as.data.frame(table(cell_group$cell_group_name))
  colnames(group_cellcount)=c("cell_group_name","cellcount")
  group_cellcount$cellratio=group_cellcount$cellcount / sum(group_cellcount$cellcount)
  
  group_cnvtype=read.table(paste0("./",tt,"/",r,"/CNVgroup_and_CNVtype_in_sampleA.txt"),header = T,sep = "\t",stringsAsFactors = F)
  group_cnvtype=group_cnvtype%>%inner_join(group_cellcount,by="cell_group_name")
  
  cellpercent=c()
  for (i in alltype) {
    if(i %in% unique(group_cnvtype$cnv_type)){
      tmp=sum(group_cnvtype[group_cnvtype$cnv_type == i,"cellratio"])
      cellpercent=append(cellpercent,tmp)
    }else{
      cellpercent=append(cellpercent,0)
    }
  }
  names(cellpercent)=alltype
  one.sample.stat=as.data.frame(cellpercent)

  some.sample.stat[[tt]]=one.sample.stat
  
  
}

some.sample.stat=as.data.frame(some.sample.stat)
rownames(some.sample.stat)=alltype
colnames(some.sample.stat)=sam
library(RColorBrewer)
library(scales)
library(pheatmap)
pdf(paste0("heat_cnv1.pdf"),3.8,12)
pheatmap(some.sample.stat,
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(brewer.pal(9,"PuRd"))(100),
         show_rownames=T,
         
         #rev(brewer.pal(n = 7, name ="RdYlBu"))
         #brewer.pal(9,"PuRd")
         #brewer.pal(9,"RdPu"),
         border_color = "grey"
         
)
dev.off()

#figS2g
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(SCENIC)

scenicOptions <- initializeScenic(
  org="mgi", 
  dbDir="./cisTarget_databases", # RcisTarget databases location
  dbs="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",# file name of motif database
  datasetTitle="SCENIC on Mouse Cell Atlas", # choose a name for your analysis
  nCores=1
)

genesKept <- geneFiltering(
  exprMat=dge, 
  scenicOptions=scenicOptions,
  minCountsPerGene = 1,
  minSamples = 20
)
dge <- dge[genesKept, ]
dim(dge)
if(!dir.exists("output")) {
  dir.create("output")
}
saveRDS(cell.info, "output/s1_cell.info.rds")
source("utils/add_cellAnnotation.R")
saveLoom <- function(exprMat, output){
  cellInfo2 <- data.frame(
    row.names = rownames(cell.info),
    cellType =cell.info$fincell
  )
  loom <- build_loom(output, dgem=exprMat)
  loom <- add_cellAnnotation(loom, cellInfo2)
  close_loom(loom)
}
saveLoom(dge, "output/s1_avg20_rep1.loom")
loom <- build_loom("output/s1_exprMat.loom", dgem=dge)
loom <- add_cellAnnotation(loom, cell.info)
close_loom(loom)



library(ggtern)
library(scales)
library(ggplot2)
library(reshape2)
library("Seurat")
library("pheatmap")
library("viridis")
avgData <- sigScores %>% 
  apply(1, function(x){
    tapply(x, meta, mean) # ExpMean
  }) %>% t

avgData=avgData[,group_choose]
ggtern(data=avgData, aes(x=A,
                           y=D,
                           z=B_C,color=tf))+geom_point(aes(size=0.5), alpha=0.5)+# scale_color_manual(values =sample_color)+
   scale_color_gradientn(colors = colorRampPalette(c("#1E3163","#00C1D4","#FFED99","#FF7600"))(10))+
  stat_density_tern(n = 100,h=1, alpha=0.7,color="black")+
  theme(tern.panel.background = element_rect(fill = "white"), 
        tern.panel.grid.minor = element_line(color = "gray90"), 
        tern.axis.arrow.show = TRUE, 
        tern.axis.arrow.text.L = element_text(color = 'black'),  
        tern.axis.arrow.text.T = element_text(color = 'black'),
        tern.axis.arrow.text.R = element_text(color = 'black'),
        tern.axis.arrow.sep = 0.1,
        tern.panel.grid.major.T = element_line(color = 'gray92', linetype = 1, linewidth = 0.8), 
        tern.panel.grid.major.L = element_line(color = 'gray92', linetype = 1, linewidth = 0.8),
        tern.panel.grid.major.R = element_line(color = 'gray92', linetype = 1, linewidth = 0.8),
        tern.panel.grid.minor.T = element_line(color = 'gray94', linetype = 1, linewidth = 0.8), 
        tern.panel.grid.minor.L = element_line(color = 'gray94', linetype = 1, linewidth = 0.8),
        tern.panel.grid.minor.R = element_line(color = 'gray94', linetype = 1, linewidth = 0.8),
        tern.axis.title.L = element_text(color = '#1CC5FE', size = 11),
        tern.axis.title.T = element_text(color = '#FB7D80', size = 11),
        tern.axis.title.R = element_text(color = '#FBA27D', size = 11),
        tern.axis.text.L = element_text(size = 10,face = 'bold'),
        tern.axis.text.R = element_text(size = 10,face = 'bold'),
        tern.axis.text.T = element_text(size = 10,face = 'bold'),
        tern.axis.vshift = 0.04,
        tern.axis.line.T = element_line(linewidth = 0.8),
        tern.axis.line.R = element_line(linewidth = 0.8),
        tern.axis.line.L = element_line(linewidth = 0.8))
ggsave('p.pdf',width = 4,height = 4)


#figS2h
library(data.table)
library(pbapply)
library(plyr)
library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)
library(RColorBrewer)
library(dplyr)

rasMat <- readRDS("./output/s5_avg20_rep1.rasMat.rds")
s3_avg20_rep1.regulons <- read.delim("./output/s3_avg20_rep1.regulons.txt", header=FALSE, row.names=1)
rasMat2=t(rasMat)
rasMat_sub=rasMat2[,colnames(rasMat2)%in%rownames(data.combined@meta.data)]
row.names(rasMat_sub)=row.names(s3_avg20_rep1.regulons) 
fib.cell=colnames(rasMat_sub)
sub_CT<- subset(data.combined, name %in% fib.cell)
data.combined=sub_CT
data.combined$fincell=Idents(data.combined)
sigScores <-as.matrix(rasMat_sub)
Idents(data.combined) <- "group"
data.combined@assays$RNA@counts=sigScores
data.combined@assays$RNA@data=sigScores
data.combined@assays$RNA@scale.data=sigScores
data.combined.markers <- FindAllMarkers(data.combined, 
                                        only.pos = TRUE, 
                                        min.pct = 0.25,
                                        assay ="RNA" ,
                                        logfc.threshold = 0)

saveRDS(data.combined.markers, "degs.rds")
top10 <- data.combined.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
RotatedAxis2 <- function (...) 
{
  rotated.theme <- theme(axis.text.x = element_text(angle = 90, 
                                                    hjust = 1), validate = TRUE, ...)
  return(rotated.theme)
}
pdf("p.pdf",8,3.5)
DotPlot(data.combined, 
        assay ="RNA",
        features = unique(top10$gene))+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  RotatedAxis2()+scale_x_discrete("")+scale_y_discrete("")
dev.off()


