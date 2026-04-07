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

#fig6a
dotplot <- DotPlot(data.combined_fib, features =markers,split.by = "group",cols =col )+coord_flip()
data <- dotplot$data
colnames(data)
data <- separate(data = data, col = id, into = c("Celltype", "group"), sep = "_",remove = F)
data$Celltype=factor(data$Celltype,levels = unique(tmp$cell_type_name))
data=data[order(data$Celltype),]
data$id=factor(data$id,levels = unique(data$id))
p2 <-ggplot(data=data,aes(x=id,y= features.plot)) + 
  geom_tile(data=data,aes(x=id,y= features.plot,  fill=avg.exp.scaled),color="white",size=5)+
  scale_fill_gradient2(high="#CC3333", low="white")+
  #scale_color_gradient2(high='#CC3333', low="white")+
  theme_bw()+
  theme(axis.text =element_text(size = 10, color = "black"),
        axis.text.x=element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.margin = margin(-0.2,-0.2,0,0,'cm'),
        legend.key.height = unit(0.3,'cm'),
        legend.key.width = unit(0.5,'cm'),
  theme(plot.margin=unit(c(2.5, 2.5, 2.5, 2.5),'cm'))+ 
  guides(size=guide_legend(title = "Exp",order = 1))+
  coord_cartesian(clip = 'off')

p3<-annoSegment(object = p2,
            annoPos = 'top',
            aesGroup = T,
            aesGroName = 'Celltype',
            yPosition = 6.7,
            segWidth = 0.1,
            pCol=celltype_col,
            addText=T,
            textSize = 8,
            textCol = rep("black",length(celltype)),
            textRot = 90,
            vjust = 0.5,
            hjust = 0)

annoPoint2(object = p3,
           annoPos = 'top',
           xPosition = c(1:37),
           yPosition = 0.1,
           ptSize=0.9,
           pCol=rep(c("#1CC5FE","#FB7D80"),13)
)

#fig6b
 mat=data.combined_fib[["RNA"]]@data
 sigScores <- as.matrix(data.combined_fib[["RNA"]]@data[rownames(mat)%in%genes,])
 data=as.data.frame(t(sigScores))
 id=data.combined_fib@meta.data$name
 type=factor(data.combined_fib@meta.data$fincell,levels=celltype)
 group=as.vector(data.combined_fib@meta.data$group)
 
 cellratio1 <- meta%>%
   group_by(fincell, group) %>%
   summarise(avg.freq = mean(express),
             sd.freq = sd(express),
             se.freq = sd(express)/sqrt(n())) 
 cellratio1$fincell=factor(cellratio1$fincell,levels = celltype)
 
 p<- ggplot(data = cellratio1, mapping = aes(x = fincell, y = avg.freq,
                                         color =group, fill = group)) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
   scale_fill_manual(values=sample_color)+
   scale_color_manual(values=sample_color)+
   geom_bar(stat = "identity", position=position_dodge(),width = 0.8)+
   geom_errorbar(data = cellratio1,
                 aes(ymin=avg.freq-se.freq, ymax=avg.freq+se.freq),
                 width=0.3,
                 size = 0.5,
                 color = "black",
                 position=position_dodge(0.9))

#fig6c
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
                        data= dat) 
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit               = fit,
                conf.int          = F, 
                risk.table        = T, 
                risk.table.col    = "strata",
                palette           = mycol, 
                data              = dat,
                # xlim              = c(0,120), 
                size              = 1,
                break.time.by     = 5, 
                legend.title      = "",
                xlab              = "Time (years)",
                ylab              = "Overall survival",
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
p

#fig6d
PCA <-princomp(prop_ratio_fin)
PCA_fin=PCA$scores[,c(1,2)]
PCA_fin=as.data.frame(PCA_fin)
group=strsplit(as.vector(rownames(PCA_fin)),"_")
sam<-time<-origin<-c()
for (i in 1:length(group)) {
  group_sub=group[[i]]
  sam=c(sam,group_sub[1])
  time=c(time,group_sub[2])
  origin=c(origin,group_sub[3])
}
PCA_fin$sam=sam
PCA_fin$time=time
PCA_fin$origin=origin
PCA_fin$time=factor(PCA_fin$time,levels = rev(unique(PCA_fin$time)))
PCA_fin$time=factor(PCA_fin$time,levels = unique(PCA_fin$time))
ggplot(PCA_fin, aes(x = Comp.1, y = Comp.2,fill=sam,starshape = origin,color=sam,size=time)) +
  geom_star(alpha=0.7) +
  scale_starshape_manual(c(28,14,14,6,15))+
  scale_fill_manual(values = col)+
 scale_color_manual(values = col)+
  theme_classic()+
  scale_size_manual(values = c(3,6,9))

#fig6e
LR_sec=tile
receptor <- unique(LR_sec$receptor)
receptor
LR_sec$Target <- LR_sec$receptor#新建一列
for (i in 1:length(receptor)){
  LR_sec[, "receptor"][LR_sec[, "receptor"]==receptor[i]] = i
  
}
ligand <- unique(LR_sec$ligand)
LR_sec$ligand <- factor(LR_sec$ligand, levels = c(1:length(unique(LR_sec$ligand))))
LR_sec$receptor <- factor(LR_sec$receptor, levels = c(1:length(unique(LR_sec$receptor))))

if (length(unique(LR_sec$source))-length(unique(LR_sec$Target))>0) {
  p11 <- LR_sec %>% 
    ggparcoord(
      columns = 1:2,
      groupColumn = 3, 
      alphaLines = 0.6,scale="globalminmax",
      splineFactor=T)+
    scale_y_continuous(breaks = c(1:length(unique(LR_sec$source))),
                       labels =unique(LR_sec$source),
                       sec.axis = sec_axis(~./1, 
                                           breaks=c(1:length(unique(LR_sec$source))),
                                           labels=ylab_right))+
    theme(legend.position = 'none')+
    theme(plot.margin = margin(t = 0,  
                               r = 0,  
                               b = 0,  
                               l = 0,  
                               unit = "cm"))+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          legend.position = 'none',
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.left = element_text(colour =col))+
    scale_color_manual(values = col_unique)  
  
}else{
  p11 <- LR_sec %>% 
    ggparcoord(
      columns = 1:2, 
      groupColumn = 3, 
      showPoints = F, 
      alphaLines = 0.6,scale="globalminmax",
      splineFactor=T)+
    scale_y_continuous(breaks = c(1:length(unique(LR_sec$Target))),
                       labels =ylab_right,
                       sec.axis = sec_axis(~./1, 
                                           breaks=c(1:length(unique(LR_sec$Target))),
                                           labels=unique(LR_sec$Target)))+
    theme(legend.position = 'none')+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          legend.position = 'none',
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.left = element_text(colour =col))+
    scale_color_manual(values =col_unique)#修改下连线的分组颜色，dittoColors()是dittoSeq包的一个颜色函数
  
}

p11

p1 <- LR_sec %>% 
  ggparcoord(
    columns = 1:2, 
    groupColumn = 3, 
    showPoints = F,
    alphaLines = 0.6,
    scale="globalminmax"
    )+
  annotate('text', x = 2, 
           y = c(1:length(unique(LR_sec$Target))) , 
           label = unique(LR_sec$Target),size=2 ,family="B",hjust=0)+
  annotate('text', x = 1, 
           y = c(1:length(unique(LR_sec$source))) , 
           label = unique(LR_sec$source),size=2 ,family="B",hjust=1,col=col)+
  theme(legend.position = 'none')+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  scale_color_manual(values = col_unique )
library(Seurat)
data.combined=data.combined
Idents(data.combined)<-"group"
HC <-data.combined
ligand_dot <- DotPlot(HC, features = ligand,assay = "RNA")
ligand_exp <- ligand_dot$data
p2 <- ggplot(ligand_exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  theme_bw()+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "right")+
  theme(axis.text.x=element_text(angle=90,hjust = 0,vjust=1,colour = 'black',size = 9),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.3,'cm'),
        legend.key.width = unit(0.3,'cm'),
        legend.title = element_text(size=8,vjust = 1),
        legend.text = element_text(size = 5))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#1A5592','white',"#B83D3D"),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"),
                        name = "Scaled expr")+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(title = "",order = 1,nrow = 2))+
  theme(plot.margin=unit(c(1, 1, 1, 1),'cm'))
data.combined=data.combined_fib
Idents(data.combined)<-"group"
immune=data.combined

Target=LR_sec$Target
Target=strsplit(Target,"_")
target_fin=c()
for (i in 1:length(Target)) {
  tmp=Target[[i]]
  target_fin=c(target_fin,tmp[1])
}
p5 <- ggplot(receptor_exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  theme_bw()+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "left")+
  theme(axis.text.x=element_text(angle=90,hjust = 0,vjust=1,colour = 'black',size = 9),
        axis.text.y = element_blank(),
        title = element_text(colour = 'black',size = 8),
        plot.title = element_text(hjust=0.5),
        legend.key.height = unit(0.3,'cm'),
        legend.key.width = unit(0.3,'cm'),
        legend.title = element_text(size=8,vjust = 1),
        legend.text = element_text(size = 6))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#1A5592','white',"#B83D3D"),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"),
                        name = "Scaled expr")+
  labs(x=NULL,y=NULL)+
  ggtitle("")+
  theme(plot.margin=unit(c(1, 1, 1, 1),'cm'))

p7 <- p5
p7

#fig6g,h

library(ggvoronoi)
library(ggplot2)
library(dplyr)
library(ggforce)
celltype=sort(unique(data_fin$celltype))
names(col)= celltype
celltype_col=col
sam=gsub(".txt","",sam)

for (i in sam) {
data_fin_sub=data_fin[data_fin$study==i,]
data_fin_sub=data_fin_sub[order(data_fin_sub$celltype),]
point<-ggplot(data=data_fin_sub,aes(x=Centroid_X,y=Centroid_y))+ geom_point(
  aes(colour=celltype, 
      fill=celltype), size = 3,alpha=1)+
  scale_color_manual(values=celltype_col) +
  scale_fill_manual(values=celltype_col)+
  theme()+theme(panel.background = element_rect(fill = 'white'),
                        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave(paste0("point1_",i,".pdf"),plot =point, width = 60, height = 60,limitsize = F)

point<-ggplot(data=data_fin_sub,aes(x=Centroid_X,y=Centroid_y))+ geom_point(
  aes(colour=celltype, 
      fill=celltype), size = 3,alpha=1)+
  scale_color_manual(values=celltype_col) +
  scale_fill_manual(values=celltype_col)+
  theme()+theme(panel.background = element_rect(fill = 'black'),
                panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave(paste0("point2_",i,".pdf"),plot =point, width = 60, height = 60,limitsize = F)
data_fin3=data_fin_sub

map<-  ggplot(data_fin3, aes(Centroid_X, Centroid_y, group = -1L)) +
  geom_voronoi_tile(aes(fill = celltype,colour = celltype), normalize = F,expand = unit(-1.5, 'pt'), radius = unit(0, 'pt'), max.radius = 10)+
scale_color_manual(values=celltype_col) +
  scale_fill_manual(values=celltype_col)+
  theme()+
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
ggsave(paste0("voronoi1_",i,".pdf"),plot =map, width = 250, height = 250,limitsize = F)


map<-  ggplot(data_fin3, aes(Centroid_X, Centroid_y, group = -1L)) +
  geom_voronoi_tile(aes(fill = celltype,colour = celltype), normalize = F,expand = unit(-1.5, 'pt'), radius = unit(0, 'pt'), max.radius = 10)+
  scale_color_manual(values=celltype_col) +
  scale_fill_manual(values=celltype_col)+
  theme()+
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
ggsave(paste0("voronoi2_",i,".pdf"),plot =map, width = 250, height = 250,limitsize = F)


}



p1<-ggplot(data = count2, aes(x = sam, y = Freq, fill = cell)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=celltype_col) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  rotate_x_text()+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))
p1


data_for_ratio1=data_fin[,c("study","celltype")]
count=as.data.frame(table(data_for_ratio1$study))
count = count %>% 
  separate(Var1, into = c('sam','group'), sep = '_')
count$sam=factor(count$sam,levels = sam_sort)

p2<-ggplot(data=count)+
  geom_bar(aes(x=sam, y=Freq, fill=group), stat='identity')+
  scale_fill_manual(values=group_cols) +theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="count")+
  theme(axis.text.x = element_blank())+
theme(axis.ticks.x = element_line(color =  NA))
p2 
library(patchwork)
p <- p2 + p1 + plot_layout(ncol = 1, heights = c(1.2, 2))

p <- ggplot(count_sub, aes(x = group, y = Freq, color = group)) +
    geom_jitter(size = 7, width = 0.15, show.legend = FALSE,alpha=0.6) +  
    scale_color_manual(values = group_cols) +  
    theme(panel.grid = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(color = 'black')) + 
    labs(x = 'Group', y = 'Cell Freq') + 
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, 
                 geom = 'crossbar', width = 0.3, size = 0.25) +  
    stat_summary(fun.data = function(x) median_hilow(x, 0.5), 
                 geom = 'errorbar', width = 0.25, size = 0.4)+ 
    stat_compare_means(label = "p.signif",comparisons = my_comparisons,method = "t.test")
  ggsave(paste0(i,"_ratio.pdf"), plot = p, width = 3.2, height = 3.5)

#fig6j,i,figS4e
library(ggridges)
library(RColorBrewer)
library(ggpubr)

G01_all=data.frame()
for (i in unique(data_fin$study)) {
  G01_sub=load(paste0("Gcross",i,"_",cell1,"_",cell2,".RData"))
  G01_sub=data.frame(G01$r,G01$rs)
  G01_sub$sam=i
  G01_all=rbind(G01_all,G01_sub)
}
library(tidyverse)
colnames(G01_all)<-c("r","value","sam_group")
mydata1 = G01_all %>% 
  separate(sam_group, into = c('sam','group'), sep = '_')
sample_color <- c("#1CC5FE","#A48AD3","#FB7D80")
ggplot(mydata1, aes(x=r,y=value,color=group,fill=group)) +
  labs(y = "Border corrected estimate of G")+
  labs(x = "Distance")+
  scale_fill_manual(values = sample_color) +
  scale_color_manual(values = sample_color) +
  geom_smooth(span = 1)+theme_classic()


mydata1=data.frame(G01$r,G01$theo,G01$rs)# border corrected estimate of "G"
colnames(mydata1)<-c("r","theo","rs")
library(ggridges)
library(RColorBrewer)
mydata1=mydata1[1:100,]
p<-ggplot(mydata1, aes(x=r)) +
  geom_ridgeline_gradient( aes(y=ymin, height = ymax-ymin,  fill = rs-theo)) +
  geom_line(aes(y=theo),color="blue",size=1)+#color="black",
  geom_line(aes(y=rs),color="red",size=1)+#color="black",
  
  scale_fill_gradientn(colours= colorRampPalette(c("#1E3163","#00C1D4","#FFED99","#FF7600"))(10),
                       name = "Value", 
                       limits=c(-max(mydata1$ymax-mydata1$ymin),max(mydata1$ymax-mydata1$ymin)), 
                       breaks = c(-max(mydata1$ymax-mydata1$ymin),0,max(mydata1$ymax-mydata1$ymin)))+
  theme_classic2()+
  theme(#legend.position = c(0.15,0.8),
    legend.background = element_blank()) 
p

library(ggpubr)
library(gghalves)
library(tidyverse)
plot_dat <- read.csv("./distance2_100.csv", row.names=1)
plot_dat$sam_group=rownames(plot_dat)
plot_dat = plot_dat %>% 
  separate(sam_group, into = c('sam','group'), sep = '_')
plot_dat$value=plot_dat$treg_to_tfam_mean
ggplot(plot_dat,aes(x = group,y=value,fill=group,
                    color=group))+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  geom_point(aes(x = group,y=value,
                 color = group),
             position = position_jitter(width =0.03),size =2, shape = 20)+
  geom_boxplot(outlier.shape = NA, width =0.12,alpha=0.55)+
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=2,trim=F,color=NA,alpha=0.55)+
  theme_bw()+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(color =cols,size=12 ))+
  stat_compare_means(label = "p.signif",comparisons = my_comparisons,method = "t.test")