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

#figS5a
#refer to fig1c

#figS5b
#refer to fig3d

#figS5c
#refer to fig4a

#figS5d
#refer to fig3h,i

#figS5e
#refer to fig1c

#figS5f
vln.dat=FetchData(data.combined,c(gene,"fincell"))
vln.dat$fincell=factor(vln.dat$fincell,levels = celltype)
vln.dat=vln.dat[order(vln.dat$fincell),]
vln.dat$fincell=factor(vln.dat$fincell,levels = unique(vln.dat$fincell))
vln.dat.melt=vln.dat %>% 
  reshape2::melt(,heatmap_gene) %>%
  dplyr::rename("Gene"="variable") %>%
  group_by(fincell,Gene) %>%
  mutate(fillcolor=mean(value))
pal=colorRampPalette(RColorBrewer::brewer.pal(n = 9,name="YlOrRd"))(100)
vln.dat.melt[vln.dat.melt$fillcolor>2.5,]$fillcolor <-2.5
p1 = ggplot(vln.dat.melt,aes(x=fincell,y=value,fill=fillcolor))+
  geom_violin(linetype="blank",scale = "width")+
  scale_fill_gradientn(colors=pal,name="Average Expression")+
  facet_grid(Gene~.)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        #axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        axis.text.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "left")
p1

#figS5g
df=data.combined@meta.data
A <- prop.table(table(df$fincell, df$group), margin = 2)
A <- as.data.frame(A)
colnames(A) <- c("celltype", "group", "Freq")
ggplot(A,aes(x = group,y =Freq,
             group=celltype))+
  stat_summary(geom = 'line',fun='mean',cex=1,col='white')+#先要有折线
  geom_area(data = A,aes(fill=celltype))+#折线下面积，填充用celltype
  scale_fill_manual(values=celltype_col)+
  labs(x=NULL,y=NULL)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text = element_text(color = "black",size = 10))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('area_treat.pdf', width = 3.5, height = 4)

#figS5h
#refer to fig1c

#figS5i
#refer to figS3a

#figS5i
#refer to figS3a

#figS5j
fin_ratio=as.data.frame(fin_ratio)
fin_ratio$cell=rownames(fin_ratio)
fin_ratio$fin_ratio2=log2(fin_ratio$fin_ratio)
fin_ratio=fin_ratio[order(fin_ratio$fin_ratio2),]
fin_ratio$cell = factor(fin_ratio$cell, levels = fin_ratio$cell)
ggplot(fin_ratio,aes(x=cell,y=fin_ratio2))+geom_segment(aes(x=cell,xend=cell,y=0,yend=fin_ratio2),
                   size=1.5,color="#C9CACA",linetype="solid")+
  geom_point(size=5,shape=21,aes(color=cell,fill=cell))+
  theme_light()+
  ylim(-3, 4.5) +
  scale_color_manual(values = celltype_col) +
  scale_fill_manual(values = celltype_col) +
  theme(legend.position = "none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 0.8,  vjust = 0.6),
        axis.ticks.x=element_blank()
        #axis.ticks.y=element_blank()
        )+
  xlab("")+ylab("logFC of Ratio")


