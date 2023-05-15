args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(gridExtra)
  library(cowplot)
  library("ggsci")
})


#get protein-length accuracy values for QX1410 (automated and curated) and AF16
curated_QX <- read.table("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/plot/predhits-CURATION-VF", header = T, sep = "", dec = ".") 
merger_QX <- read.table("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/plot/predhits-QX1410-MERGER-ND", header = T, sep = "", dec = ".") 
AF16_WS280 <- read.table("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/plot/predhits-AF16-WS280-ND", header = T, sep = "", dec = ".") 

#subset based on shared matches between QX1410 and AF16
AFQX_matches <- AF16_WS280 %>% dplyr::filter(Subject_parent %in% curated_QX$Subject_parent)
QXC_subset <- curated_QX %>% dplyr::filter(Subject_parent %in% AFQX_matches$Subject_parent) %>% dplyr::select(delta)
QXA_subset <- merger_QX %>% dplyr::filter(Subject_parent %in% AFQX_matches$Subject_parent) %>% dplyr::select(delta)
AF_subset <- AFQX_matches %>% dplyr::select(delta)

#get counts for identical and 5% off matches
nrow(QXC_subset %>% dplyr::filter(delta==1))
nrow(QXC_subset %>% dplyr::filter(delta < 1.05 & delta >0.95 & !(delta ==1)))
nrow(AF_subset %>% dplyr::filter(delta==1))
nrow(AF_subset %>% dplyr::filter(delta < 1.05 & delta >0.95 & !(delta ==1)))
nrow(QXA_subset %>% dplyr::filter(delta==1))
nrow(QXA_subset %>% dplyr::filter(delta < 1.05 & delta >0.95 & !(delta ==1)))

#transoform data for histogram plot
colnames(QXC_subset) <- c('CURATION.VF')
colnames(QXA_subset) <- c('Merger')
colnames(AF_subset) <- c('AF16.WS280')
QXC_subset <- data.frame(as.table(as.matrix(QXC_subset))) %>% dplyr::mutate(Freq=log10(Freq))
QXC_subset <- QXC_subset[,2:3]
QXA_subset <- data.frame(as.table(as.matrix(QXA_subset))) %>% dplyr::mutate(Freq=log10(Freq))
QXA_subset <- QXA_subset[,2:3]
AF_subset <- data.frame(as.table(as.matrix(AF_subset))) %>% dplyr::mutate(Freq=log10(Freq))
AF_subset <- AF_subset[,2:3]
hs_QX <- hist(QXC_subset$Freq,breaks=seq(-1.5,1.5,by=0.01),plot = FALSE)
ds_QX <- data.frame(x = hs_QX$breaks,y = c(hs_QX$counts,NA))
hm_QX <- hist(QXA_subset$Freq,breaks=seq(-1.5,1.5,by=0.01),plot = FALSE)
dm_QX <- data.frame(x = hm_QX$breaks,y = c(hm_QX$counts,NA))
h280_AF16 <- hist(AF_subset$Freq,breaks=seq(-1.5,1.5,by=0.01),plot = FALSE)
d280_AF16 <- data.frame(x = h280_AF16$breaks,y = c(h280_AF16$counts,NA))
dbind <- bind_rows(list(`QX1410 (CURATED)`=ds_QX,`QX1410 (AUTOMATED)`=dm_QX,`AF16 (WS280)`=d280_AF16), .id="Query")
dbind$Query <- factor(dbind$Query, levels=c("AF16 (WS280)","QX1410 (AUTOMATED)","QX1410 (CURATED)")) 


#plot histogram (Figure 4)
bot <- ggplot() + 
  geom_step(data = dbind,aes(x = x,y = y,color=Query),stat = "identity") + coord_cartesian(ylim = c(0, 1500), xlim = c(-0.15, 0.15)) + 
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=14),
        axis.title.x=element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + scale_colour_manual(values=c("grey","blue","chartreuse2","darkgoldenrod1","red"))

top <- ggplot() + 
  geom_step(data = dbind,aes(x = x,y = y,color=Query),stat = "identity") + coord_cartesian(ylim = c(1550, 7000), xlim = c(-0.15, 0.15)) +
  theme(axis.line.x=element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background=element_rect(colour ="white", fill ="grey"),
        legend.margin = margin(6, 6, 6, 6)) + scale_colour_manual(values=c("grey","blue","chartreuse2","darkgoldenrod1","red"))

figureQX <- plot_grid(top,bot,rel_heights = c(1,3),align = 'v',ncol=1)
y.grob <- grid::textGrob("counts",
                         gp=gpar( col="black", fontsize=12), rot=90)
x.grob <- grid::textGrob(expression(paste("log"[10],"(Query protein length / N2 protein length)")),
                         gp=gpar(col="black", fontsize=12))


figure_final <- grid.arrange(arrangeGrob(figureQX,ncol=1, left = y.grob, bottom = x.grob))






