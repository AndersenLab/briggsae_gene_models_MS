library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ape)
library(genemodel)
library(cowplot)
library(gridExtra)
library(gridtext)
library(showtext)
showtext_auto()


##This script parses GFF features of AF16 and QX1410 attf-3 and pat-4 orthologs (Supplemental Fig 1)

fourgene <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/busco_dup_figure/four_new.gff") %>%
  dplyr::filter(type=="gene") %>%
  dplyr::select(start,end)

four <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/busco_dup_figure/four_new.gff") %>%
  dplyr::select(type,start,end,attributes) %>%
  dplyr::mutate(tair=ifelse(type=="CDS","coding_region",
                            ifelse(type=="five_prime_UTR","5' utr",
                                   ifelse(type=="three_prime_UTR","3' utr",as.character(type))))) %>%
  dplyr::filter(!(tair=="mRNA" | tair =="gene" | tair =="stop_codon" | tair =="start_codon")) %>%
  tidyr::separate(attributes,into = c("ID","Parent"), sep = ";") %>%
  dplyr::select(-ID,-type) %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(intron=ifelse(any(tair=="intron"),"wi","ni")) %>%
  dplyr::ungroup()

two <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/busco_dup_figure/two_new.gff") %>%
  dplyr::select(type,start,end,attributes) %>%
  dplyr::mutate(tair=ifelse(type=="CDS","coding_region",
                            ifelse(type=="five_prime_UTR","5' utr",
                                   ifelse(type=="three_prime_UTR","3' utr",as.character(type))))) %>%
  dplyr::filter(!(tair=="mRNA" | tair =="gene" | tair =="stop_codon" | tair =="start_codon")) %>%
  tidyr::separate(attributes,into = c("ID","Parent"), sep = ";") %>%
  dplyr::select(-ID,-type) %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(intron=ifelse(any(tair=="intron"),"wi","ni")) %>%
  dplyr::ungroup()

intronless_four <- four %>%
  dplyr::filter(intron=="ni" & tair == "exon") %>% 
  dplyr::group_by(Parent) %>%
  dplyr::mutate(intstart=end+1) %>%
  dplyr::mutate(intend=lead(start)-1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!(is.na(intend))) %>%
  dplyr::mutate(tair="intron") %>%
  dplyr::select(intstart,intend,Parent,tair) %>%
  dplyr::rename(start=intstart) %>%
  dplyr::rename(end=intend)
  
intronless_two <- two %>%
  dplyr::filter(intron=="ni" & tair == "exon") %>% 
  dplyr::group_by(Parent) %>%
  dplyr::mutate(intstart=end+1) %>%
  dplyr::mutate(intend=lead(start)-1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!(is.na(intend))) %>%
  dplyr::mutate(tair="intron") %>%
  dplyr::select(intstart,intend,Parent,tair) %>%
  dplyr::rename(start=intstart) %>%
  dplyr::rename(end=intend)

cds_four <- four %>%
  dplyr::filter(tair=="coding_region") %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(cdsstart=min(start)) %>%
  dplyr::mutate(cdsend=max(end)) %>%
  dplyr::distinct(cdsstart,cdsend) %>%
  dplyr::mutate(tair="ORF") %>%
  dplyr::select(cdsstart,cdsend,Parent,tair) %>%
  dplyr::ungroup() %>%
  dplyr::rename(start=cdsstart) %>%
  dplyr::rename(end=cdsend)

cds_two <- two %>%
  dplyr::filter(tair=="coding_region") %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(cdsstart=min(start)) %>%
  dplyr::mutate(cdsend=max(end)) %>%
  dplyr::distinct(cdsstart,cdsend) %>%
  dplyr::mutate(tair="ORF") %>%
  dplyr::select(cdsstart,cdsend,Parent,tair) %>%
  dplyr::ungroup() %>%
  dplyr::rename(start=cdsstart) %>%
  dplyr::rename(end=cdsend)

tair_four <- rbind(cds_four,four %>% dplyr::select(-intron),intronless_four) %>%
  dplyr::arrange(Parent) %>%
  dplyr::select(tair,start,end,Parent) %>%
  dplyr::rename(type=tair) %>%
  dplyr::mutate(start=start-fourgene$start) %>%
  dplyr::mutate(end=end-fourgene$start) %>%
  tidyr::unite("coordinates",start,end,sep = "-")

tair_two <- rbind(cds_two,two %>% dplyr::select(-intron),intronless_two) %>%
  dplyr::arrange(Parent) %>%
  dplyr::select(tair,start,end,Parent) %>%
  dplyr::rename(type=tair) %>%
  dplyr::mutate(start=start-fourgene$start) %>%
  dplyr::mutate(end=end-fourgene$start) %>%
  tidyr::unite("coordinates",start,end,sep = "-")

onegene <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/busco_dup_figure/one_new.gff") %>%
  dplyr::filter(type=="gene") %>%
  dplyr::select(start,end)

one <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/busco_dup_figure/one_new.gff") %>%
  dplyr::select(type,start,end,attributes) %>%
  dplyr::mutate(tair=ifelse(type=="CDS","coding_region",
                            ifelse(type=="five_prime_UTR","5' utr",
                                   ifelse(type=="three_prime_UTR","3' utr",as.character(type))))) %>%
  dplyr::filter(!(tair=="mRNA" | tair =="gene" | tair =="stop_codon" | tair =="start_codon")) %>%
  tidyr::separate(attributes,into = c("ID","Parent"), sep = ";") %>%
  dplyr::select(-ID,-type) %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(intron=ifelse(any(tair=="intron"),"wi","ni")) %>%
  dplyr::ungroup()


intronless_one <- one %>%
  dplyr::filter(intron=="ni" & tair == "exon") %>% 
  dplyr::mutate(intstart=end+1) %>%
  dplyr::mutate(intend=lead(start)-1) %>%
  dplyr::filter(!(is.na(intend))) %>%
  dplyr::mutate(tair="intron") %>%
  dplyr::select(intstart,intend,Parent,tair) %>%
  dplyr::rename(start=intstart) %>%
  dplyr::rename(end=intend)

cds_one <- one %>%
  dplyr::filter(tair=="coding_region") %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(cdsstart=min(start)) %>%
  dplyr::mutate(cdsend=max(end)) %>%
  dplyr::distinct(cdsstart,cdsend) %>%
  dplyr::mutate(tair="ORF") %>%
  dplyr::select(cdsstart,cdsend,Parent,tair) %>%
  dplyr::ungroup() %>%
  dplyr::rename(start=cdsstart) %>%
  dplyr::rename(end=cdsend)

tair_one <- rbind(cds_one,one %>% dplyr::select(-intron),intronless_one)  %>%
  dplyr::arrange(Parent) %>%
  dplyr::select(tair,start,end) %>%
  dplyr::rename(type=tair) %>%
  dplyr::arrange(start) %>%
  dplyr::mutate(start=start-onegene$start) %>%
  dplyr::mutate(end=end-onegene$start) %>%
  tidyr::unite("coordinates",start,end,sep = "-")



tair_four_clean <- rbind(tair_two,tair_four) %>%
  tidyr::separate(coordinates, into = c("start","end"), sep="-") %>%
  dplyr::filter(!(type=="ORF")) %>%
  dplyr::filter(!(type=="exon")) %>%
  dplyr::mutate(start=as.numeric(start)+1663398) %>%
  dplyr::mutate(end=as.numeric(end)+1663398) %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(ngroup=cur_group_id()*1.5) %>%
  dplyr::ungroup()%>%
  dplyr::mutate(ngroup=ifelse(ngroup==7.5,1.5,ifelse(ngroup==9,3,ngroup)))

tips <-tair_four_clean %>% 
  dplyr::filter(!(type=="intron")) %>%
  dplyr::group_by(Parent) %>%
  dplyr::arrange(desc(start)) %>%
  dplyr::filter(row_number()==n()) %>%
  dplyr::mutate(xpol=list(c((end-(end-start)/2),(end-(end-start)/2),start))) %>%
  dplyr::mutate(ypol=list(c(ngroup+0.5,ngroup-0.5,ngroup))) %>%
  dplyr::ungroup() 

ids <- factor(c("t1", "t2", "t3", "t4","t5","t6"))
t2 <- data.frame(x=unlist(dplyr::pull(tips,xpol)),y=unlist(dplyr::pull(tips,ypol)),z=rep(ids,each=3))



p1 <- ggplot() + geom_rect(data = tair_four_clean %>% 
                       dplyr::filter(!(type=="intron")) %>%
                       dplyr::group_by(Parent) %>%
                       dplyr::arrange(desc(start)) %>%
                       dplyr::filter(!(row_number()==n())) %>%
                       dplyr::ungroup(), aes(xmin = start ,xmax = end ,ymin =ngroup-0.5 , ymax= ngroup+0.5),fill="lightblue") +
           geom_rect(data=tair_four_clean %>% 
                       dplyr::filter(!(type=="intron")) %>%
                       dplyr::group_by(Parent) %>%
                       dplyr::arrange(desc(start)) %>%
                       dplyr::filter(row_number()==n()) %>%
                       dplyr::ungroup(), aes(xmin=end,xmax=end-(end-start)/2,ymin = ngroup-0.5 , ymax=ngroup+0.5),fill="lightblue") +
           geom_polygon(data = t2,aes(x=x,y=y,group=z),fill="lightblue") +
           geom_segment(data = tair_four_clean %>% dplyr::filter(type=="intron"), aes(x=start,xend=start+((end-start)/2),y=ngroup,yend=ngroup+0.5),color="darkgrey") +
           geom_segment(data = tair_four_clean %>% dplyr::filter(type=="intron"), aes(x=start+((end-start)/2),xend=end,y=ngroup+0.5,yend=ngroup),color="darkgrey") + 
           annotate("text",x=1670750,y=7,label="QX1410.11877 locus",family="helvetica") +         
           theme(axis.title  = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              panel.background = element_blank(),
              plot.title = element_text(hjust = 0,size = 12)) +
           xlim(1663000,1678500) + ggtitle("A")

tair_one$Parent <- "Parent=transcript:QX1410.11877.3"

tair_one_clean <- rbind(tair_one,tair_two) %>%
  tidyr::separate(coordinates, into = c("start","end"), sep="-") %>%
  dplyr::filter(!(type=="ORF")) %>%
  dplyr::filter(!(type=="exon")) %>%
  dplyr::mutate(start=as.numeric(start)+1663398) %>%
  dplyr::mutate(end=as.numeric(end)+1663398) %>%  
  dplyr::group_by(Parent) %>%
  dplyr::mutate(ngroup=cur_group_id()*1.5) %>%
  dplyr::ungroup()

tips <-tair_one_clean %>% 
  dplyr::filter(!(type=="intron")) %>%
  dplyr::group_by(Parent) %>%
  dplyr::arrange(desc(start)) %>%
  dplyr::filter(row_number()==n()) %>%
  dplyr::mutate(xpol=list(c((end-(end-start)/2),(end-(end-start)/2),start))) %>%
  dplyr::mutate(ngroup=ifelse(ngroup==4.5,3,ngroup)) %>%
  dplyr::mutate(ypol=list(c(ngroup+0.5,ngroup-0.5,ngroup))) %>%
  dplyr::ungroup() 

tair_one_clean2 <- tair_one_clean %>% dplyr::mutate(ngroup=ifelse(ngroup==4.5,3,ngroup))


ids <- factor(c("t1","t2","t3"))
t3 <- data.frame(x=unlist(dplyr::pull(tips,xpol)),y=unlist(dplyr::pull(tips,ypol)),z=rep(ids,each=3))



text1 <- "paste(italic(pat-4))"
p2 <- ggplot() + geom_rect(data = tair_one_clean2 %>% 
                             dplyr::filter(!(type=="intron")) %>%
                             dplyr::group_by(Parent) %>%
                             dplyr::arrange(desc(start)) %>%
                             dplyr::filter(!(row_number()==n())) %>%
                             dplyr::ungroup(), aes(xmin = start ,xmax = end ,ymin =ngroup-0.5 , ymax= ngroup+0.5),fill="pink") +
  geom_rect(data=tair_one_clean2 %>% 
              dplyr::filter(!(type=="intron")) %>%
              dplyr::group_by(Parent) %>%
              dplyr::arrange(desc(start)) %>%
              dplyr::filter(row_number()==n()) %>%
              dplyr::ungroup(), aes(xmin=end,xmax=end-(end-start)/2,ymin = ngroup-0.5 , ymax=ngroup+0.5),fill="pink") +
  geom_polygon(data = t3,aes(x=x,y=y,group=z),fill="pink") +
  geom_segment(data = tair_one_clean2 %>% dplyr::filter(type=="intron"), aes(x=start,xend=start+((end-start)/2),y=ngroup,yend=ngroup+0.5),color="darkgrey") +
  geom_segment(data = tair_one_clean2 %>% dplyr::filter(type=="intron"), aes(x=start+((end-start)/2),xend=end,y=ngroup+0.5,yend=ngroup),color="darkgrey") + 
  annotate(geom="text",x=1666500,y=4,label="italic('pat-4')",parse=T,family="helvetica",size=4) +
  annotate(geom="text",x=1675500,y=4,label="italic('attf-3')",parse=T,family="helvetica",size=4) +
  theme(axis.title.y  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        plot.title = element_text(hjust = 0,size = 12)) +
  xlim(1663000,1678500) + ggtitle("B") + xlab("Genomic position")



plot_grid(p1,p2,ncol=1,nrow=2,rel_heights = c(1.66,1))


