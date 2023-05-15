library(dplyr)
library(tidyr)
library(readr)
library(Biostrings)
library(ape)
library(stringr)
library(lubridate)

##This script parses and classifies user-created curation descriptions from Apollo

raw_curated_descriptions<- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/backup_anno/Curation-VF-230214.gff3") %>%
  dplyr::filter(type=="gene") %>%
  tidyr::separate(attributes,into=c("preID","DescP"),sep=";description=") %>%
  tidyr::separate(DescP,into=c("description","Other"),sep=";ID") %>%
  dplyr::select(seqid,start,end,description) %>%
  dplyr::filter(!is.na(description)) %>% 
  dplyr::mutate(newdesc=ifelse(str_detect(description,"Low|low"),"LC",
                               ifelse(str_detect(description,"xon|erminal"),"ICR",
                                      ifelse(str_detect(description,"plit"),"GS",
                                             ifelse(str_detect(description,"usio"),"GF",
                                                    ifelse(str_detect(description,"imple|asic|High"),"UN",
                                                           ifelse(str_detect(description,"ultip|forms"),"MI",
                                                                             ifelse(str_detect(description,"ew"),"NG",
                                                                                    ifelse(str_detect(description,"UTR"),"UE",
                                                                                           ifelse(str_detect(description,"ncompl"),"ICR",
                                                                                                  ifelse(str_detect(description,"raker"),"NB",
                                                                                                         ifelse(str_detect(description,"RAKER|REAKER"),"NB",
                                                                                                                ifelse(str_detect(description,"Imple"),"UN",
                                                                                                                       ifelse(str_detect(description,"xten"),"UE","NaN"))))))))))))))

unk_descriptions<- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/backup_anno/Curation-VF-230214.gff3") %>%
  dplyr::filter(type=="gene") %>%
  tidyr::separate(attributes,into=c("preID","DescP"),sep=";description=") %>%
  tidyr::separate(DescP,into=c("description","Other"),sep=";ID") %>%
  dplyr::select(seqid,start,end,description) %>%
  dplyr::filter(is.na(description))


c1 <- nrow(unk_descriptions)+2
c2 <- nrow(raw_curated_descriptions %>% dplyr::filter(newdesc=="ICR"))
c3 <- nrow(raw_curated_descriptions %>% dplyr::filter(newdesc=="GF"))
c4 <- nrow(raw_curated_descriptions %>% dplyr::filter(newdesc=="GS"))
c5 <- nrow(raw_curated_descriptions %>% dplyr::filter(newdesc=="NB"))
c6 <- nrow(raw_curated_descriptions %>% dplyr::filter(newdesc=="UN"))
c7 <- nrow(raw_curated_descriptions %>% dplyr::filter(newdesc=="UE"))
c8 <- nrow(raw_curated_descriptions %>% dplyr::filter(newdesc=="MI"))
c9 <- nrow(raw_curated_descriptions %>% dplyr::filter(newdesc=="LC"))
c10 <-nrow(raw_curated_descriptions %>% dplyr::filter(newdesc=="NG"))


SE <- c2+c3+c4+c8
ICR <- c2+c8

raw_dates <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/backup_anno/Curation-VF-230214.gff3") %>%
  dplyr::filter(type=="gene") %>%
  tidyr::separate(attributes,into=c("Other", "date_created"),sep= ";date_creation=") %>%
  dplyr::mutate(day_num = as.numeric(difftime(date_created, min(date_created), units = "days"))) %>%
  dplyr::group_by(day_num) %>%
  dplyr::mutate(nGene=n()) %>%
  dplyr::distinct(day_num,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(day_num) %>%
  dplyr::mutate(cumSum=cumsum(nGene)) %>%
  tidyr::separate(date_created,into=c("year","month","day"),sep="-")%>%
  dplyr::filter(!(year==2023) & day_num <418)


ggplot(raw_dates,aes(x=day_num,y=cumSum)) + geom_smooth(span=0.14) + geom_vline(xintercept = 221,linetype='dotted',color="red") + xlim(0,400)+
  xlab("Days") + ylab("Number of genes curated")+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=14),
        axis.title=element_text(size=16),
        panel.grid.major = element_line(color="lightgrey",linetype = "dashed"))


QXC_genes_BED <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/clean_anno/Curation-VF-230214.PC.clean.renamed.longest.gff3") %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(seqid,start,end,attributes,score,strand) %>%
  dplyr::mutate(score = coalesce(score, 100))



QXA_genes_BED <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/final_annotations/check_dup_prot/QX1410.SB.final_anno.dedup.renamed.longest.gff") %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(seqid,start,end,attributes,score,strand) %>%
  dplyr::mutate(score = coalesce(score, 100))


write.table(QXC_genes_BED,file="/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/stats/intersect/QX_CUR_genes.bed",quote = F, sep="\t",col.names = F,row.names = F)

write.table(QXA_genes_BED,file="/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/stats/intersect/QX_AUTO_genes.bed",quote = F, sep="\t",col.names = F,row.names = F)


intersect_res <- read.delim("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/stats/intersect/intersect.out",header = F)
colnames(intersect_res) <- c("A_seqid","A_start","A_end","A_name","A_score","A_strand","C_seqid","C_start","C_end","C_name","C_score","C_strand","overlap")

splits <- intersect_res %>% 
  dplyr::mutate(A_len=abs(A_start-A_end)) %>%
  dplyr::mutate(perc_overlap=overlap/A_len) %>%
  tidyr::separate(A_name,into = c("A_ID","A_Parent","A_other"),sep=";",extra = 'merge') %>%
  dplyr::select(-A_other) %>%
  dplyr::filter(perc_overlap > 0.05) %>%
  tidyr::separate(C_name,into = c("C_ID","C_Parent"),sep=";") %>%
  dplyr::mutate(isF=ifelse(perc_overlap==1,1,0)) %>%
  dplyr::group_by(A_Parent) %>%
  dplyr::mutate(sumPO=sum(perc_overlap)) %>%
  dplyr::mutate(nGroup=n()) %>%
  dplyr::mutate(sumF=sum(isF)) %>%
  dplyr::distinct(A_Parent,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nGroup > 1 & sumF < 1)

fusion <- intersect_res %>% 
  dplyr::mutate(C_len=abs(C_start-C_end)) %>%
  dplyr::mutate(perc_overlap=overlap/C_len) %>%
  tidyr::separate(A_name,into = c("A_ID","A_Parent","A_other"),sep=";",extra = 'merge') %>%
  dplyr::select(-A_other) %>%
  dplyr::filter(perc_overlap > 0.05) %>%
  tidyr::separate(C_name,into = c("C_ID","C_Parent"),sep=";") %>%
  dplyr::mutate(isF=ifelse(perc_overlap==1,1,0)) %>%
  dplyr::group_by(C_Parent) %>%
  dplyr::mutate(sumPO=sum(perc_overlap)) %>%
  dplyr::mutate(nGroup=n()) %>%
  dplyr::mutate(sumF=sum(isF)) %>%
  dplyr::distinct(C_Parent,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nGroup > 1 & sumF < 1)

isoforms <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/clean_anno/Curation-VF-230214.PC.clean.renamed.gff3") %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("ID","Parent"),sep = ";") %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(nGroup=n()) %>%
  dplyr::distinct(Parent,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nGroup > 1) %>%
  dplyr::mutate(isoN=nGroup-1)

sum(isoforms$isoN)
nrow(isoforms)

L1L2_relationships <-  ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/clean_anno/Curation-VF-230214.PC.clean.renamed.gff3") %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(attributes) %>%
  tidyr::separate(attributes, into=c("ID_tran","ID_gene"),sep=";")

utrs <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/clean_anno/Curation-VF-230214.PC.clean.renamed.gff3") %>%
  dplyr::filter(type=="five_prime_UTR" | type=="three_prime_UTR") %>%
  tidyr::separate(attributes, into=c("ID_utr","ID_tran"),sep=";") %>%
  dplyr::mutate(ID_tran=gsub("Parent=","ID=",ID_tran)) %>%
  dplyr::left_join(.,L1L2_relationships,by="ID_tran")

fiveP <- utrs %>%
  dplyr::filter(type=="five_prime_UTR") %>%
  dplyr::group_by(ID_gene) %>%
  dplyr::distinct(ID_gene,.keep_all = T) %>%
  dplyr::ungroup()

threeP <- utrs %>%
  dplyr::filter(type=="three_prime_UTR") %>%
  dplyr::group_by(ID_gene) %>%
  dplyr::distinct(ID_gene,.keep_all = T) %>%
  dplyr::ungroup()
