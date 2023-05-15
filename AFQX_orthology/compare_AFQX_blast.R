library(dplyr)
library(tidyr)
library(readr)
library(Biostrings)
library(ape)
library(ggplot2)
library(stringr)
library(ggsci)
library(cowplot)
library(grid)
library(gridExtra)

#read curated QX1410 protein-length accuracy values
curated <- readr::read_tsv("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/plot/predhits-CURATION-VF")

#read AF16 protein-length accuracy values
af16 <- readr::read_tsv("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/plot/predhits-AF16-WS280-ND")

#read AF16 protein sequence FASTA
af16_seq <- Biostrings::readAAStringSet("/projects/b1059/projects/Nicolas/c.briggsae/wormbase/WS280/c_briggsae.PRJNA10731.WS280.protein_coding.prot.fa") 

#read AF16 GFF
af16_gff <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/wormbase/WS280/c_briggsae.PRJNA10731.WS280.protein_coding.gff") %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("preID","Other"),sep=";Parent") %>%
  tidyr::separate(Other,into=c("Parent","Other"),sep=";Name") %>%
  dplyr::mutate(ID=gsub("ID=","",preID)) %>%
  dplyr::mutate(Gene=gsub("=Gene:","",Parent)) %>%
  dplyr::select(Gene,ID,seqid,start,end)

#read curated GFF and get transcript features
curated_gff <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/clean_anno/Curation-VF-230214.PC.clean.renamed.gff3") %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("preID","Other"),sep=";Parent") %>%
  dplyr::mutate(ID=gsub("ID=","",preID)) %>%
  dplyr::select(ID,seqid,start,end)
colnames(curated_gff) <- c("ID","QX1410_chr","QX1410_start","QX1410_end")

#read curated GFF and get gene features
predGFF <- read.gff("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/clean_anno/Curation-VF-230214.PC.clean.renamed.gff3")
predGenes <- predGFF %>% dplyr::filter(type=="gene")
predGenes$attributes <- gsub(".*=","",predGenes$attributes)
predGenes$len <- predGenes$end - predGenes$start

#read curated GFF and get exon features
predExon <- predGFF %>% dplyr::filter(type=="exon") %>%
  tidyr::separate(attributes,into=c("Exon_ID","Exon_parent"),sep=";",extra="merge")
predExon$Exon_ID <- gsub(".*=","",predExon$Exon_ID)
predExon$len <- predExon$end - predExon$start
predBED <- predExon %>% dplyr::select(seqid,start,end) %>%
  dplyr::arrange(match(seqid,c("I","II","III","IV","V","X")),start)

#write exons to BED file for featureCounts (coverage)
write.table(predBED,file="/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/featureCounts/bedtools/exon_pos.bed",quote = F, sep="\t",col.names = F,row.names = F)

#transfor AF16 FASTA into 2-column data frame of protein lengths
cbprotnames <- names(af16_seq)
cbprotseq <- paste(af16_seq)
cbprotdf <- data.frame(cbprotnames,cbprotseq) %>% 
  dplyr::mutate(first=substr(cbprotseq,1,1)) %>%
  dplyr::mutate(len=(nchar(cbprotseq))) %>%
  dplyr::mutate(transcript=gsub(":","_",cbprotnames)) %>%
  dplyr::select(transcript,first,len) 

#read C. briggsae QX1410 (C.b.) protein FASTA into AA stringset object
pred <- readAAStringSet("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/blastp/Curation-VF-230214.PC.clean.renamed.prot.fa")

#read C. elegans N2 (C.e.) protein FASTA into stringset object
elegansprot <- readAAStringSet("/projects/b1059/projects/Nicolas/c.elegans/N2/wormbase/WS279/c_elegans.PRJNA13758.WS279.protein_coding.prot.fa")

#transfor QX1410 FASTA into 2-column data frame of protein lengths
prednames <- names(pred)
predseq <- paste(pred)
preddf <- data.frame(prednames,predseq) %>%
  dplyr::mutate(first=substr(predseq,1,1)) %>%
  dplyr::mutate(len=(nchar(predseq))) %>%
  dplyr::mutate(transcript=gsub(":","_",prednames)) %>%
  dplyr::select(transcript,first,len) 

#transfor N2 FASTA into 2-column data frame of protein lengths
eprotnames <- names(elegansprot)
eprotseq <- paste(elegansprot)
eprotdf <- data.frame(eprotnames,eprotseq) %>%
  dplyr::mutate(first=substr(eprotseq,1,1)) %>%
  dplyr::mutate(len=(nchar(eprotseq))) %>%
  dplyr::mutate(transcript=gsub(":","_",eprotnames)) %>%
  dplyr::select(transcript,first,len) 

#aggregate protein-length accuracy values of gene models with a shared ortholog between QX1410 and AF16
comp2 <- dplyr::left_join(af16,curated,by="Subject_parent") %>% 
  dplyr::select(Subject_parent, Transcript_ID.x,Transcript_ID.y,delta.x,delta.y,eval.x,eval.y) %>%
  dplyr::filter(!is.na(Transcript_ID.y)) %>%
  dplyr::mutate(perc_PLchange=round((1-delta.x)-(1-delta.y),digits=4)*100) %>%
  dplyr::mutate(oppDir=ifelse((delta.x>1.05 & delta.y<0.95) | (delta.x<0.95 & delta.y>1.05)  ,"OPP","NOPP")) %>%
  dplyr::left_join(.,curated_gff,by=c("Transcript_ID.y"="ID"))

#classify as revised or unrevised based on PL accuracy cutoffs
categ_prob <- comp2 %>% 
  dplyr::mutate(color=ifelse( (perc_PLchange < (-10) | perc_PLchange > 10 ) & (delta.y < 0.9 | delta.y > 1.1) ,"REVISED","UNREVISED")) %>%
  dplyr::arrange(QX1410_chr, perc_PLchange)

#read exon coverage values
exonCov <- read.table("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/featureCounts/bedtools/bedcov.oe", sep="\t",header = F) 
colnames(exonCov) <- c("seqid","start","end","counts","base_matches","length","base_fraction")
exonCovDF <- exonCov %>% tidyr::unite("coords",seqid,start,end,sep = "-")

#read PLacc values with gene tag
predWtag <- read.delim("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/plot/predhits-CURATION-VF",header = T) %>% dplyr::select(Transcript_parent,delta)
colnames(predWtag) <- c("gene","delta")

#aggregate exon counts by gene and calculate RPK
predExonDF <- predExon %>% dplyr::select(Exon_ID,seqid,start,end) %>% 
  tidyr::unite("coords",seqid,start,end,sep = "-") %>%
  dplyr::left_join(.,exonCovDF,by="coords") %>%
  tidyr::separate(coords,into=c("seqid","start","end"),sep="-") %>% 
  tidyr::separate(Exon_ID, into = c("strain","geneid","L2"),sep="\\.", remove = F) %>%
  tidyr::separate(L2, into=c("tranid","exonNo"),sep="-") %>%
  tidyr::unite("transcript",geneid,tranid,sep=".",remove=F) %>%
  dplyr::mutate(eCov=counts*150/length) %>%
  dplyr::group_by(transcript) %>%
  #dplyr::mutate(groupn=n()) %>%
  dplyr::mutate(sumLen=sum(length)) %>%
  dplyr::mutate(sumBase_match=sum(base_matches)) %>%
  dplyr::mutate(averageDepth=mean(counts)) %>%
  dplyr::mutate(RPK=sum(counts)/(sumLen/1000)) %>%
  dplyr::mutate(coverage=mean(counts)*(150)/sumLen)%>%
  dplyr::mutate(avgCov=mean(eCov)) %>%
  dplyr::mutate(totalCounts=sum(counts)) %>%
  dplyr::mutate(avgCounts=mean(counts)) %>%
  dplyr::distinct(transcript,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(geneid) %>%
  dplyr::arrange(desc(sumBase_match),.by_group = T) %>%
  dplyr::distinct(geneid,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percCov=sumBase_match/sumLen) %>%
  tidyr::unite("gene",strain,geneid,sep=".") %>%
  dplyr::left_join(.,predWtag,by="gene") %>%
  dplyr::select(gene,percCov,averageDepth,RPK,delta,coverage,totalCounts,avgCounts,sumLen,avgCov) %>% 
  dplyr::filter(!is.na(delta))

#calculate TPM
perMil <- sum(predExonDF$RPK)/1000000
predExonDF$TPM <- predExonDF$RPK/perMil

#add common gene tag to dataframe containing PLacc values
categ_probWtag <- categ_prob %>% tidyr::separate(Transcript_ID.y, into=c("strain","gene","tranid"),sep="\\.",remove = F) %>% 
  dplyr::select(-tranid) %>%
  tidyr::unite("QX1410_gID",strain,gene,sep=".",remove = T) 

#join PLacc and TPM values
categ_probFinal <- dplyr::left_join(categ_probWtag,(predExonDF %>% dplyr::select(-delta)),by=c("QX1410_gID"="gene")) %>% dplyr::mutate(covtag=ifelse(coverage<1,"low","high"))

#plot revised and unrevised genes protein length accuracy (PLacc) values (Fig 5)
ggplot(categ_probFinal %>% dplyr::mutate(covtag=ifelse(coverage<1,"low","high"))) + geom_point(aes(x=delta.x,y=delta.y,color=color),shape=19) +
  scale_x_continuous(name="AF16-N2 protein length accuracy",limits=c(0,4)) +
  scale_y_continuous(name="QX1410-N2 protein length accuracy",limits=c(0,4)) +
  theme(legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=14),
        axis.title=element_text(size=16),
        panel.grid.major = element_line(color="lightgrey",linetype = "dashed"),
        legend.position = c(.91, .23),
        legend.justification = c("right", "top"), 
        legend.box.just = "right", 
        legend.background = element_blank(),
        legend.box.background=element_rect(color ="lightgrey", fill ="white"),
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=14),
        legend.key = element_rect(color="white",fill="white")) + scale_colour_manual(values=c("red","darkgoldenrod1")) 

#plot TPM x PLacc (Supplemental Fig 3)
ggplot(predExonDF) + geom_point(aes(x=TPM,y=abs(1-delta)*100),shape=19) + ylab("Percent protein sequence length difference from N2") + xlab("Transcripts per million") + xlim(0,1000) +
  theme(legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=14),
        axis.title=element_text(size=16),
        panel.grid.major = element_line(color="lightgrey",linetype = "dashed"),
        legend.position = c(.91, .20),
        legend.justification = c("right", "top"), 
        legend.box.just = "right", 
        legend.background = element_blank(),
        legend.box.background=element_rect(color ="lightgrey", fill ="white"),
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=14),
        legend.key = element_rect(color="white",fill="white")) + scale_colour_manual(values=c("red","darkgoldenrod1")) 

#get counts of genes that are over 10% off from N2 length, but under 20 TPM
all10 <- nrow(predExonDF %>% dplyr::filter((delta > 1.1 | delta <0.9)))
u20_all10 <- nrow(predExonDF %>% dplyr::filter(TPM<20 & (delta > 1.1 | delta <0.9)))
o20_all10 <- nrow(predExonDF %>% dplyr::filter(TPM>20 & (delta > 1.1 | delta <0.9)))

#print revision set
revision_set <- categ_probFinal %>% dplyr::filter(color=="REVISED") %>% dplyr::select(Transcript_ID.y,QX1410_chr,QX1410_start,QX1410_end,color,TPM)
colnames(revision_set) <- c("ID","Chromosome","Start","End","Status","TPM")
write.table(revision_set,"/projects/b1059/projects/Nicolas/github/revised_genes.tsv",quote = F, sep="\t",col.names = F,row.names = F)

#read orthofinder data
orthogroups <- readr::read_tsv("/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/orthology/VF/OrthoFinder/Results_Feb14/Orthogroups/Orthogroups.tsv")
colnames(orthogroups) <- c("Orthogroup","QX1410","AF16","N2")

#identify orthogroups with multiple genes
orthogroups_wCounts <- orthogroups %>% dplyr::mutate(nQX=str_count(QX1410,",")+1) %>%
  dplyr::mutate(nAF=str_count(AF16,",")+1) %>%
  dplyr::mutate(nN2=str_count(N2,",")+1) 
  
#identify genes in AF16 missing in QX1410
N2_missing_QX <- orthogroups %>% dplyr::filter(is.na(QX1410) & (!is.na(AF16) & !is.na(`N2`))) %>%
  dplyr::mutate(groupS_AF=str_count(AF16,",")+1) %>%
  dplyr::mutate(groupS_N2=str_count(N2,",")+1) %>%
  dplyr::filter(groupS_N2 ==1 & groupS_AF == 1) %>%
  dplyr::select(-groupS_N2,-groupS_AF) %>%
  dplyr::left_join(.,cbprotdf,by=c("AF16"="transcript")) %>%
  dplyr::select(-first) %>%
  dplyr::rename(AF16_len = len) %>%
  dplyr::left_join(.,eprotdf,by=c("N2"="transcript")) %>%
  dplyr::select(-first) %>%
  dplyr::rename(N2_len = len) %>% 
  dplyr::mutate(rat=AF16_len/N2_len)


#identify genes in QX1410 missing in AF16
N2_missing_AF <- orthogroups %>% dplyr::filter(!is.na(QX1410) & (is.na(AF16) & !is.na(`N2`))) %>%
  dplyr::mutate(groupS_QX=str_count(QX1410,",")+1) %>%
  dplyr::mutate(groupS_N2=str_count(N2,",")+1) %>%
  dplyr::filter(groupS_N2 ==1 & groupS_QX == 1) %>%
  dplyr::select(-groupS_N2,-groupS_QX) %>%
  dplyr::left_join(.,preddf,by=c("QX1410"="transcript")) %>%
  dplyr::select(-first) %>%
  dplyr::rename(QX_len = len) %>%
  dplyr::left_join(.,eprotdf,by=c("N2"="transcript")) %>%
  dplyr::select(-first) %>%
  dplyr::rename(N2_len = len) %>% 
  dplyr::mutate(rat=QX_len/N2_len)

#get transcript IDs for AF16
af16_coords <- af16_gff %>% 
  dplyr::mutate(ID=gsub(":","_",ID))
colnames(af16_coords) <- c("AF16_parent","AF16_ID","AF16_chr","AF16_start","AF16_end")

#identify 2:1:1 orthologs
expansions_1tm <- orthogroups_wCounts %>% 
  dplyr::filter(nN2 == 1 & nAF == 2 & nQX == 1) %>%
  dplyr::left_join(.,curated_gff,by=c("QX1410"="ID")) %>%
  tidyr::unite("coords1",QX1410_chr,QX1410_start,sep = ":") %>%
  tidyr::unite("coords",coords1,QX1410_end,sep = "..") 

#unnest expanded orthologous groups
expansions_1tm_unnested <- expansions_1tm %>%
  dplyr::mutate(AF16 = strsplit(as.character(AF16), ", ")) %>% 
  tidyr::unnest(AF16) %>%
  dplyr::left_join(.,af16_coords,by=c("AF16"="AF16_ID")) %>%
  dplyr::group_by(QX1410) %>%
  dplyr::arrange(desc(AF16_start)) %>%
  dplyr::mutate(distance=ifelse(AF16_chr==lead(AF16_chr),(AF16_start-lead(AF16_start))/1000,NA)) %>%
  dplyr::mutate(diffTag=ifelse(AF16_chr==lead(AF16_chr),"SC",paste0("DC",", ",lead(AF16_chr)))) %>%
  dplyr::ungroup()

#identify genes unique to QX1410
uniqueQX <- orthogroups %>% dplyr::filter(!is.na(QX1410) & (is.na(AF16) & is.na(`N2`)))
#identify genes unique to AF16
uniqueAF <- orthogroups %>% dplyr::filter(is.na(QX1410) & (!is.na(AF16) & is.na(`N2`)))

#set relationships between AF16 and QX1410 gene models
AFQX_relationships <- orthogroups_wCounts %>% dplyr::select(QX1410,AF16,nQX,nAF) %>%
  dplyr::filter(!is.na(QX1410) & !is.na(AF16)) 

#get parent-child relationships of AF16 gene and transcript GFF features
AF_L1L2 <- af16_gff %>% dplyr::mutate(ID2=gsub(":","_",ID)) %>% dplyr::select(Gene,ID2)

#get 1:1 orthologs between AF16 and QX1410
N2QX_11 <- orthogroups_wCounts %>% dplyr::filter(nN2 == 1 & nQX == 1& !(nAF==1))
N2AF_11 <- orthogroups_wCounts %>% dplyr::filter(nN2 == 1 & nAF == 1 & !(nQX==1))
AFQX_11 <- AFQX_relationships %>% dplyr::filter(nAF == 1 & nQX == 1) %>%
  dplyr::left_join(.,categ_probFinal[,1:4],by=c("QX1410"="Transcript_ID.y")) %>%
  dplyr::mutate(AF16_ss=gsub(":","_",Transcript_ID.x)) %>%
  dplyr::left_join(.,AF_L1L2,by=c("AF16"="ID2")) %>%
  dplyr::rename(.,Parent = Gene) %>%
  dplyr::left_join(.,AF_L1L2,by=c("AF16_ss"="ID2")) %>%
  dplyr::rename(.,Parent_ss = Gene) %>%
  dplyr::select(QX1410_gID,QX1410,Parent,Parent_ss,AF16,AF16_ss) %>%
  dplyr::mutate(concordanceT=ifelse(Parent==Parent_ss,"C","NC")) 

#get 1:1:1 relationships between QX1410, N2, and AF16
all_111 <-orthogroups_wCounts %>% dplyr::filter(nAF == 1 & nQX == 1 & nN2==1)
all_111_relationships <- all_111 %>% 
  dplyr::left_join(.,categ_probFinal %>% dplyr::select(Transcript_ID.x,Transcript_ID.y),by=c("QX1410"="Transcript_ID.y")) %>%
  dplyr::mutate(AF16_ss=gsub(":","_",Transcript_ID.x)) %>% 
  dplyr::select(-Transcript_ID.x) %>%
  dplyr::left_join(.,AF_L1L2,by=c("AF16"="ID2")) %>%
  dplyr::rename(.,Parent = Gene) %>%
  dplyr::left_join(.,AF_L1L2,by=c("AF16_ss"="ID2")) %>%
  dplyr::rename(.,Parent_ss = Gene) %>%
  dplyr::mutate(concordanceT=ifelse(Parent==Parent_ss,"C","NC")) %>%
  dplyr::left_join(.,preddf,by=c("QX1410"="transcript")) %>%
  dplyr::rename(.,QX_len = len) %>%
  dplyr::left_join(.,eprotdf,by=c("N2"="transcript")) %>%
  dplyr::rename(.,N2_len = len) %>%
  dplyr::left_join(.,cbprotdf,by=c("AF16"="transcript")) %>%
  dplyr::rename(.,AF_len = len) %>%
  dplyr::mutate(QXN2_ratio=QX_len/N2_len) %>%
  dplyr::mutate(AFN2_ratio=AF_len/N2_len)

#use orthology relationships to build new QX1410 and AF16 protein-length accuracy histogram data
QXN2_rat <- data.frame(as.table(as.matrix(all_111_relationships %>% dplyr::select(QXN2_ratio))))
QXN2_rat <- QXN2_rat[,2:3]
AFN2_rat <- data.frame(as.table(as.matrix(all_111_relationships %>% dplyr::select(AFN2_ratio))))
AFN2_rat <- AFN2_rat[,2:3]
h280_AF16 <- hist(AFN2_rat$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
d280_AF16 <- data.frame(x = h280_AF16$breaks,y = c(h280_AF16$counts,NA))
hm_QX <- hist(QXN2_rat$Freq,breaks=seq(0,200,by=0.02),plot = FALSE)
dm_QX <- data.frame(x = hm_QX$breaks,y = c(hm_QX$counts,NA))
dbind <- bind_rows(list(`QX1410 (CURATED)`=dm_QX,`AF16 (WS280)`=d280_AF16), .id="Query")
dbind$Query <- factor(dbind$Query, levels=c("AF16 (WS280)","QX1410 (CURATED)")) 

#plot histogram (Supplemental Figure 4)
bot <- ggplot() + 
  geom_step(data = dbind,aes(x = x,y = y,color=Query),stat = "identity") +
  xlim(0,2) + coord_cartesian(ylim = c(0, 1500), xlim = c(0.5, 1.5)) + 
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=14),
        axis.title.x=element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + scale_colour_manual(values=c("grey","blue","chartreuse2","darkgoldenrod1","red"))

top <- ggplot() + 
  geom_step(data = dbind,aes(x = x,y = y,color=Query),stat = "identity") +
  xlim(0,2) + coord_cartesian(ylim = c(1550, 7000), xlim = c(0.5, 1.5)) +
  theme(axis.line.x=element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=14),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"), 
        legend.box.just = "right", legend.box.background=element_rect(colour ="white", fill ="grey"),
        legend.margin = margin(6, 6, 6, 6)) + scale_colour_manual(values=c("grey","blue","chartreuse2","darkgoldenrod1","red"))

figureQX <- plot_grid(top,bot,rel_heights = c(1,3),align = 'v',ncol=1)
y.grob <- grid::textGrob("counts", gp=gpar( col="black", fontsize=17), rot=90)
x.grob <- grid::textGrob("Query protein length / N2 protein length", gp=gpar(col="black", fontsize=17))
figure_final <- grid.arrange(arrangeGrob(figureQX,ncol=1, left = y.grob, bottom = x.grob))

#count number of identical and 5% off matches in new PLacc values derived from orthogroup sets
af_ones <- nrow(all_111_relationships %>% dplyr::filter(AFN2_ratio == 1))
af_fives <- nrow(all_111_relationships %>% dplyr::filter((AFN2_ratio < 1.05 | AFN2_ratio > 0.95) & !(AFN2_ratio == 1)))
qx_ones <- nrow(all_111_relationships %>% dplyr::filter(QXN2_ratio == 1))
qx_fives <- nrow(all_111_relationships %>% dplyr::filter((QXN2_ratio < 1.05 | QXN2_ratio > 0.95) & !(QXN2_ratio == 1)))

