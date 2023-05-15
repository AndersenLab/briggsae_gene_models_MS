library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(showtext)
library(ape)
showtext_auto()

transformed_coords <- readr::read_tsv("/projects/b1059/projects/Nicolas/c.briggsae/genome_alignments_AFvQX/old/AFQX/AF16_transformed") %>% 
  tidyr::separate(`[TAGS]`,into = c("QX1410_chr","AF16_chr"),sep="\t") %>%
  dplyr::mutate(misplaced=ifelse(QX1410_chr==AF16_chr,"CP","MISPLACED"))

cp_coords <- transformed_coords %>% dplyr::filter(misplaced=="CP")
misplaced <- transformed_coords %>% dplyr::filter(misplaced=="MISPLACED")
II_coords <- cp_coords %>% dplyr::filter(QX1410_chr == "II")

vector1 <- c(12415648,12424577,11935983,11946973)

vector_polg <- c(11778449,11788541)


p2 <- ggplot(II_coords) + geom_rect(aes(xmin=vector1[1],xmax=vector1[2],ymin=11540000,ymax=12220000),fill="lightpink") +
  geom_rect(aes(xmin=vector1[3],xmax=vector1[4],ymin=11540000,ymax=12220000),fill="lightpink") +
  geom_rect(aes(xmin=11830000,xmax=12530000,ymin=vector_polg[1],ymax=vector_polg[2]),fill="lightblue") +
  geom_rect(aes(xmin=vector1[1],xmax=vector1[2],ymin=vector_polg[1],ymax=vector_polg[2]),fill="#8077D5") +
  geom_rect(aes(xmin=vector1[3],xmax=vector1[4],ymin=vector_polg[1],ymax=vector_polg[2]),fill="#8077D5") +
  geom_segment(aes(x=`[S2]`,xend=`[E2]`,y=`[S1]`,yend=`[E1]`)) + coord_cartesian(xlim=c(11930000,12430000),ylim=c(11650000,12110000)) +
  annotate("text",x=12180000,y=11798495,size=4,label="paste('QX1410 ',italic('polg-1'))", color="blue", family="helvetica",fontface="bold", parse=T) +
  annotate("text",x=12405648,y=11859495+1e4,size=4,label="paste('AF16 ',italic('polg-1'))", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +
  annotate("text",x=11956973-3e4,y=11859495+1e4,size=4,label="paste('AF16 ','CBG30119.1')", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +xlab("AF16 physical position") + ylab("QX1410 physical position")+theme(panel.background = element_blank(),
                                                                                                                                                                                                                                axis.title = element_blank(),
                                                                                                                                                                                                                                panel.border = element_rect(fill = NA))

IV_coords <- cp_coords %>% dplyr::filter(QX1410_chr == "IV")
vector2 <- c(7225388,7228637,8365205,8368291)
vector_best7 <- c(8223266,8226494)


p1 <- ggplot(IV_coords) + geom_rect(aes(xmin=vector2[1],xmax=vector2[2],ymin=8023266,ymax=8426494),fill="lightpink") +
  geom_rect(aes(xmin=vector2[3],xmax=vector2[4],ymin=8023266,ymax=8426494),fill="lightpink") +
  geom_rect(aes(xmin=7025388,xmax=8568291,ymin=vector_best7[1],ymax=vector_best7[2]),fill="lightblue") +
  geom_rect(aes(xmin=vector2[1],xmax=vector2[2],ymin=vector_best7[1],ymax=vector_best7[2]),fill="#8077D5") +
  geom_rect(aes(xmin=vector2[3],xmax=vector2[4],ymin=vector_best7[1],ymax=vector_best7[2]),fill="#8077D5") +
  geom_segment(aes(x=`[S2]`,xend=`[E2]`,y=`[S1]`,yend=`[E1]`)) + coord_cartesian(xlim=c(min(vector2)-1e5,max(vector2)+5e4),ylim=c(vector_best7[1]-1e5,vector_best7[2]+1e5)) +
  annotate("text",x=min(vector2) + (max(vector2)-min(vector2))/2,y=vector_best7[2]+3e3,size=4,label="paste('QX1410 ',italic('best-7'))", color="blue", family="helvetica",fontface="bold", parse=T) +
  annotate("text",y=vector_best7[2]+4e4,x=vector2[1]-3e4,size=4,label="paste('AF16 ',italic('best-7.2'))", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +
  annotate("text",y=vector_best7[2]+4e4,x=vector2[3]-3e4,size=4,label="paste('AF16 ',italic('best-7.1'))", color="red", family="helvetica",fontface="bold", parse=T,angle=90) +xlab("AF16 physical position") + ylab("QX1410 physical position")+theme(panel.background = element_blank(),
                                                                                                                                                                                                                                                       axis.title = element_blank(),
                                                                                                                                                                                                                                                panel.border = element_rect(fill = NA))

y.grob <- textGrob("QX1410 genome physical position",
                   gp=gpar(family="helvetica",fontsize=14),rot=90)
x.grob <- textGrob("AF16 genome physical position",
                   gp=gpar(family="helvetica",fontsize=14))
a.grob <- textGrob("C", gp=gpar(family="helvetica",fontsize=18),hjust=17)
b.grob <- textGrob("D", gp=gpar(family="helvetica",fontsize=18),hjust=17)

p12 <- plot_grid(arrangeGrob(p2,top=a.grob),arrangeGrob(p1,top=b.grob),ncol=2,nrow=1)

grid.arrange(arrangeGrob(p12,bottom=x.grob,left=y.grob))




af1 <- c(12230000,12242000)
afg1 <- c(12234639,12235897)
afg2 <- c(12237087,12238007)

qx2 <- c(12415000,12432950)
qxg1<- c(12421686,12422950)

V_coords <- cp_coords %>% dplyr::filter(QX1410_chr == "V")

p3 <- ggplot(V_coords) +
  geom_rect(aes(xmin=af1[1]-1e4,xmax=af1[2]+1e4,ymin=qxg1[1],ymax=qxg1[2]),fill="lightblue")+
  geom_rect(aes(xmin=afg1[1],xmax=afg1[2],ymin=qx2[1]-1e4,ymax=qx2[2]+1e4),fill="pink") +
  geom_rect(aes(xmin=afg2[1],xmax=afg2[2],ymin=qx2[1]-1e4,ymax=qx2[2]+1e4),fill="pink") +
  geom_rect(aes(xmin=afg1[1],xmax=afg1[2],ymin=qxg1[1],ymax=qxg1[2]),fill="#8077D5") +
  geom_rect(aes(xmin=afg2[1],xmax=afg2[2],ymin=qxg1[1],ymax=qxg1[2]),fill="#8077D5") +
  geom_segment(aes(x=`[S2]`,xend=`[E2]`,y=`[S1]`,yend=`[E1]`)) + 
  coord_cartesian(xlim=c(af1[1],af1[2]-1e2),ylim=c(qx2[1],qx2[2])) + 
  annotate("text",x=afg1[1]-3e2 ,y=qxg1[2]+5e3 ,label="paste('AF16 ',italic('ubh-1'))",parse=T,angle=90,color='red',family="helvetica") +
  annotate("text",x=afg2[1]-3e2 ,y=qxg1[2]+5e3 ,label="AF16 CBG18955",angle=90,color='red',family="helvetica") +
  annotate("text",x=afg2[2]+3e3 ,y=qxg1[2]+5e2 ,label="paste('QX1410 ',italic('ubh-1'))",parse=T,color='blue',family="helvetica") +
  xlab("AF16 physical position") + 
  ylab("QX1410 physical position") +
  theme(panel.background = element_blank(),
              axis.title = element_blank(),
              panel.border = element_rect(fill = NA))
  

af11 <- c(6706000,6710000)
afg11 <- c(6707063,6707454)
afg22 <- c(6708660,6709051)

qx22 <- c(6696000,6702000)
qxg11<- c(6698086,6698601)

I_coords <- cp_coords %>% dplyr::filter(QX1410_chr == "I")

p4 <- ggplot(I_coords) +
  geom_rect(aes(xmin=af11[1]-1e4,xmax=af11[2]+1e4,ymin=qxg11[1],ymax=qxg11[2]),fill="lightblue")+
  geom_rect(aes(xmin=afg11[1],xmax=afg11[2],ymin=qx22[1]-1e4,ymax=qx22[2]+1e4),fill="pink") +
  geom_rect(aes(xmin=afg22[1],xmax=afg22[2],ymin=qx22[1]-1e4,ymax=qx22[2]+1e4),fill="pink") +
  geom_rect(aes(xmin=afg11[1],xmax=afg11[2],ymin=qxg11[1],ymax=qxg11[2]),fill="#8077D5") +
  geom_rect(aes(xmin=afg22[1],xmax=afg22[2],ymin=qxg11[1],ymax=qxg11[2]),fill="#8077D5") +
  geom_segment(aes(x=`[S2]`,xend=`[E2]`,y=`[S1]`,yend=`[E1]`)) + 
  coord_cartesian(xlim=c(af11[1],af11[2]),ylim=c(qx22[1],qx22[2])) + 
  annotate("text",x=afg11[1]-1e2 ,y=qxg11[2]+2e3 ,label="paste('AF16 ',italic('mzt-1.1'))",parse=T,angle=90,color='red',family="helvetica") +
  annotate("text",x=afg22[1]-1e2 ,y=qxg11[2]+2e3 ,label="paste('AF16 ',italic('mzt-1.2'))",parse=T,angle=90,color='red',family="helvetica") +
  annotate("text",x=afg22[2]+5e2 ,y=qxg11[2]+1e2 ,label="paste('QX1410 ',italic('mzt-1'))",parse=T,color='blue',family="helvetica") +
  xlab("AF16 physical position") + 
  ylab("QX1410 physical position") +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(fill = NA))

y.grob <- textGrob("QX1410 genome physical position",
                   gp=gpar(family="helvetica",fontsize=14),rot=90)
x.grob <- textGrob("AF16 genome physical position",
                   gp=gpar(family="helvetica",fontsize=14))
a.grob <- textGrob("A", gp=gpar(family="helvetica",fontsize=18),hjust=17)
b.grob <- textGrob("B", gp=gpar(family="helvetica",fontsize=18),hjust=20)

p34 <- plot_grid(arrangeGrob(p3,top=a.grob),arrangeGrob(p4,top=b.grob),ncol=2,nrow=1)

grid.arrange(arrangeGrob(p34,p12,bottom=x.grob,left=y.grob))

mpChrom <- misplaced %>%
  dplyr::filter(!grepl("cb",AF16_chr))

AF16_genelist <- ape::read.gff("/projects/b1059/projects/Nicolas/c.briggsae/wormbase/WS280/c_briggsae.PRJNA10731.WS280.protein_coding.gff") %>%
  dplyr::filter(type=="gene") %>%
  dplyr::select(seqid,start,end, attributes) %>%
  dplyr::filter(!grepl("cb",seqid)) %>%
  tidyr::separate(attributes, into = c("pre","post"), sep = ";sequence_name=") %>%
  tidyr::separate(post, into = c("name","post2"), sep = ";biotype") %>%
  tidyr::separate(post2, into = c("pre2","aliases"),sep = ";Alias=") %>%
  tidyr::separate(aliases, into= c("alias","alt"), sep=",") %>%
  dplyr::select(-pre,-pre2,-alt)

hits <- data.frame(seqid=factor(),start=integer(),end=integer())
colnames(hits) <- c("seqid","start","end")

## slow AF
for (i in 1:nrow(mpChrom)){
  hit <- AF16_genelist %>%
    dplyr::filter(seqid == mpChrom$AF16_chr[i] & start >= mpChrom$`[S2]`[i] & end <= mpChrom$`[E2]`[i])
  hits <- bind_rows(hits,hit)
}

clean_hits <- hits %>% tidyr::unite("coords",start,end,sep="-") %>%
  dplyr::distinct(coords,.keep_all = T)


duplications <- data.frame(S1=double(),E1=double(),S2=double(),E2=double(),chr=character(),IDY=double())
cp_coords_short <- cp_coords %>% dplyr::select(`[S1]`,`[E1]`,`[S2]`,`[E2]`,AF16_chr,`[% IDY]`)
colnames(cp_coords_short) <- c("S1","E1","S2","E2","chr","IDY")

#### VERY SLOW - 30 MINUTES
for (i in 1:nrow(cp_coords_short)) {
  glab <- paste0(cp_coords_short$chr[i],":",cp_coords_short$S1[i],"-",cp_coords_short$E1[i],"|",cp_coords_short$chr[i],":",cp_coords_short$S2[i],"-",cp_coords_short$E2[i])
  subset <- cp_coords_short %>% 
    dplyr::filter(!(S1 > cp_coords_short$E1[i] & S1 > cp_coords_short$S1[i] & E1 > cp_coords_short$E1[i] & E1 > cp_coords_short$S1[i]) & 
                  !(S1 < cp_coords_short$S1[i] & S1 < cp_coords_short$E1[i] & E1 < cp_coords_short$S1[i] & E1 < cp_coords_short$E1[i]) & 
                  !(S1==cp_coords_short$S1[i] & E1==cp_coords_short$E1[i]) &
                    chr == cp_coords_short$chr[i]) %>%
    dplyr::mutate(group=glab)
  #print(subset)
  duplications <- bind_rows(duplications,subset)
}

write_tsv(duplications, "/projects/b1059/projects/Nicolas/c.briggsae/genome_alignments_AFvQX/g2g_AF16_duplicated_contigs.tsv")

clean_dup <- duplications %>% tidyr::separate(group, into = c("SE1","SE2"), sep = "\\|") %>%
  tidyr::separate(SE1,into = c("Rchr","Rcoords"), sep = ":") %>%
  dplyr::select(-Rchr) %>%
  tidyr::separate(Rcoords, into = c("RS1","RE1"), sep="-",convert = T) %>%
  tidyr::separate(SE2, into = c("Rchr","Rcoords"), sep=":") %>%
  dplyr::select(-Rchr) %>%
  tidyr::separate(Rcoords, into = c("RS2","RE2"),sep="-",convert = T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sort_id=paste(sort(c(S1,E1,S2,E2,RS1,RE1,RS2,RE2)),collapse = "-")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sort_id) %>%
  #dplyr::mutate(ngroup=cur_group_id()) %>%
  #dplyr::mutate(gsize=n()) %>%
  dplyr::distinct(sort_id,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::select(-sort_id) #%>%
  #dplyr::mutate(fstring=paste(chr,S1,E1,S2,E2,sep="-")) %>%
  #dplyr::mutate(rstring=paste(chr,RS1,RE1,RS2,RE2,sep="-")) %>% 
  #dplyr::group_by(fstring) %>%
  #dplyr::mutate(ng1=cur_group_id(),ngs1=n()) %>%
  #dplyr::ungroup() %>%
  #dplyr::group_by(rstring) %>%
  #dplyr::mutate(ng2=cur_group_id(),ngs2=n()) %>%
  #dplyr::ungroup() %>%
  #dplyr::mutate(overlaps=ifelse(ngs1==1 & ngs2==1,"single","multi")) 

#simple_overlaps <- clean_dup %>% dplyr::filter(overlaps=="single")
#multi_overlaps <- clean_dup  %>% dplyr::filter(overlaps=="multi")


# clasif: [++ / top: S1 < RS1 < E1 | S1 < RE1 > E1 ] \\ form: RS2 : RS2 + (E1 - RS1)              | form2: E2 - (E1 - RS1) : E2
# clasif: [++ / con: S1 < RS1 < E1 | S1 < RE1 < E1 ] \\ form: RS2 : RE2                           | form2: S2 + (RS1 - S1) : E2 - (E1 - RE1)
# clasif: [++ / bot: S1 > RS1 < E1 | S1 < RE1 < E1 ] \\ form: RE2 - (RE1 - S1) : RE2              | form2: S2 : S2 + (RE1 - S1)
# clasif: [++ / cco: S1 > RS1 < E1 | S1 < RE1 > E1 ] \\ form: RS2 + (S1 - RS1) : RE2 - (RE1 - E1) | form2: S2 : E2

# clasif: [-+ / top: S1 < RS1 < E1 | S1 < RE1 > E1 ] \\ form: RS2 : RS2 + (E1 - RS1)              | form2: E2 : E2 + (E1 - RS1)
# clasif: [-+ / con: S1 < RS1 < E1 | S1 < RE1 < E1 ] \\ form: RS2 : RE2                           | form2: E2 + (E1 - RE2) : S2 - (RS1 - S1)
# clasif: [-+ / bot: S1 > RS1 < E1 | S1 < RE1 < E1 ] \\ form: RE2 - (RE1 - S1) : RE2              | form2: S2 - (RE1 - S1) : S2
# clasif: [-+ / cco: S1 > RS1 < E1 | S1 < RE1 > E1 ] \\ form: RS2 + (S1 - RS1) : RE2 - (RE1 - E1) | form2: E2 : S2

# clasif: [+- / top: S1 < RS1 < E1 | S1 < RE1 > E1 ] \\ form: RS2 - (E1 - RS1) : RS2              | form2: E2 - (E1 - RS1) : E2
# clasif: [+- / con: S1 < RS1 < E1 | S1 < RE1 < E1 ] \\ form: RE2 : RS2                           | form2: S2 + (RS1 - S1) : E2 - (E1 - RE1)
# clasif: [+- / bot: S1 > RS1 < E1 | S1 < RE1 < E1 ] \\ form: RE2 : RE2 + (RE1 - S1)              | form2: S2 : S2 + (RE1 - S1)
# clasif: [+- / cco: S1 > RS1 < E1 | S1 < RE1 > E1 ] \\ form: RE2 + (RE1 - E1) : RS2 - (S1 - RS1) | form2: S2 : E2

# clasif: [-- / top: S1 < RS1 < E1 | S1 < RE1 > E1 ] \\ form: RS2 - (E1 - RS1) : RS2              | form2: E2 : E2 + (E1 - RS1)
# clasif: [-- / con: S1 < RS1 < E1 | S1 < RE1 < E1 ] \\ form: RE2 : RS2                           | form2: E2 + (E1 - RE2) : S2 - (RS1 - S1)
# clasif: [-- / bot: S1 > RS1 < E1 | S1 < RE1 < E1 ] \\ form: RE2 : RE2 + (RE1 - S1)              | form2: S2 - (RE1 - S1) : S2
# clasif: [-- / cco: S1 > RS1 < E1 | S1 < RE1 > E1 ] \\ form: RE2 + (RE1 - E1) : RS2 - (S1 - RS1) | form2: E2 : S2



af16_overlaps <- clean_dup %>% 
  dplyr::mutate(slopeF=(E2-S2)/(E1-S1)) %>%
  dplyr::mutate(slopeR=(RE2-RS2)/(RE1-RS1)) %>%
  dplyr::mutate(spat=ifelse(slopeR > 0 & slopeF > 0,"PP",
                            ifelse(slopeR < 0 & slopeF < 0,"MM",
                                   ifelse(slopeR >0 & slopeF <0,"MP","PM")))) %>%
  #dplyr::select(-slopeR,-slopeF) %>%
  dplyr::mutate(position=ifelse(S1<RS1,ifelse(E1<RE1,"top","con"),ifelse(E1>RE1,"bot","cco"))) %>%
  dplyr::mutate(overlap=ifelse(slopeR>0,ifelse(position=="top",paste0(RS2,":",RS2+(E1-RS1)),
                                               ifelse(position=="con",paste0(RS2,":",RE2),
                                                      ifelse(position=="bot",paste0(RE2-(RE1-S1),":",RE2), 
                                                             paste0(RS2+(S1-RS1),":",RE2-(RE1-E1))))),
                                        ifelse(position=="top",paste(RS2-(E1-RS1),":",RS2),
                                               ifelse(position=="con",paste0(RE2:RS2),
                                                      ifelse(position=="bot",paste0(RE2,":",RE2+(RE1-S1)),
                                                             paste0(RE2+(RE1-E1),":",RS2-(S1-RS1))))))) %>%
  dplyr::mutate(overlap2=ifelse(slopeF>0,ifelse(position=="top",paste0(E2-(E1-RS1),":",E2),
                                               ifelse(position=="con",paste0(S2+(RS1-S1),":",E2-(E1-RE1)),
                                                      ifelse(position=="bot",paste0(S2,":",S2+(RE1-S1)), 
                                                             paste0(S2,":",E2)))),
                               ifelse(position=="top",paste(E2,":",E2+(E1-RS1)),
                                      ifelse(position=="con",paste0(E2+(E1-RE2),":",S2-(RS1-S1)),
                                             ifelse(position=="bot",paste0(S2-(RE1-S1),":",S2),
                                                    paste0(E2,":",S2)))))) %>%
  tidyr::separate(overlap,into=c("overlap_start","overlap_end"),sep=":",convert = T) %>%
  tidyr::separate(overlap2,into=c("overlap2_start","overlap2_end"),sep=":",convert = T) %>%
  dplyr::mutate(overlap_len=overlap2_end-overlap2_start) %>%
  dplyr::mutate(seg_over=ifelse((overlap2_start > overlap_start & overlap2_start < overlap_end) |
                                                                                         (overlap2_end > overlap_start & overlap2_end < overlap_end) |
                                                                                         (overlap_start > overlap2_start & overlap_start < overlap2_end) |
                                                                                         (overlap_end > overlap2_start & overlap_end < overlap2_end),"Y","N")) %>%
  dplyr::filter(seg_over=="N")
  

test <- af16_overlaps %>% dplyr::filter(overlap_start==12529827)
test
ggplot(test) + 
   geom_segment(aes(x=S2,xend=E2,y=S1,yend=E1),color="red") +
   geom_segment(aes(x=RS2,xend=RE2,y=RS1,yend=RE1),color="blue",linetype="dashed") +
   annotate("text",x=test$S2,y=test$S1,label=test$slopeF,color="red") +
   annotate("text",x=test$S2,y=test$S1+1e3,label=paste0(test$spat,"+",test$position), color="green") +
   annotate("text",x=test$RS2,y=test$RS1,label=test$slopeR,color="blue") + xlab("AF16") + ylab("QX1410") + theme_bw() +
   geom_vline(xintercept = test$overlap_start,color="lightblue") +
   geom_vline(xintercept = test$overlap_end,color="lightblue") +
   geom_vline(xintercept = test$overlap2_start,color="pink") +
   geom_vline(xintercept = test$overlap2_end,color="pink") 


overlap_coords <- af16_overlaps %>%
  dplyr::select(chr,overlap_start,overlap_end,overlap2_start,overlap2_end,IDY)

affected_genes_f1 <- data.frame(chr=character(),overlap_start=integer(),overlap_end=integer(),overlap2_start=integer(),overlap2_end=integer(),set1=character(),set2=character())
affected_genes_f2 <- data.frame(chr=character(),overlap_start=integer(),overlap_end=integer(),overlap2_start=integer(),overlap2_end=integer(),set1=character(),set2=character())

for (i in 1:nrow(overlap_coords)) {
  #print(i)
  reg_genes <- data.frame(chr=character(),overlap_start=integer(),overlap_end=integer(),overlap2_start=integer(),overlap2_end=integer(),set1=character(),set2=character())
  reg_genes_f2 <- data.frame(chr=character(),overlap_start=integer(),overlap_end=integer(),overlap2_start=integer(),overlap2_end=integer(),set1=character(),set2=character())
  reg <- overlap_coords[i,]
  genes_set1 <- AF16_genelist %>%
    dplyr::filter(start > reg$overlap_start & end < reg$overlap_end) %>%
    dplyr::select(alias) 
  genes_set2 <- AF16_genelist %>%
    dplyr::filter(start > reg$overlap2_start & end < reg$overlap2_end) %>%
    dplyr::select(alias)
  if (nrow(genes_set1)==0 & nrow(genes_set2)==0){
    next
  }
  if(nrow(genes_set1)>0 & nrow(genes_set2) == 0) {
    #print("X0")
    reg_genes <- cbind(reg,data.frame(set1=genes_set1$alias,set2=rep(NA,nrow(genes_set1))))
    reg_genes_f2 <- cbind(reg,data.frame(set1=paste(genes_set1$alias,collapse = ","),set2=NA))
  } else if (nrow(genes_set1)==0 & nrow(genes_set2) > 0) {
    #print("0X")
    reg_genes <- cbind(reg,data.frame(set1=rep(NA,nrow(genes_set2)),set2=genes_set2$alias))  
    reg_genes_f2 <- cbind(reg,data.frame(set1=NA,set2=paste(genes_set2$alias,collapse = ",")))
  } else if (nrow(genes_set1) > nrow(genes_set2)){
    #print("Xx")
    reg_genes <- cbind(reg,data.frame(set1=genes_set1$alias,set2=c(genes_set2$alias,rep(NA,nrow(genes_set1)-nrow(genes_set2)))))
    reg_genes_f2 <- cbind(reg,data.frame(set1=paste(genes_set1$alias,collapse = ","),set2=paste(genes_set2$alias,collapse = ",")))
  } else if (nrow(genes_set1) < nrow(genes_set2)){
    #print("xX")
    reg_genes <- cbind(reg,data.frame(set1=c(genes_set1$alias,rep(NA,nrow(genes_set2)-nrow(genes_set1))),set2=genes_set2$alias))
    reg_genes_f2 <- cbind(reg,data.frame(set1=paste(genes_set1$alias,collapse = ","),set2=paste(genes_set2$alias,collapse = ",")))
  } else {
    #print("XX")
    reg_genes <- cbind(reg,data.frame(set1=genes_set1$alias,set2=genes_set2$alias)) 
    reg_genes_f2 <- cbind(reg,data.frame(set1=paste(genes_set1$alias,collapse = ","),set2=paste(genes_set2$alias,collapse = ",")))
  }
  affected_genes_f1 <- rbind(affected_genes_f1,reg_genes)
  affected_genes_f2 <- rbind(affected_genes_f2,reg_genes_f2)
}

affected_genes_table <- affected_genes_f2 %>%
  dplyr::filter(!(is.na(set1)) & !(is.na(set2))) %>%
  tidyr::unite("preseg1",chr,overlap_start,sep=":",remove = F) %>%
  tidyr::unite("segment A",preseg1,overlap_end,sep="-") %>%
  dplyr::select(-overlap_start) %>%
  tidyr::unite("preseg2",chr,overlap2_start,sep=":") %>%
  tidyr::unite("segment B",preseg2,overlap2_end,sep="-") %>%
  dplyr::rename("overlap identity"=IDY) %>%
  dplyr::rename("genes in segement A"=set1) %>%
  dplyr::rename("genes in segement B"=set2)

write_delim(affected_genes_table, "/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/busco_dup_figure/genes_affected_byDup.tsv",delim="\t")


affected_genes_wsize <- affected_genes_f1 %>%
  dplyr::mutate(size1=overlap_end-overlap_start,size2=overlap2_end-overlap2_start) %>%
  dplyr::group_by(overlap_start) %>%
  dplyr::mutate(group=cur_group_id()) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(remove=ifelse(gsize==1 & (is.na(set1) | is.na(set2)),"Y","N")) %>%
  dplyr::filter(remove=="N")


set1 <- affected_genes$set1
set1 <- set1[!is.na(set1)]
set2 <- affected_genes$set2
set2 <- set2[!is.na(set2)]

set3 <- c(set1,set2)
set4 <- unique(set3)


