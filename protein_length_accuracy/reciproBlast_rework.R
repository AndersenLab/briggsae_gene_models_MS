args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library("Biostrings")
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(ape)
})

#set output file suffix
suffix <- "QX1410-MERGER-ND"
#set output directory
outdir <- "/projects/b1059/projects/Nicolas/c.briggsae/apollo-analysis/plot/"

#read C. briggsae (C.b.) protein FASTA into AA stringset object
pred <- readAAStringSet("/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/final_annotations/check_dup_prot/QX1410.SB.final_anno.dedup.renamed.prot.fa")

#read C. elegans (C.e.) protein FASTA into stringset object
elegansprot <- readAAStringSet("/projects/b1059/projects/Nicolas/c.elegans/N2/wormbase/WS279/c_elegans.PRJNA13758.WS279.protein_coding.prot.fa")

#read blast results (forward and reciprocal)
blast = read.table("/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/Plots/blast_out/forward/QX1410.SB.final_anno.dedup.renamed.prot.pb.out", header = FALSE, sep = "", dec = ".")
recipro = read.table("/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/Plots/blast_out/reciprocals/QX1410.SB.final_anno.dedup.renamed.prot.recipro.pb.out", header = FALSE, sep = "", dec = ".")
colnames(blast) <- c("query","subject","ident","length","mismatch","gapop","qstart","qend","sstart","send","eval","bitscore")
colnames(recipro) <- c("query","subject","ident","length","mismatch","gapop","qstart","qend","sstart","send","eval","bitscore")

#read C.b. GFF 
predGFF <- read.gff("/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/final_annotations/check_dup_prot/QX1410.SB.final_anno.dedup.renamed.gff")
predGFF <- predGFF %>% dplyr::filter(type=="mRNA") %>%
                       tidyr::separate(attributes,into=c("Transcript_ID","Transcript_parent","Name","Other"),sep=";",extra="merge")
predGFF$Transcript_ID <- gsub(".*=","",predGFF$Transcript_ID)
predGFF$Transcript_parent <- gsub(".*=","",predGFF$Transcript_parent)

#read C.e. GFF 
celeGFF <- read.gff("/projects/b1059/projects/Nicolas/c.elegans/N2/wormbase/WS279/c_elegans.PRJNA13758.WS279.protein_coding.gff")
celeGFF <- celeGFF %>% dplyr::filter(type=="mRNA") %>% 
                       tidyr::separate(attributes,into=c("Subject_ID","Subject_parent","Name","Other"),sep=";",extra="merge")
celeGFF$Subject_ID <- gsub(".*=","",celeGFF$Subject_ID)
celeGFF$Subject_parent <- gsub(".*=","",celeGFF$Subject_parent)


#transform AA stringset object to 2-column dataframe
prednames <- names(pred)
predseq <- paste(pred)
preddf <- data.frame(prednames,predseq)
eprotnames <- names(elegansprot)
eprotseq <- paste(elegansprot)
eprotdf <- data.frame(eprotnames,eprotseq)

#calculate protein length for each FASTA entry
preddf$qlen <- str_count(preddf$predseq)
eprotdf$slen <- str_count(eprotdf$eprotseq)


#select best hit(s) per transcript for C.b.
orderedhits <- blast[order(blast$query,-blast$bitscore, blast$eval),]
zeroeval <- orderedhits[orderedhits$eval == 0,]
nonzeroeval <- orderedhits[orderedhits$eval > 0,]
filterednzhits <- nonzeroeval[!duplicated(nonzeroeval$query),]
filteredhits <- bind_rows(zeroeval,filterednzhits)
filteredhits <- filteredhits[order(filteredhits$query),]
filteredhits_qse <- filteredhits[,c(1:2,11,12)]
colnames(filteredhits_qse) <- c('Transcript_ID','Subject_ID','eval','bitscore')


#select best hit(s) per transcript for C.e.
rorderedhits <- recipro[order(recipro$query, -recipro$bitscore, recipro$eval),]
rzeroeval <- rorderedhits[rorderedhits$eval == 0,]
rnonzeroeval <- rorderedhits[rorderedhits$eval > 0,]
rfilterednzhits <- rnonzeroeval[!duplicated(rnonzeroeval$query),]
rfilteredhits <- bind_rows(rzeroeval,rfilterednzhits)
rfilteredhits <- rfilteredhits[order(rfilteredhits$query),]
rfilteredhits_qse <- rfilteredhits[,c(2,1,11,12)]
colnames(rfilteredhits_qse) <- c('Transcript_ID','Subject_ID','eval','bitscore')

#append gene-transcript relationships for every BLAST hit
forwardHits<- dplyr::left_join(filteredhits_qse,predGFF %>% dplyr::select("Transcript_ID","Transcript_parent"),by="Transcript_ID",keep=F)
forwardHits<- dplyr::left_join(forwardHits,celeGFF %>% dplyr::select("Subject_ID","Subject_parent"),by="Subject_ID",keep=F)
forwardHits <- forwardHits[,c(6,2,5,1,3,4)]
reciproHits<- dplyr::left_join(rfilteredhits_qse,predGFF %>% dplyr::select("Transcript_ID","Transcript_parent"),by="Transcript_ID",keep=F)
reciproHits<- dplyr::left_join(reciproHits,celeGFF %>% dplyr::select("Subject_ID","Subject_parent"),by="Subject_ID",keep=F)
reciproHits <- reciproHits[,c(6,2,5,1,3,4)]

#merge forward and reciprocal hits
aggregates <- dplyr::bind_rows(forwardHits,reciproHits, .id = "set")
aggregates <- aggregates[order(aggregates$Transcript_parent),]

#identify which transcripts have matching forward and reciprocal hits
common <- aggregates %>% 
  group_by(Transcript_ID,Subject_ID) %>% 
  mutate(dupe = n()>1) %>%
  dplyr::ungroup()

#identify best hit per locus
common <- common[order(common$Transcript_parent,-common$bitscore,common$eval),]
common_dupes <- common[common$dupe == TRUE,]
common_ndg <- common_dupes[!duplicated(common_dupes$Transcript_parent),]
common_ndg <- common_ndg[!duplicated(common_ndg$Subject_parent),]
commonNoDup <- common_ndg[,1:6]

#append protein lengths to best hit per locus and estimate protein length ratio (C.b. protein length / C.e. protein length)
predhits <- dplyr::left_join(commonNoDup,preddf %>% dplyr::select(prednames,qlen),by=c("Transcript_ID"="prednames"))
predhits <- dplyr::left_join(predhits,eprotdf %>% dplyr::select(eprotnames,slen),by=c("Subject_ID"="eprotnames"))
predhits$delta <- (predhits$qlen / predhits$slen)
predhits$normDiff <- (predhits$qlen - predhits$slen) / predhits$slen

#get list of 1:1 reciprocal hits
ones <- predhits[predhits$qlen == predhits$slen,]
#get list of 5% off reciprocal hits
fives <- predhits[(abs((predhits$qlen - predhits$slen)/predhits$slen) < 0.05) & (abs((predhits$qlen - predhits$slen)/predhits$slen) > 0),]

#select relevant columns from final output
outFile <- predhits %>% select(Transcript_ID,Subject_ID,delta)
colnames(outFile) <- c("query","subject","ratio")

#store counts for each hit type
totalcount <- paste("Total Hits:",nrow(predhits),sep=" ")
onescount <- paste("1:1 Hits:",nrow(ones),sep=" ")
fivescount <- paste("5% Off Hits:",nrow(fives),sep=" ")

#set output file dir and name
filenameRatio <- paste0(outdir,"/ratio-",suffix)
filenameNormDiff <- paste0(outdir,"/normDiff-",suffix)
filenameReport <- paste0(outdir,"/report-",suffix)
filenamePreds <- paste0(outdir,"/predhits-",suffix)

#write output tables
write.table(predhits$delta, sep=",",  col.names=FALSE, row.names = FALSE,file = filenameRatio,quote = F)
write.table(predhits$normDiff, sep=",",  col.names=FALSE, row.names = FALSE,file = filenameNormDiff,quote = F)
write.table(totalcount, append=FALSE, col.names=FALSE, row.names = FALSE,file = filenameReport)
write.table(onescount, append=TRUE, col.names=FALSE, row.names = FALSE,file = filenameReport)
write.table(fivescount, append=TRUE, col.names=FALSE, row.names = FALSE,file = filenameReport)
write.table(predhits %>% dplyr::select(-set,-qlen,-slen) , sep="\t",quote = F,  col.names=T, row.names = FALSE,file = filenamePreds)
