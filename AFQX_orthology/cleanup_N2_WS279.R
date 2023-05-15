library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(zoo)

#read gff
ws279 <- read.gff("/projects/b1059/projects/Nicolas/c.elegans/N2/wormbase/WS279/c_elegans.PRJNA13758.WS279.annotations.gff3")

#extract WormBase features only
anno_WB_only <- ws279 %>% dplyr::filter(source=="WormBase") %>%
  dplyr::filter(grepl("gene|mRNA|exon|intron|CDS|start_codon|stop_codon|five_prime_UTR|three_prime_UTR",type))

#extract all L1s
L1_features <- anno_WB_only %>% dplyr::filter(type=="gene") 

#keep only protein coding L1s
L1_features_pc <- L1_features %>% dplyr::filter(grepl("biotype=protein_coding",attributes)) #%>%
  # tidyr::separate(attributes,c("Feature_ID","Name","Other"),sep=";",extra="merge") #<- this was used to match L2 features to its respective parent

#extract L2s and separate features
L2_features <- anno_WB_only %>% dplyr::filter(type=="mRNA") %>%
   tidyr::separate(attributes,c("Feature_ID","Parent_ID","Other"),sep=";",extra="merge") 
L2_features$Parent_ID <- gsub(".*=","ID=",L2_features$Parent_ID)

#L2_check <- L2_features %>% dplyr::filter(Parent_ID %in% L1_features_pc$Feature_ID) #<- this was used to check that mRNA parents are protein-coding genes (duh!)

#extract nonCDS L3 features
L3_features_nonCDS <- anno_WB_only %>% dplyr::filter(grepl("exon|intron|start_codon|stop_codon|five_prime_UTR|three_prime_UTR",type)) %>%
  tidyr::separate(attributes,c("Parent_ID","Other"),sep=";",extra="merge")
L3_features_nonCDS$Parent_ID <- gsub(".*=","ID=",L3_features_nonCDS$Parent_ID)


#L3_nonCDScheck <- L3_features_nonCDS %>% dplyr::filter(!Parent_ID %in% L2_features$Feature_ID) #<- this was used to check if there were L3 features that did not belong to mRNA
# filter out L3s not found in mRNAs
L3_features_nonCDS_clean <- L3_features_nonCDS %>% dplyr::filter(Parent_ID %in% L2_features$Feature_ID)


#extract CDS
L3_CDSonly <- anno_WB_only %>% dplyr::filter(type=="CDS") #%>%
#   tidyr::separate(attributes,c("Feature_ID","Parent_ID","Other"),sep=";",extra="merge")
# L3_CDSonly$Parent_ID <- gsub(".*=","ID=",L3_CDSonly$Parent_ID)
# L3_CDScheck <- L3_CDSonly %>% dplyr::filter(!Parent_ID %in% L2_features$Feature_ID) #<- this was used to check if CDS parents were mRNA (duh!)

#reset GFF format (reset Parent prefix and unite attributes)
L3_features_nonCDS_clean$Parent_ID <- gsub(".*=","Parent=",L3_features_nonCDS_clean$Parent_ID)
L3_features_nonCDS_cleanGFF <- L3_features_nonCDS_clean %>% tidyr::unite("attributes",Parent_ID,Other,sep=";",na.rm = T)
L2_features$Parent_ID <- gsub(".*=","Parent=",L2_features$Parent_ID)
L2_features_cleanGFF <- L2_features %>% tidyr::unite("attributes",Feature_ID,Parent_ID,Other,sep=";",na.rm = T)

#merge all features and sort by coordinate
clean_GFF <- bind_rows(L1_features_pc,L2_features_cleanGFF,L3_features_nonCDS_cleanGFF,L3_CDSonly) %>% dplyr::arrange(seqid,start)

#write out clean GFF
write.table(clean_GFF,file="/projects/b1059/projects/Nicolas/c.elegans/N2/wormbase/WS279/c_elegans.PRJNA13758.WS279.protein_coding.gff",quote = F,sep = "\t",row.names = F,col.names = F,na = ".")
