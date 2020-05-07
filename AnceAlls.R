
ancealls<- function(gr,n,cr,fr=1,dip=2) {
library(dplyr)
library(ggplot2)
library(stringr)
#######################################################################################################
#Getting ancestral alleles by all outgroup
Ancestral_Allele <- data.frame(Chr=0,Pos=0,Alleles_n=0,AA=0,Freq=0)
Ancestral_Allele <- Ancestral_Allele[FALSE,]
options(warn=-1)
for (chr in cr) {
  #Input first data frequency file from dir
  Data4 <- read.table(paste0(gr,"_Chr_",chr,".FRQ"), header = T, sep = "\t", fill = T, col.names = paste0("V", seq_len(7)))
  #Changing col names
  colnames(Data4) <- c("Chr", "Pos", paste0(gr,"_Alleles"), paste0(gr, "_N_Chr"),paste0(gr, "_A1"), paste0(gr, "_A2"), paste0(gr, "_A3") )
  #Fill empty string with 'NA'
  Data4[Data4==""]<-NA
  ########################################################################################################
  AA <- subset.data.frame(Data4, Data4[,4]==(n*dip),select = c(1,2,3,5,6,7))#Chr,Pos,paste0(gr,"_Alleles"),paste0(gr,"_A1"),paste0(gr,"_A2"),paste0(gr,"_A3")))
  AA$d <- substring(AA[,4], 3)
  AA$e <- substring(AA[,5], 3)
  AA$f <- substring(AA[,6], 3)
  k <- as.character(c("d","e","f"))
  AA$i <- k[max.col(replace(AA[k], is.na(AA[k]), -Inf), "first")]
  AA$AA <- ifelse(AA$i=="d", substring(AA[,4], 1, 1), 
                  ifelse(AA$i=="e", substring(AA[,5], 1, 1), substring(AA[,6], 1, 1)))
  AA$d <- as.numeric(AA$d)
  AA$e <- as.numeric(AA$e)
  AA$f <- as.numeric(AA$f)
  AA$Chr <- as.factor(as.character(AA$Chr))
  
  #Getting allele frequency for the AA 
  AA$Freq <- apply(AA[,7:9], 1, max, na.rm=TRUE) 
  AA <- subset.data.frame(AA,Freq>0.5,select = c(1,2,3,11,12))
  colnames(AA) <- c("Chr", "Pos", "Alleles_n","AA","Freq")
  AA <- na.omit(AA)
  AA <- subset(AA, AA!="*")
  #Ancestral_Allele <- AA #(Activate only for the first time, when setting up the first table)
  Ancestral_Allele <- rbind(Ancestral_Allele,AA)
}
#Filter ancestral alleles with 100% frequency
Ancestral_Allele1 <- subset.data.frame(Ancestral_Allele, Ancestral_Allele$Freq>=fr)
Summary_AA_Content <- data.frame(A=sum(str_count(Ancestral_Allele$AA, "A")),
                                 G=sum(str_count(Ancestral_Allele$AA, "G")), 
                                 T=sum(str_count(Ancestral_Allele$AA, "T")), 
                                 C=sum(str_count(Ancestral_Allele$AA, "C")))
return(list(Ances_Allele = Ancestral_Allele1, summary=Summary_AA_Content))
}

##################### End of Ancestral Alleles by All outgroup Species############################################
