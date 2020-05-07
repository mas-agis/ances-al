persist.ancesA <- function(AncesAA,gr,n,cr,dip=2,scan=10000,o=0.1, p=0.01, q=0.001) {
  library(dplyr)
  library(ggplot2)
  library(stringr)
  options(warn=-1)
  treshold0 <- data.frame(Chr=0,bar1=0,bar2=0,bar3=0)
  treshold0 <- treshold0[FALSE,]
  akhir0 <- data.frame(Chr=0,Start=0,End=0,Window=0,Score=0,Ancestral_count=0)
  akhir0 <- akhir0[FALSE,]
  for (chr in cr) {
    #Input first data frequency file from dir
    Data1 <- read.table(paste0(gr,"_Chr_",chr,".frq"), header = T, sep = "\t", fill = T, col.names = paste0("V", seq_len(7)))
    #Changing col names
    colnames(Data1) <- c("Chr", "Pos", paste0(gr,"_Alleles"), paste0(gr, "_N_Chr"),paste0(gr, "_A1"), paste0(gr, "_A2"), paste0(gr, "_A3") )
    #Fill empty string with 'NA'
    Data1[Data1==""]<-NA
    Data1 <- subset.data.frame(Data1, Data1[,4]==dip*n) 
    
    Data1$a <- substring(Data1[,5], 3)
    Data1$b <- substring(Data1[,6], 3)
    Data1$c <- substring(Data1[,7], 3)
    
    i <- as.character(c("a","b","c"))
    Data1$g <- i[max.col(replace(Data1[i], is.na(Data1[i]), -Inf), "first")]
    
    Data1$AA <- ifelse(Data1$g=="a", substring(Data1[,5], 1, 1), 
                       ifelse(Data1$g=="b", substring(Data1[,6], 1, 1), substring(Data1[,7], 1, 1)))
    Data1$Freq <- apply(Data1[,8:10], 1, max, na.rm=TRUE) 
    Data1$Chr <- as.factor(as.character(Data1$Chr))
    Data1$Pos <- as.factor(as.character(Data1$Pos))
    Data1$AA <- (as.character(Data1$AA))
    Data1 <- subset(Data1, AA!="*")
    Data1$AA <- as.factor(Data1$AA)
    Data1 <- Data1[,c(1,2,3,4,12,13)] #keep only wanted columns
    b <- inner_join(AncesAA, Data1, by = c("Chr","Pos"), copy =FALSE, keep=T)
    b[,9] <- as.numeric(b[,9])
    b$a <- ifelse(b$AA.x==b$AA.y, b$Freq.y, 0)
    b$b <- ifelse(b$AA.x!=b$AA.y, 1-b$Freq.y, 0)
    b$ChangeOfAncestral <- ifelse(b$a!=0, b$a-b$Freq.x, b$b-b$Freq.x)
    b$Pos <- as.numeric(b$Pos)
    
    ###########################################################################################################################
    #For alleles with 0 changes
    b$Poin <- ifelse(b$ChangeOfAncestral==0,1,0)
    akhir1 <- b %>%  mutate(Cut = cut(Pos, breaks = seq(0, max(Pos, na.rm = T), by = scan))) %>% 
      group_by(Cut,Poin,.drop = FALSE) %>%
      summarise(Ancestral_count = n()) %>%
      filter(Poin == 1)
    akhir1 <- as.data.frame(akhir1)
    f <- str_split_fixed(akhir1$Cut, ",", 2)
    h <-as.data.frame((f))
    h$V1 <- as.character(h$V1)
    h$V2 <- as.character(h$V2)
    h$V1 <- substr(h$V1,2,nchar(h$V1))
    h$V2 = substr(h$V2,1,nchar(h$V2)-1)
    h$V1 <- as.integer(h$V1)
    h$V2 <- as.integer(h$V2)
    akhir1$Chr <- as.factor(chr)
    akhir1 <- cbind(akhir1,h)
    colnames(akhir1) <- c("Window","Score","Ancestral_count","Chr","Start","End")
    akhir1$Window <- as.factor(seq.int(nrow(akhir1))) 
    akhir1 <- akhir1[,c(4:6,1,2,3)]
    
    #Setting treshold
    bar1 <- quantile(akhir1$Ancestral_count, prob=1-(o/100))
    bar2 <- quantile(akhir1$Ancestral_count, prob=1-(p/100))
    bar3 <- quantile(akhir1$Ancestral_count, prob=1-(q/100))
    
    #For statistical purpose
    treshold <- data.frame(Chr=chr,bar1=bar1,bar2=bar2,bar3=bar3)
    treshold0 <- rbind(treshold0,treshold)
    akhir0 <- rbind(akhir0,akhir1)
  }
    return(list(AACounts= akhir0, treshold=treshold0))
  
}
  
