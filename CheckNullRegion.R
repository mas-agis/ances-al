#Checking whether AA exist in windows with no count of AA
check.region <- function(AncesAA,Region_WO_AA,cr) {
  Ancestral_Allele1 <- AncesAA
  main1 <- Region_WO_AA 
  str(Ancestral_Allele1)
  Ancestral_Allele1$Pos <- as.numeric(Ancestral_Allele1$Pos)
  Gabungan <- data.frame(Chr=0,Start=0,End=0,AA_Count=0,Actual_AA_Sites=0)
  Gabungan <- Gabungan[FALSE,]
  for (i in cr) {
    data_i <- filter(Ancestral_Allele1,Chr==i)
    data_j <- filter(main1,Chr==i)
    Hitung <- data_i %>% mutate(Cut = cut(Pos, breaks = seq(0, max(Pos, na.rm = T), by = 10000))) %>%
      group_by(Cut,.drop = TRUE) %>%
      summarise(AA_count = n())
    f <- str_split_fixed(Hitung$Cut, ",", 2)
    h <-as.data.frame((f))
    h$V1 <- as.character(h$V1)
    h$V2 <- as.character(h$V2)
    h$V1 <- substr(h$V1,2,nchar(h$V1))
    h$V2 = substr(h$V2,1,nchar(h$V2)-1)
    h$V1 <- as.integer(h$V1)
    h$V2 <- as.integer(h$V2)
    Hitung$Chr <- as.factor(i)
    Hitung <- cbind(Hitung,h)
    Gabung <- left_join(data_j,Hitung, by=c("Start"="V1"))
    Gabung <- Gabung[,c(1,2,3,4,7)]
    colnames(Gabung) <- c("Chr","Start","End","AA_Count","Actual_AA_Sites")
    Gabungan <- rbind(Gabungan,Gabung)
  }
  return(Gabungan)
}
