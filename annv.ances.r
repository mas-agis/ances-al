annv.ances<-function(AAcount, treshold, bar=bar1){
  #getting regions higher than bar1
  treshold$Chr <- as.factor(treshold$Chr)
  cobalagi<- left_join(AAcount,treshold,by="Chr")
  cobadehlagi<- subset.data.frame(cobalagi,subset= Ancestral_count>bar1)
  cobadehlagi$Ref <- 0 
  cobadehlagi$min <- "-"
  output <- cobadehlagi[,c(1:3,10,11,6)]
  #getting regions with null count AA
  main <- AAcount
  main <- na.omit(main)
  uwa <- data.frame(Chr=0,Start=0,End=0,Ref=0,min=0)
  uwa <- uwa[FALSE,]
  for (i in unique(main$Chr)) {
    data_i <- filter(main,Chr==i)
    auu <- setdiff(seq(0, max(data_i$Start, na.rm = FALSE), by = 10000),data_i$Start)
    aau <- auu+10000
    khatam <- data.frame(Chr=i, Start=auu, End=aau, Ref=0,min=0)
    #khatam <- uwa
    uwa <- rbind(uwa,khatam)
  }
  #Ratio of windows with no count vs. ancestral allele 
  Simpul_Main <- AAcount %>% group_by(Chr) %>% summarise(Jumlah=n())
  Simpul_Main1 <- uwa %>% group_by(Chr) %>% summarise(Jumlah=n())
  Simpulan <- cbind.data.frame(Simpul_Main$Chr,Simpul_Main$Jumlah,Simpul_Main1$Jumlah)
  Simpulan$ratio <- Simpulan$`Simpul_Main1$Jumlah`/(Simpulan$`Simpul_Main$Jumlah`+Simpulan$`Simpul_Main1$Jumlah`)
  colnames(Simpulan) <- c("Chr","Windows_wt_AA","Windows_wo_AA","ratio")
  return(list(Above=output,Zero=uwa,Ratio=Simpulan))
}
