annv.ances<-function(AAcount, treshold, bar=bar1){
  #getting regions higher than bar1
  treshold$Chr <- as.factor(treshold$Chr)
  cobalagi<- left_join(AAcount,treshold,by="Chr")
  cobadehlagi<- subset.data.frame(cobalagi,subset= Ancestral_count>bar1)
  cobadehlagi$Ref <- 0 
  cobadehlagi$min <- "-"
  output <- cobadehlagi[,c(1:3,10,11,6)]
  #getting regions with null count AA
  auu <- setdiff(seq(0, max(AAcount$Start, na.rm = T), by = 10000),AAcount$Start)
  aau <- auu+10000
  uwa <- data.frame(Chr=as.factor((AAcount[1,1])), Start=auu, End=aau)
  uwa$Ref <- 0 
  uwa$min <- "-"
  #Ratio of windows with no count vs. ancestral allele 
  Simpul_Main <- AAcount %>% group_by(Chr) %>% summarise(Jumlah=n())
  Simpul_Main1 <- uwa %>% group_by(Chr) %>% summarise(Jumlah=n())
  Simpulan <- cbind.data.frame(Simpul_Main$Chr,Simpul_Main$Jumlah,Simpul_Main1$Jumlah)
  Simpulan$ratio <- Simpulan$`Simpul_Main1$Jumlah`/(Simpulan$`Simpul_Main$Jumlah`+Simpulan$`Simpul_Main1$Jumlah`)
  colnames(Simpulan) <- c("Chr","Windows_wt_AA","Windows_wo_AA","ratio")
  return(list(Above=output,Zero=uwa,Ratio=Simpulan))
}
