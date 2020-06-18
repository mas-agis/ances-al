#Getting average ancestral count of window per chromosome
average <- function(AAcount,treshold){
  main <- read.table("taurusCount_null_changes_alleles.txt", header = T)
  main <- taurus$AACounts
  treshold <- read.table("taurustreshold_null_changes.txt", header=T)
  treshold <- taurus$treshold
  str(main)
  str(treshold)
  rata2 <- main %>%
    group_by(Chr) %>%
    summarise(Average = mean(Ancestral_count))
  return(rata2)
  }
