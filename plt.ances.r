plt.ances <- function(AAcounts,treshold,cr){
    akhir1 <- subset.data.frame(AAcounts,Chr==cr)
    treshol <- subset.data.frame(treshold,Chr==cr)
    akhir1$Window <- as.numeric(akhir1$Window)
    breaks <- seq(0, nrow(akhir1), 1200) 
    p <- ggplot(akhir1, aes(x = Window, y = Ancestral_count)) + geom_point(col="wheat4") + scale_x_continuous(breaks = breaks ) +
      geom_smooth(method = "loess", span = 0.9, col = "red") + 
      geom_smooth(method = "loess", span = 0.2, col = "blue")  +  
      labs(title = paste0("Ancestral count in chr ",cr), x = "Window", y= "Ancestral Count", subtitle = paste0("Mean ", round(mean(akhir1$Ancestral_count),2), "; SD ", round(sd(akhir1$Ancestral_count),2))) +
      geom_hline(aes(yintercept = treshol$bar1, color="Treshold 1"),linetype=2, size=1) +
      geom_hline(aes(yintercept = treshol$bar2, color="Treshold 2"),linetype=2, size=1) + 
      geom_hline(aes(yintercept = treshol$bar3, color="Treshold 3"),linetype=2, size=1) +
      scale_color_manual(name = "Quantile", values = c( "goldenrod1","orangered2", "magenta3"))
      last_plot()
  }
