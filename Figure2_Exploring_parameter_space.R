
# Authors: 
#         Alfredo Sanchez-Tojar, MPIO (Seewiesen) and ICL (Silwood Park), alfredo.tojar@gmail.com
#         Damien Farine, MPIO (Radolfzell) and EGI (Oxford), dfarine@orn.mpg.de

# Script first created on the 27th of October, 2016

###############################################################################
# Description of script and Instructions
###############################################################################

# This script aims to generate a figure that shows the variety of hierarchy
# scenarios explored in Figure 2 from: 

# Sanchez-Tojar, A., Schroeder, J., Farine, D.R. (In preparation) A practical 
# guide for inferring reliable dominance hierarchies and estimating their 
# uncertainty

# more info at the Open Science Framework: http://doi.org/10.17605/OSF.IO/9GYEK


###############################################################################
# Packages needed
###############################################################################

# packages needed this script

library(RColorBrewer)


# Clear memory
rm(list=ls())


###############################################################################
# Functions needed
###############################################################################

plot_winner_prob <- function(diff.rank, a, b,coline) {
  
  diff.rank.norm <- diff.rank/max(diff.rank)
  
  lines(diff.rank, 0.5+0.5/(1+exp(-diff.rank.norm*a+b)),col=coline)
  
}


###############################################################################
# Multi-pannel plot to explore parameter space: minimal version = 3 pannels
###############################################################################

avalues <- seq(0,30,5)
bvalues <- seq(-5,5,5)


a1 <- c("(a)","(b)","(c)")
colours <- brewer.pal(7,"Set1")


tiff("plots/Figure2_Exploring_parameter_space_25ind_REDUCED.tiff",
  #"plots/supplements/after_revision/FigureS01_Exploring_parameter_space_10ind_REDUCED.tiff",
  #"plots/supplements/after_revision/FigureS08_Exploring_parameter_space_50ind_REDUCED.tiff",
  height=21, width=29.7,
  units='cm', compression="lzw", res=600)


op <- par(oma = c(6,6,1,1) + 0.1,
          mar = c(0.5,0.5,0,0) + 0.1,
          cex.lab=2.5)


par(mfrow=c(1,3))


# For each value of b plot all the curves generated by changing the value of a

for (i in 1:length(bvalues)){
  
  plot(#c(1,10),
       c(1,25),
       #c(1,50),
       c(0.3,1),type="n",
       ylab="",
       xlab="",
       xaxt="n",
       yaxt="n")
  
  for (j in 1:length(avalues)){
    coline <- colours[j]
    #plot_winner_prob(1:10,a=avalues[j],b=bvalues[i],coline)
    plot_winner_prob(1:25,a=avalues[j],b=bvalues[i],coline)
    #plot_winner_prob(1:50,a=avalues[j],b=bvalues[i],coline)
    
  }
  
  axis(1,
       #at=seq(1,10,1),
       at=seq(1,25,2),
       #at=seq(1,50,6),
       cex.axis=1,tck=0.015)
  
  if(i==1){
    
    axis(2,at=seq(0.3,1,0.1),cex.axis=1.2,las=2,tck=0.015)
    
  } else {
    
    axis(2,at=seq(0.3,1,0.1),cex.axis=1.2,las=2,tck=0.015,
         labels=FALSE)
    
  }
  
  #lines(c(0,12),c(0.5,0.5),col="red",lty=3,lwd=1.5) # line at randomness, i.e. P=0.5
  lines(c(0,27),c(0.5,0.5),col="red",lty=3,lwd=1.5) # line at randomness, i.e. P=0.5
  #lines(c(0,57),c(0.5,0.5),col="red",lty=3,lwd=1.5) # line at randomness, i.e. P=0.5
  
  
  btext <- paste("\nb = ",bvalues[i])
  # text(9.6,0.29,a1[i],adj = 0,cex=1.5)
  # text(6,0.42,btext,adj = 0,cex=2.5)
  text(24,0.29,a1[i],adj = 0,cex=1.5)
  text(15,0.42,btext,adj = 0,cex=2.5)
  # text(47,0.29,a1[i],adj = 0,cex=1.5)
  # text(30,0.42,btext,adj = 0,cex=2.5)
  
  par(xpd=TRUE)
  legend(4,
         0.47,
         c("a = 0","a = 5","a = 10","a = 15","a = 20","a = 25","a = 30"),
         col=colours[1:7],
         cex=1.5,bty='n',
         pch=rep(19,4),
         inset=c(0,0))
  
}

title(xlab = "Difference in rank",
      ylab = "Probability of higher rank winning",
      outer = TRUE, line = 3)


dev.off()


# ###############################################################################
# # Multi-pannel plot to explore parameter space (i.e. hierarchy scenarios): preprint version
# ###############################################################################
# 
# avalues <- seq(0,30,5)
# bvalues <- seq(-5,35,5)
# 
# 
# a1 <- c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)")
# colours <- brewer.pal(7,"Set1")
# 
# 
# tiff(#"plots/Figure2_Exploring_parameter_space.tiff", 
#      "plots/FigureS1_Exploring_parameter_space_10ind.tiff", 
#      height=29.7, width=21,
#      units='cm', compression="lzw", res=600)
# 
# 
# op <- par(mfrow = c(3,2),
#           oma = c(6,6,1,1) + 0.1,
#           mar = c(0.5,0.5,0,0) + 0.1,
#           cex.lab=2.5)
# 
# 
# par(mfrow=c(3,3))
# 
# 
# # For each value of b plot all the curves generated by changing the value of a
# 
# for (i in 1:length(bvalues)){
#   
#   plot(#c(1,50),
#        c(1,10),
#        c(0.2,1),type="n",
#        #ylab="P of higher rank winning", 
#        #xlab="Difference in rank",
#        ylab="",
#        xlab="",
#        xaxt="n",
#        yaxt="n")
#   
#   for (j in 1:length(avalues)){
#     coline <- colours[j]
#     #plot_winner_prob(1:50,a=avalues[j],b=bvalues[i],coline)
#     plot_winner_prob(1:10,a=avalues[j],b=bvalues[i],coline)
# 
#   }
#   
#   if(i<7){
#     
#     axis(1,
#          #at=seq(0,50,10),
#          at=seq(0,10,1),
#          cex.axis=1,tck=0.015,
#          labels=FALSE)
#     
#   } else {
#     
#     axis(1,
#          #at=seq(0,50,10),
#          at=seq(0,10,1),
#          #labels=as.character(N.obs.values),
#          cex.axis=1,tck=0.015)
#     
#   }
#   
#   
#   if(i==1 | i==4 | i==7){
#     
#     axis(2,at=seq(0.2,1,0.1),cex.axis=1.2,las=2,tck=0.015)
#     
#   } else {
#     
#     axis(2,at=seq(0.2,1,0.1),cex.axis=1.2,las=2,tck=0.015,
#          labels=FALSE)
#     
#   }
# 
#   #lines(c(0,52),c(0.5,0.5),col="red",lty=3,lwd=1.5) # line at randomness, i.e. P=0.5
#   lines(c(0,12),c(0.5,0.5),col="red",lty=3,lwd=1.5) # line at randomness, i.e. P=0.5
#   
#     
#   btext <- paste("\nb = ",bvalues[i])
#   #text(46,0.22,a1[i],adj = 0,cex=1.5)
#   #text(34,0.37,btext,adj = 0,cex=2)
#   text(9,0.22,a1[i],adj = 0,cex=1.5)
#   text(7,0.37,btext,adj = 0,cex=2)
#   
#   par(xpd=TRUE)
#   legend(#3,
#          1,
#          0.4,
#          c("a = 0","a = 5","a = 10","a = 15"),
#          col=colours[1:4],
#          cex=1,bty='n',
#          pch=rep(19,4),
#          inset=c(0,0))
#   
#   legend(#17,
#          3.5,
#          0.4,
#          c("a = 20","a = 25","a = 30"),
#          col=colours[5:7],
#          cex=1,bty='n',
#          pch=rep(19,4),
#          inset=c(0,0))
#   
#   
# }
# 
# title(xlab = "Difference in rank",
#       ylab = "Probability of higher rank winning",
#       outer = TRUE, line = 3)
# 
# dev.off()