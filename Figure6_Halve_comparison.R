
# Authors: 
#         Alfredo Sanchez-Tojar, MPIO (Seewiesen) and ICL (Silwood Park), alfredo.tojar@gmail.com
#         Damien Farine, MPIO (Radolfzell) and EGI (Oxford), dfarine@orn.mpg.de

# Script first created on the 27th of October, 2016

###############################################################################
# Description of script and Instructions
###############################################################################

# This script aims to generate a figure that shows how halve comparison
# changes depending on the steepness of the hierarchy and the number of
# interactions recorded. Figure 5 for:

# Sanchez-Tojar, A., Schroeder, J., Farine, D.R. (In preparation) A practical 
# guide for inferring reliable dominance hierarchies and estimating their 
# uncertainty

# more info at the Open Science Framework: https://osf.io/9gyek/


###############################################################################
# Packages needed
###############################################################################

# packages needed for this script

library(aniDom)
library(doBy)
library(RColorBrewer)


# Clear memory and get to know where you are
rm(list=ls())


###############################################################################
# Functions needed
###############################################################################

plot_winner_prob <- function(diff.rank, a, b,coline) {
  
  diff.rank.norm <- diff.rank/max(diff.rank)
  
  lines(diff.rank, 0.5+0.5/(1+exp(-diff.rank.norm*a+b)),col=coline)
  
}


estimate_uncertainty_by_splitting_2 <- 
  
  
  function (winners, losers, identities = NULL, sigmoid.param = 1/100, 
            K = 200, init.score = 0, randomise = FALSE, n.rands = 1000) 
  {
    if (is.null(identities)) {
      identities <- unique(c(winners, losers))
    }
    n.inds <- length(identities)
    if (randomise == FALSE) {
      n.rands <- 1
    }
    n.observations <- length(winners)
    #     obs1 <- seq(1,n.observations,2)
    #     obs2 <- seq(2,n.observations,2)
    obs1 <- seq(1,n.observations/2,1)
    obs2 <- seq(max(obs1)+1,n.observations,1)
    if (randomise == FALSE) {
      scores1 <- elo_scores(winners[obs1], losers[obs1], identities, 
                            sigmoid.param, K, init.score, randomise, n.rands)
      scores2 <- elo_scores(winners[obs2], losers[obs2], identities, 
                            sigmoid.param, K, init.score, randomise, n.rands)
      scores.cor <- cor(scores1, scores2, use="complete.obs", method = "spearman")
    }
    else {
      scores.cor <- rep(NA, n.rands)
      for (i in 1:n.rands) {
        new.order <- sample(1:n.observations)
        winners.tmp <- winners[new.order]
        losers.tmp <- losers[new.order]
        scores1 <- elo_scores(winners.tmp[obs1], losers.tmp[obs1], 
                              identities, sigmoid.param, K, init.score, randomise = FALSE)
        scores2 <- elo_scores(winners.tmp[obs2], losers.tmp[obs2], 
                              identities, sigmoid.param, K, init.score, randomise = FALSE)
        scores.cor[i] <- cor(scores1, scores2, use="complete.obs", method = "spearman")
      }
      CIs <- quantile(scores.cor, c(0.025, 0.975), na.rm = TRUE)
      scores.cor <- c(mean(scores.cor, na.rm = T), CIs[1], 
                      CIs[2])
      names(scores.cor) <- c("Mean", names(CIs))
    }
    invisible(scores.cor)
  }


###############################################################################
# Spliting and comparing halves of the data
###############################################################################

ptm <- proc.time()


db.split <- data.frame(Ninds=integer(),
                       Nobs=integer(),
                       poiss=logical(),
                       dombias=logical(),
                       alevel=integer(),
                       blevel=integer(),
                       elo.rand.split=numeric(),
                       lower=numeric(),
                       upper=numeric(),
                       stringsAsFactors=FALSE)


avalues <- c(10,15,30,15,10,5,0)
bvalues <- c(-5,0,5,5,5,5,5)
#N.inds.values <- c(50)
N.inds.values <- c(10)
N.obs.values <- c(1,4,7,10,15,20,30,40,50)
poiss <- c(FALSE,FALSE,TRUE,TRUE)
dombias <- c(FALSE,TRUE,FALSE,TRUE)


for (typ in 1:length(poiss)){
  
  for (j in 1:length(avalues)){
    
    for (p in 1:length(N.inds.values)){
      
      for (o in 1:length(N.obs.values)){
        
        for (numsim in 1:100){
          
          output <- generate_interactions(N.inds.values[p],
                                          N.inds.values[p]*N.obs.values[o],
                                          a=avalues[j],
                                          b=bvalues[j],
                                          id.biased=poiss[typ],
                                          rank.biased=dombias[typ])
          
          winner <- output$interactions$Winner
          loser <- output$interactions$Loser
          
          
          split.res <- estimate_uncertainty_by_splitting_2(winner,
                                                           loser, 
                                                           identities=c(1:N.inds.values[p]),
                                                           init.score=1000,
                                                           randomise = TRUE, 
                                                           n.rands = 1000)
          
          db.split<-rbind(db.split,c(N.inds.values[p],N.obs.values[o],
                                     poiss[typ],
                                     dombias[typ],
                                     avalues[j],
                                     bvalues[j],
                                     split.res[[1]],
                                     split.res[[2]],
                                     split.res[[3]]))
          
          
        }
      }
    }
  }
}

names(db.split) <- c("Ninds","Nobs",
                     "poiss","dombias",
                     "alevel","blevel",
                     "elo.rand.split",
                     "lower","upper")


proc.time() - ptm


# write.csv(db.split,
#           "databases_package/Fig6_db_split_elorand_100sim_fixed_biases.csv",row.names=FALSE)

write.csv(db.split,
          "databases_package/final_data_for_Figures_backup/Fig6_db_split_elorand_100sim_fixed_biases_10ind.csv",row.names=FALSE)



###############################################################################
# Plotting estimated rank/2 ~ estimated rank/2: spearman correaltion and
# 95% CI intervals 
###############################################################################

db_split100sim0 <- 
  read.table("databases_package/final_data_for_Figures_backup/Fig6_db_split_elorand_100sim_fixed_biases_full.csv",header=TRUE,sep=",")


#db<-db_split100sim0[db_split100sim0$poiss==1 & db_split100sim0$dombias==0,]
#db<-db_split100sim0[db_split100sim0$poiss==0 & db_split100sim0$dombias==0,]
#db<-db_split100sim0[db_split100sim0$poiss==0 & db_split100sim0$dombias==1,]
db<-db_split100sim0[db_split100sim0$poiss==1 & db_split100sim0$dombias==1,]


avalues <- c(15,15,15,10,10,10,5,5,5)
bvalues <- c(5,5,5,5,5,5,5,5,5)
N.inds.values <- c(50)
#N.inds.values <- c(10)
N.obs.values <- c(1,4,7,10,15,20,30,40,50)


a <- c("(a)","x","x","(b)","x","x","(c)","x","x")


tiff(#"plots/Figure6_Halve_comparison_poisson.tiff",
     #"plots/supplements/FigureS6_Halve_comparison_uniform.tiff", 
     #"plots/supplements/FigureS13_Halve_comparison_dombias.tiff",
     "plots/supplements/FigureS20_Halve_comparison_poiss+dombias.tiff",
     height=29.7, width=21,
     units='cm', compression="lzw", res=600)


for (p in 1:length(N.inds.values)){
  
  m <- rbind(c(1,1,2),c(1,1,3),
             c(4,4,5),c(4,4,6),
             c(7,7,8),c(7,7,9))
  
  layout(m)
  
  op <- par(oma = c(6,3,1,1) + 0.1,
            mar = c(0.5,5,1,0) + 0.1,
            cex.lab=2.5)
  
  db.2 <- db[db$Ninds==N.inds.values[p],]
  
  for (i in 1:length(avalues)){
    
    if(i %in% c(1,4,7)){
      
      db.3 <- db.2[db.2$alevel==avalues[i]
                   & db.2$blevel==bvalues[i]
                   ,]
      
      db.4 <-summaryBy(elo.rand.split ~ Nobs, 
                       data = db.3, 
                       FUN = function(x) { c(m = mean(x),
                                             q = quantile(x,probs=c(0.025,0.975))) })
      
      names(db.4) <- c("Nobs",
                       "elo.rand.split.m","lower","upper")
      
      plot(db.4$elo.rand.split.m~db.4$Nobs,0.5,type="n",
           ylab="",
           xlab="",
           xaxt="n",
           yaxt="n",
           ylim=c(-0.4,1))
      
      if(i<7){
        
        axis(1,at=N.obs.values,
             cex.axis=1,tck=0.015,
             labels=FALSE)
        
        
      } else {
        
        axis(1,at=N.obs.values,
             labels=as.character(N.obs.values),
             cex.axis=1,tck=0.015)
        
        mtext(" ratio of interactions to individuals",
              side=1, adj=0, line=4, cex=1.8); 
        
      }
      
      axis(2,at=seq(-0.4,1,0.2),cex.axis=1.2,las=2)
      
      
      #adding points for the means and shadowed areas for the 95% CI
      points(db.4$Nobs,db.4$elo.rand.split.m,type="b",col="blue",pch=19)
      polygon(c(db.4$Nobs,rev(db.4$Nobs)),
              c(db.4$lower,rev(db.4$upper)),
              border=NA,col=rgb(0,0,1, 0.15))
      
      lines(c(0,51),c(0.35,0.35),col="red",lty=3,lwd=1.5)
      
    
      text(48,-0.37,a[i],adj = 0 ,cex=1.5)

    }else {
      
      if(i %in% c(2,5,8)){
        
        plot(c(1,50),c(0.5,1),type="n",
             ylab="", 
             xlab="",
             xaxt="n",
             yaxt="n",
             cex=1.5)
        
        axis(1,at=seq(0,50,10),
             cex.axis=1,tck=0.015)
        
        axis(2,at=seq(0.5,1,0.1),cex.axis=1.2,las=2,tck=0.015) 
        
        plot_winner_prob(1:50,a=avalues[i],b=bvalues[i],"black")
        
        mtext("P (dominant wins)",
              side=2, adj=0, line=3, cex=1.10); 
        
      } else {
        
        plot(c(1,50),c(0,1),type="n",
             ylab="", 
             xlab="",
             xaxt="n",
             yaxt="n",
             frame.plot=FALSE)    
        
        
        atext <- paste("\na = ",avalues[i])
        btext <- paste("\nb = ",bvalues[i])
        ttext <- paste0(atext,btext,sep="\n")
        text(15,0.45,ttext,adj = 0,cex=2.25)
        
        
        mtext("Difference in rank    ",
              side=3, adj=1, line=-2, cex=1.15); 
        
      }
      
    }
    
  }
  
  
  title(ylab = "Spearman rank correlation coefficient",
        outer = TRUE, line = 0)
  
  par(mfrow=c(1,1))
  
}

dev.off()


###############################################################################
# Plotting: SUPPLEMENTARY MATERIAL: steep and flat scenarios
###############################################################################

db_split100sim0 <- 
  read.table("databases_package/final_data_for_Figures_backup/Fig6_db_split_elorand_100sim_fixed_biases_full.csv",header=TRUE,sep=",")


#db<-db_split100sim0[db_split100sim0$poiss==1 & db_split100sim0$dombias==0,]
#db<-db_split100sim0[db_split100sim0$poiss==0 & db_split100sim0$dombias==0,]
#db<-db_split100sim0[db_split100sim0$poiss==0 & db_split100sim0$dombias==1,]
db<-db_split100sim0[db_split100sim0$poiss==1 & db_split100sim0$dombias==1,]


avalues <- c(10,10,10,15,15,15,30,30,30,0,0,0)
bvalues <- c(-5,-5,-5,0,0,0,5,5,5,5,5,5)
N.inds.values <- c(50)
#N.inds.values <- c(10)
N.obs.values <- c(1,4,7,10,15,20,30,40,50)


a <- c("(a)","x","x","(b)","x","x","(c)","x","x","(d)","x","x")


tiff(#"plots/supplements/FigureS3_Halve_comparison_steep_and_flat_poisson.tiff",
  #"plots/supplements/FigureS7_Halve_comparison_steep_and_flat_uniform.tiff", 
  #"plots/supplements/FigureS14_Halve_comparison_steep_and_flat_dombias.tiff",
  "plots/supplements/FigureS21_Halve_comparison_steep_and_flat_poiss+dombias.tiff",
  #"plots/supplements/FigureS13_Halve_comparison_steep_and_flat_poisson_10ind.tiff",
  #"plots/supplements/FigureS15_Halve_comparison_steep_and_flat_uniform_10ind.tiff",
  height=29.7, width=21,
  units='cm', compression="lzw", res=600)


for (p in 1:length(N.inds.values)){
  
  m <- rbind(c(1,1,2),c(1,1,3),
             c(4,4,5),c(4,4,6),
             c(7,7,8),c(7,7,9),
             c(10,10,11),c(10,10,12))
  
  layout(m)
  
  op <- par(oma = c(6,3,1,1) + 0.1,
            mar = c(0.5,5,1,0) + 0.1,
            cex.lab=2.5)
  
  db.2 <- db[db$Ninds==N.inds.values[p],]
  
  for (i in 1:length(avalues)){
    
    if(i %in% c(1,4,7,10)){
      
      db.3 <- db.2[db.2$alevel==avalues[i]
                   & db.2$blevel==bvalues[i]
                   ,]
      
      db.4 <-summaryBy(elo.rand.split ~ Nobs, 
                       data = db.3, 
                       FUN = function(x) { c(m = mean(x),
                                             q = quantile(x,probs=c(0.025,0.975))) })
      
      names(db.4) <- c("Nobs",
                       "elo.rand.split.m","lower","upper")
      
      plot(db.4$elo.rand.split.m~db.4$Nobs,0.5,type="n",
           ylab="",
           xlab="",
           xaxt="n",
           yaxt="n",
           ylim=c(-0.4,1))
      
      if(i<10){
        
        axis(1,at=N.obs.values,
             cex.axis=1,tck=0.015,
             labels=FALSE)
        
        
      } else {
        
        axis(1,at=N.obs.values,
             labels=as.character(N.obs.values),
             cex.axis=1,tck=0.015)
        
        mtext(" ratio of interactions to individuals",
              side=1, adj=0, line=4, cex=1.8); 
        
      }
      
      axis(2,
           at=seq(-0.4,1,0.2),
           #at=seq(0,1,0.1),
           cex.axis=1.2,las=2)
      
      
      #adding points for the means and shadowed areas for the 95% CI
      points(db.4$Nobs,db.4$elo.rand.split.m,type="b",col="blue",pch=19)
      polygon(c(db.4$Nobs,rev(db.4$Nobs)),
              c(db.4$lower,rev(db.4$upper)),
              border=NA,col=rgb(0,0,1, 0.15))
      
      lines(c(0,51),c(0.35,0.35),col="red",lty=3,lwd=1.5)
      
      
      text(48,-0.35,a[i],adj = 0 ,cex=1.5)
      #text(48,0.03,a[i],adj = 0 ,cex=1.5)
      
      
    }else {
      
      if(i %in% c(2,5,8,11)){
        
        plot(c(1,50),
             #c(1,10),
             c(0.5,1),type="n",
             ylab="", 
             xlab="",
             xaxt="n",
             yaxt="n",
             cex=1.5)
        
        axis(1,
             at=seq(0,50,10),
             #at=seq(0,10,1),
             cex.axis=1,tck=0.015)
        
        axis(2,at=seq(0.5,1,0.1),cex.axis=1.2,las=2,tck=0.015) 
        
        plot_winner_prob(1:50,a=avalues[i],b=bvalues[i],"black")
        #plot_winner_prob(1:10,a=avalues[i],b=bvalues[i],"black")
        
        mtext("P (dominant wins)",
              side=2, adj=0, line=3, cex=0.85); 
        
      } else {
        
        plot(c(1,50),c(0,1),type="n",
             ylab="", 
             xlab="",
             xaxt="n",
             yaxt="n",
             frame.plot=FALSE)    
        
        
        atext <- paste("\na = ",avalues[i])
        btext <- paste("\nb = ",bvalues[i])
        ttext <- paste0(atext,btext,sep="\n")
        text(15,0.35,ttext,adj = 0,cex=2.25)
        
        
        mtext("Difference in rank     ",
              side=3, adj=1, line=-2, cex=1); 
        
      }
      
    }
    
  }
  
  
  title(ylab = "Spearman rank correlation coefficient",
        outer = TRUE, line = 0)
  
  par(mfrow=c(1,1))
  
}

dev.off()