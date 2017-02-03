
# Authors: 
#         Alfredo Sanchez-Tojar, MPIO (Seewiesen) and ICL (Silwood Park), alfredo.tojar@gmail.com
#         Damien Farine, MPIO (Radolfzell) and EGI (Oxford), dfarine@orn.mpg.de

# Script first created on the 27th of October, 2016

###############################################################################
# Description of script and Instructions
###############################################################################

# This script aims to generate a figure that shows how the latent hierarchy
# affects the performance of the hierarchy method, Figure 2 from: 

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


###############################################################################
# Simulation: Estimating original Elo-rating and correlating it with real rank
###############################################################################

ptm <- proc.time()


avalues <- seq(0,30,5)
bvalues <- seq(-5,20,5)
N.inds.values <- c(50)
N.obs.values <- c(1,4,7,10,15,20,30,40,50)
poiss <- c(FALSE,FALSE,TRUE,TRUE)
dombias <- c(FALSE,TRUE,FALSE,TRUE)


#creating empty database
db <- data.frame(Ninds=integer(),
                 Nobs=integer(),
                 poiss=logical(),
                 dombias=logical(),
                 alevel=integer(),
                 blevel=integer(),
                 spearman=numeric(),
                 stringsAsFactors=FALSE)


# for each b value and each a value, estimate elo.rating for each number
# of interactions per individual.

for (typ in 1:length(poiss)){
  
  for (i in 1:length(bvalues)){
    
    for (j in 1:length(avalues)){
      
      for (p in 1:length(N.inds.values)){
        
        for (o in 1:length(N.obs.values)){
          
          for (simnum in 1:100){
            
            output <- generate_interactions(N.inds.values[p],
                                            N.inds.values[p]*N.obs.values[o],
                                            a=avalues[j],
                                            b=bvalues[i],
                                            id.biased=poiss[typ],
                                            rank.biased=dombias[typ])
            
            winner <- output$interactions$Winner
            loser <- output$interactions$Loser
            hierarchy <- output$hierarchy
            
            result.no.rand <- elo_scores(winner,
                                         loser,
                                         identities=c(1:N.inds.values[p]),
                                         init.score=1000,
                                         randomise=FALSE)
            
            spearman.cor<-cor(output$hierarchy$Rank,
                              rank(-result.no.rand,na.last="keep"),
                              use="complete.obs",method="spearman")
            
            
            #adding values to db
            db<-rbind(db,c(N.inds.values[p],N.obs.values[o],
                           poiss[typ],dombias[typ],
                           avalues[j],bvalues[i],
                           spearman.cor))
            
          }        
        }
      }
    }
  }
  
  #renaming variables in database
  names(db) <- c("Ninds","Nobs",
                 "poiss","dombias",
                 "alevel","blevel",
                 "spearman")
  
}


proc.time() - ptm


write.csv(db,
          "databases_package/Fig2_elo_no_rand_parameter_space_100sim_fixed_biases.csv",row.names=FALSE)


###############################################################################
# Plotting: MAIN TEXT: Elo-rating ~ real rank: spearman correlation - 95%CI
###############################################################################

eloparameterspace <- read.table("databases_package/Fig2_elo_no_rand_parameter_space_100sim_fixed_biases.csv",
                                header=TRUE,sep=",")


db<-eloparameterspace[eloparameterspace$poiss==1 & eloparameterspace$dombias==0,]


avalues <- seq(0,30,5)
bvalues <- c(-5,-5,-5,5,5,5,15,15,15)
N.inds.values <- c(50)
N.obs.values <- c(1,4,7,10,15,20,30,40,50)


a <- c("(a)","x","x","(b)","x","x","(c)","x","x")


tiff("plots/Figure2_Exploring_eloscores_and_steepness_poisson.tiff", 
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
  
  for (i in 1:length(bvalues)){
    
    db.3 <- db.2[db.2$blevel==bvalues[i],]
    
    if(i %in% c(1,4,7)){
      
      for (alev in 1:length(avalues)){
        
        assign(paste0("db.3.",alev),db.3[db.3$alevel==avalues[alev],])
        
        assign(paste0("db.4.",alev),
               summaryBy(spearman ~ Nobs, 
                         data = assign(paste0("db.3.",alev),
                                       db.3[db.3$alevel==avalues[alev],]), 
                         FUN = function(x) { c(m = mean(x),
                                               q = quantile(x,probs=c(0.025,0.975))) }))
      }
      
      names(db.4.1)<-c("Nobs","spearman","lower","upper")
      names(db.4.2)<-c("Nobs","spearman","lower","upper")
      names(db.4.3)<-c("Nobs","spearman","lower","upper")
      names(db.4.4)<-c("Nobs","spearman","lower","upper")
      names(db.4.5)<-c("Nobs","spearman","lower","upper")
      names(db.4.6)<-c("Nobs","spearman","lower","upper")
      names(db.4.7)<-c("Nobs","spearman","lower","upper")
      
      plot(db.3$spearman~db.3$Nobs,0.5,type="n",
           ylab="",
           xlab="",
           xaxt="n",
           yaxt="n",
           ylim=c(-0.6,1))
      
      if(i<7){
        
        axis(1,at=N.obs.values,
             cex.axis=1,tck=0.015,
             labels=FALSE)
        
        
      } else {
        
        axis(1,at=N.obs.values,
             labels=as.character(N.obs.values),
             cex.axis=1,tck=0.015)
        
        mtext("number of interactions/individual",
              side=1, adj=0, line=4, cex=1.8); 
        
      }
      
      axis(2,at=round(seq(-0.6,1,0.2),1),cex.axis=1.2,las=2)
      
      
      #adding points for the means and shadowed areas for the 95% CI
      colours <- brewer.pal(7,"Set1")
      
      points(db.4.1$Nobs-0.18,db.4.1$spearman,type="b",col=colours[1],pch=19)
      polygon(c(db.4.1$Nobs,rev(db.4.1$Nobs)),
              c(db.4.1$lower,rev(db.4.1$upper)),
              border=NA,col=rgb(250/255,128/255,114/255,0.15))
      
      points(db.4.2$Nobs-0.12,db.4.2$spearman,type="b",col=colours[2],pch=19)
      polygon(c(db.4.2$Nobs,rev(db.4.2$Nobs)),
              c(db.4.2$lower,rev(db.4.2$upper)),
              border=NA,col=rgb(141/255,182/255,205/255,0.15))
      
      points(db.4.3$Nobs-0.06,db.4.3$spearman,type="b",col=colours[3],pch=19)
      polygon(c(db.4.3$Nobs,rev(db.4.3$Nobs)),
              c(db.4.3$lower,rev(db.4.3$upper)),
              border=NA,col=rgb(188/255,238/255,104/255,0.15))
      
      points(db.4.4$Nobs,db.4.4$spearman,type="b",col=colours[4],pch=19)
      polygon(c(db.4.4$Nobs,rev(db.4.4$Nobs)),
              c(db.4.4$lower,rev(db.4.4$upper)),
              border=NA,col=rgb(154/255,50/255,205/255,0.15))
      
      points(db.4.5$Nobs+0.06,db.4.5$spearman,type="b",col=colours[5],pch=19)
      polygon(c(db.4.5$Nobs,rev(db.4.5$Nobs)),
              c(db.4.5$lower,rev(db.4.5$upper)),
              border=NA,col=rgb(244/255,164/255,96/255,0.15))
      
      points(db.4.6$Nobs+0.12,db.4.6$spearman,type="b",col=colours[6],pch=19)
      polygon(c(db.4.6$Nobs,rev(db.4.6$Nobs)),
              c(db.4.6$lower,rev(db.4.6$upper)),
              border=NA,col=rgb(255/255,228/255,181/255,0.15))
      
      points(db.4.7$Nobs+0.18,db.4.7$spearman,type="b",col=colours[7],pch=19)
      polygon(c(db.4.7$Nobs,rev(db.4.7$Nobs)),
              c(db.4.7$lower,rev(db.4.7$upper)),
              border=NA,col=rgb(139/255,69/255,19/255,0.15))
      
      btext <- paste("\nb = ",bvalues[i])
      text(48,-0.55,a[i],adj = 0,cex=1.5)
      text(30,-0.40,btext,adj = 0,cex=2)
      
            
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
        
        for (j in 1:length(avalues)){
          coline <- colours[j]
          plot_winner_prob(1:50,a=avalues[j],b=bvalues[i],coline)
          
        }
        
        mtext("P (dominant wins)",
              side=2, adj=0, line=3, cex=1.10); 
        
      } else {
        
        plot(c(1,50),c(0,1),type="n",
             ylab="", 
             xlab="",
             xaxt="n",
             yaxt="n",
             frame.plot=FALSE)    
        
        
        par(xpd=TRUE)
        legend(0,0.75,
               c("a = 0","a = 5","a = 10","a = 15"),
               col=colours[1:4],
               cex=1.35,bty='n',
               pch=rep(19,4),
               inset=c(0,0))
        
        legend(24,0.75,
               c("a = 20","a = 25","a = 30"),
               col=colours[5:7],
               cex=1.35,bty='n',
               pch=rep(19,4),
               inset=c(0,0))
        
        
        mtext("Difference in rank    ",
              side=3, adj=1, line=-2, cex=1.15); 
        
      }
      
    }
    
  }
  
  
  title(ylab = "spearman correlation coefficient",
        outer = TRUE, line = 0)
  
  par(mfrow=c(1,1))
  
}

dev.off()


###############################################################################
# Plotting: SUPPLEMENTARY MATERIAL: uniform, dominant bias, poisson*dominant bias
###############################################################################

# eloparameterspace <- read.table("databases_package/Fig2_elo_no_rand_parameter_space_100sim_fixed_biases.csv",
#                                 header=TRUE,sep=",")


db<-eloparameterspace[eloparameterspace$poiss==0 & eloparameterspace$dombias==0,]
#db<-eloparameterspace[eloparameterspace$poiss==0 & eloparameterspace$dombias==1,]
#db<-eloparameterspace[eloparameterspace$poiss==1 & eloparameterspace$dombias==1,]


avalues <- seq(0,30,5)
bvalues <- c(-5,-5,-5,5,5,5,15,15,15)
N.inds.values <- c(50)
N.obs.values <- c(1,4,7,10,15,20,30,40,50)


a <- c("(a)","x","x","(b)","x","x","(c)","x","x")

tiff("plots/supplements/FigureS_Exploring_eloscores_and_steepness_uniform.tiff", 
     #"plots/supplements/FigureS_Exploring_eloscores_and_steepness_dombias.tiff",
     #"plots/supplements/FigureS_Exploring_eloscores_and_steepness_poiss+dombias.tiff",
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
  
  for (i in 1:length(bvalues)){
    
    db.3 <- db.2[db.2$blevel==bvalues[i],]
    
    if(i %in% c(1,4,7)){
      
      for (alev in 1:length(avalues)){
        
        assign(paste0("db.3.",alev),db.3[db.3$alevel==avalues[alev],])
        
        assign(paste0("db.4.",alev),
               summaryBy(spearman ~ Nobs, 
                         data = assign(paste0("db.3.",alev),
                                       db.3[db.3$alevel==avalues[alev],]), 
                         FUN = function(x) { c(m = mean(x),
                                               q = quantile(x,probs=c(0.025,0.975))) }))
      }
      
      names(db.4.1)<-c("Nobs","spearman","lower","upper")
      names(db.4.2)<-c("Nobs","spearman","lower","upper")
      names(db.4.3)<-c("Nobs","spearman","lower","upper")
      names(db.4.4)<-c("Nobs","spearman","lower","upper")
      names(db.4.5)<-c("Nobs","spearman","lower","upper")
      names(db.4.6)<-c("Nobs","spearman","lower","upper")
      names(db.4.7)<-c("Nobs","spearman","lower","upper")
      
      plot(db.3$spearman~db.3$Nobs,0.5,type="n",
           ylab="",
           xlab="",
           xaxt="n",
           yaxt="n",
           ylim=c(-0.6,1))
      
      if(i<7){
        
        axis(1,at=N.obs.values,
             cex.axis=1,tck=0.015,
             labels=FALSE)
        
        
      } else {
        
        axis(1,at=N.obs.values,
             labels=as.character(N.obs.values),
             cex.axis=1,tck=0.015)
        
        mtext("number of interactions/individual",
              side=1, adj=0, line=4, cex=1.8); 
        
      }
      
      axis(2,at=round(seq(-0.6,1,0.2),1),cex.axis=1.2,las=2)
      
      
      #adding points for the means and shadowed areas for the 95% CI
      colours <- brewer.pal(7,"Set1")
      
      points(db.4.1$Nobs-0.18,db.4.1$spearman,type="b",col=colours[1],pch=19)
      polygon(c(db.4.1$Nobs,rev(db.4.1$Nobs)),
              c(db.4.1$lower,rev(db.4.1$upper)),
              border=NA,col=rgb(250/255,128/255,114/255,0.15))
      
      points(db.4.2$Nobs-0.12,db.4.2$spearman,type="b",col=colours[2],pch=19)
      polygon(c(db.4.2$Nobs,rev(db.4.2$Nobs)),
              c(db.4.2$lower,rev(db.4.2$upper)),
              border=NA,col=rgb(141/255,182/255,205/255,0.15))
      
      points(db.4.3$Nobs-0.06,db.4.3$spearman,type="b",col=colours[3],pch=19)
      polygon(c(db.4.3$Nobs,rev(db.4.3$Nobs)),
              c(db.4.3$lower,rev(db.4.3$upper)),
              border=NA,col=rgb(188/255,238/255,104/255,0.15))
      
      points(db.4.4$Nobs,db.4.4$spearman,type="b",col=colours[4],pch=19)
      polygon(c(db.4.4$Nobs,rev(db.4.4$Nobs)),
              c(db.4.4$lower,rev(db.4.4$upper)),
              border=NA,col=rgb(154/255,50/255,205/255,0.15))
      
      points(db.4.5$Nobs+0.06,db.4.5$spearman,type="b",col=colours[5],pch=19)
      polygon(c(db.4.5$Nobs,rev(db.4.5$Nobs)),
              c(db.4.5$lower,rev(db.4.5$upper)),
              border=NA,col=rgb(244/255,164/255,96/255,0.15))
      
      points(db.4.6$Nobs+0.12,db.4.6$spearman,type="b",col=colours[6],pch=19)
      polygon(c(db.4.6$Nobs,rev(db.4.6$Nobs)),
              c(db.4.6$lower,rev(db.4.6$upper)),
              border=NA,col=rgb(255/255,228/255,181/255,0.15))
      
      points(db.4.7$Nobs+0.18,db.4.7$spearman,type="b",col=colours[7],pch=19)
      polygon(c(db.4.7$Nobs,rev(db.4.7$Nobs)),
              c(db.4.7$lower,rev(db.4.7$upper)),
              border=NA,col=rgb(139/255,69/255,19/255,0.15))
      
      btext <- paste("\nb = ",bvalues[i])
      text(48,-0.55,a[i],adj = 0,cex=1.5)
      text(30,-0.40,btext,adj = 0,cex=2)
      
      
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
        
        for (j in 1:length(avalues)){
          coline <- colours[j]
          plot_winner_prob(1:50,a=avalues[j],b=bvalues[i],coline)
          
        }
        
        mtext("P (dominant wins)",
              side=2, adj=0, line=3, cex=1.10); 
        
      } else {
        
        plot(c(1,50),c(0,1),type="n",
             ylab="", 
             xlab="",
             xaxt="n",
             yaxt="n",
             frame.plot=FALSE)    
        
        
        par(xpd=TRUE)
        legend(0,0.75,
               c("a = 0","a = 5","a = 10","a = 15"),
               col=colours[1:4],
               cex=1.35,bty='n',
               pch=rep(19,4),
               inset=c(0,0))
        
        legend(24,0.75,
               c("a = 20","a = 25","a = 30"),
               col=colours[5:7],
               cex=1.35,bty='n',
               pch=rep(19,4),
               inset=c(0,0))
        
        
        mtext("Difference in rank    ",
              side=3, adj=1, line=-2, cex=1.15); 
        
      }
      
    }
    
  }
  
  
  title(ylab = "spearman correlation coefficient",
        outer = TRUE, line = 0)
  
  par(mfrow=c(1,1))
  
}

dev.off()
