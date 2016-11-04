
# Authors: 
#         Alfredo Sanchez-Tojar, MPIO (Seewiesen) and ICL (Silwood Park), alfredo.tojar@gmail.com
#         Damien Farine, MPIO (Radolfzell) and EGI (Oxford), dfarine@orn.mpg.de

# Script first created on the 27th of October, 2016

###############################################################################
# Description of script and Instructions
###############################################################################

# This script aims to test different methods for estimating individual
# dominance rank.


###############################################################################
# Packages needed
###############################################################################

# packages needed to be loaded for this script (a couple of them might be only needed in the following
# script)

library(EloRating)
library(EloChoice)
library(doBy)
library(plyr)


# Clear memory and get to know where you are
rm(list=ls())


###############################################################################
# Functions
###############################################################################


calculate_winner <- function(rank1, rank2, a, b) {
	
	diff.rank <- abs(rank1 - rank2)
	diff.rank.norm <- diff.rank/max(diff.rank)

	p.win <- 0.5+0.5/(1+exp(-diff.rank.norm*a+b))

	winner <- sample(c(1,2),1,prob=c(p.win,1-p.win))

	if (winner == 1) {
		if (rank1 < rank2) {
			winner.loser <- c(1,2)
		} else {
			winner.loser <- c(2,1)
		}
	} else {
		if (rank1 > rank2) {
			winner.loser <- c(1,2)
		} else {
			winner.loser <- c(2,1)
		}
	}
	return(winner.loser)

}



plot_winner_prob <- function(diff.rank, a, b,coline) {

	diff.rank.norm <- diff.rank/max(diff.rank)
  
	lines(diff.rank, 0.5+0.5/(1+exp(-diff.rank.norm*a+b)),col=coline)

	#plot(diff.rank, 0.5+0.5/(1+exp(-diff.rank.norm*a+b)),ylim=c(0,1), type='l',ylab="P of higher rank winning", xlab="Difference in rank")

}



select_interactants <- function(hierarchy) {

	interactants <- sample(hierarchy$ID,2)

	return(interactants)

}



generate_interactions <- function(N.inds, N.obs, a, b) {

	hierarchy <- data.frame(ID=1:N.inds, Rank=1:N.inds)

	interactions <- data.frame(Winner=rep(NA,N.obs), Loser=rep(NA,N.obs))

	for (i in 1:N.obs) {
	
		ints <- select_interactants(hierarchy)
		outcome <- calculate_winner(hierarchy$Rank[ints[1]],hierarchy$Rank[ints[2]],a,b)

		interactions$Winner[i] <- hierarchy$ID[ints[outcome[1]]]
		interactions$Loser[i] <- hierarchy$ID[ints[outcome[2]]]
	
	}
  
  interactions$Date <- "2016-10-31"
  interactions$Date[1:5] <- "2016-10-30"

	return(list(hierarchy=hierarchy,interactions=interactions))

}


elo.scores <- function(winners,losers,n.inds=NULL,sigmoid.param=1/100,
                       K=200,init.score=0,n.rands=1000,
                       return.trajectories=FALSE){
  
  if(is.null(n.inds)){
    n.inds <- max(c(unique(winners),unique(losers)))	
    #n.inds <- length(unique(c(winners,losers)))
  }
  
  T <- length(winners)
  
  if(return.trajectories){
    all.scores <- array(0,c(n.inds,T+1,n.rands))
  } else{
    all.scores <- array(0,c(n.inds,n.rands))
  }
  
  if (length(K) == 1) {
    K <- rep(K,T)
  }
  
  for(r in 1:n.rands){
    if (n.rands == 1) {
      ord <- 1:T
    } else {
      ord <- sample(1:T,T,replace=F)
    }
    winners.perm <- winners[ord]
    losers.perm <- losers[ord]
    scores<-array(NA,c(n.inds,T+1))
    scores[,1]<-init.score
    
    for(i in 1:T){
      
      scores[,i+1] <- scores[,i]
      
      winner <- winners.perm[i]
      loser <- losers.perm[i]
      p<-1/(1+exp(-sigmoid.param*(scores[winner,i]-scores[loser,i]))) #prob that winner wins
      
      if(scores[winner,i] >= scores[loser,i]){
        scores[winner,i+1] <- scores[winner,i] + (1-p)*K[i]
        scores[loser,i+1] <- scores[loser,i] - (1-p)*K[i]
      }
      else{
        scores[winner,i+1] <- scores[winner,i] + p*K[i]
        scores[loser,i+1] <- scores[loser,i] - p*K[i]
      }
    }
    
    if(return.trajectories){
      all.scores[,,r]<-scores
    } else{
      all.scores[,r]<-scores[,T+1]
    }
    
    
  }
  
  freq <- table(factor(c(winners,losers),levels=c(1:20)))
  all.scores[as.numeric(names(freq)[which(freq==0)]),] <- NA
  
  invisible(all.scores)	
}



###############################################################################
# SECTION 1: EXPLORTING PARAMETER SPACE AND ELO-RATING BEHAVIOUR WITHIN
###############################################################################

###############################################################################
# Multi-pannel plot to explore parameter space in generate_interactions
###############################################################################

# counter <- 95 # I use counter to create a grey scale
# avalues <- seq(0,36,2) #all values of a to be explored
# bvalues <- c(-5,0,5,10,15,20,25,30,35) #all values of b to be explored
counter <- 75
avalues <- seq(0,30,5)
bvalues <- c(-5,0,5,10,15,20,25,30,35)


par(mfrow=c(3,3)) #multi-pannel 3 by 3 = 9 plots at once!


# For each value of b plot all the curves generated by changing the value of a

for (i in 1:length(bvalues)){
  
  plot(c(1,50),c(0,1),type="n",
       ylab="P of higher rank winning", xlab="Difference in rank")
  btext <- paste("b = ",bvalues[i])
  text(10,0.35,btext)
  
  for (j in 1:length(avalues)){
    coline <- paste0("grey",counter)
    plot_winner_prob(1:50,a=avalues[j],b=bvalues[i],coline)
    #counter <- counter-5     
    counter <- counter-10 
  }
  #counter <- 95
  counter <- 75
  lines(c(1,50),c(0.5,0.5),col="red",lty=3) # line at randomness, i.e. P=0.5
#   legend("bottomright",c("a = 0","a = 18","a = 36"),#custom made
#          lty=c(1,1,1),col=c("grey95","grey50","grey0"),#custom made
#          cex=.5,bty='n',
#          y.intersp=0.2)
  legend("bottomleft",c("a = 0","a = 15","a = 30"),
         lty=c(1,1,1),col=c("grey75","grey45","grey15"),
         cex=.75,bty='n',
         x.intersp=0.1,
         y.intersp=0.2,
         inset=c(-0.2,-0.1),
         seg.len=0.2)
}

par(mfrow=c(1,1))


###############################################################################
# Estimating Elo-rating and correlating it with real rank
###############################################################################

library(EloRating)

ptm <- proc.time()

avalues <- seq(0,30,5)
bvalues <- seq(-5,20,5)
N.inds.values <- c(50)
N.obs.values <- c(1,3,5,7,9,
                  11,13,15,17,19,
                  21,23,25,27,29,
                  31,33,35,37,39,
                  41,43,45,47,49)
# N.obs.values <- c(1,2,3,4,5,6,7,8,9,10,
#                   11,12,13,14,15,16,17,18,19,20,
#                   30,40,50)

# #for quick checkings
# avalues<-5
# bvalues<-0
# N.inds.values <- c(50)
# N.obs.values <- c(30)


#creating empty database
db <- data.frame(Ninds=integer(),
                 Nobs=integer(),
                 alevel=integer(),
                 blevel=integer(),
                 spearman=numeric(),
                 stringsAsFactors=FALSE)


# for each b value and each a value, estimate elo.rating for each number
# of interactions per individual.

for (i in 1:length(bvalues)){
  
  for (j in 1:length(avalues)){
    
    for (p in 1:length(N.inds.values)){
      
      for (o in 1:length(N.obs.values)){
        
        output <- generate_interactions(N.inds.values[p],
                                        N.inds.values[p]*N.obs.values[o],
                                        a=avalues[j],
                                        b=bvalues[i])
        
        winner <- output$interactions$Winner
        loser <- output$interactions$Loser
        date <- output$interactions$Date
        hierarchy <- output$hierarchy
        
#         x<-elo.seq(winner=as.factor(winner),
#                    loser=as.factor(loser), 
#                    Date=date,
#                    progressbar=FALSE,
#                    k=200)
        
        result.no.rand <- elo.scores(winner,
                                     loser,
                                     n.rands=1, # estimated only once
                                     init.score=1000,
                                     n.inds=N.inds.values[p])

        spearman.cor<-cor(output$hierarchy$Rank,
                          rank(-result.no.rand),
                          use="complete.obs",method="spearman")
        
#         w<-extract.elo(x,standardize = FALSE) # extract elo from elo.seq
# 
#         
#         scores <- as.data.frame(w,
#                                 row.names = as.character(seq(1,length(w),
#                                                              1)))
#         #dataframe with Elo-rating and id
#         z <- cbind(attributes(w),
#                    scores)
#         
#         # since id equal rank I can just generate a new variable that is rank
#         z$rank <- as.numeric(as.character(z$names))
#         
#         #ranking Elo-rating
#         z$Elo.ranked <- rank(-z$w)
#         
#         spearman.cor<-cor(z$rank,z$Elo.ranked,use="complete.obs",method="spearman")
#         spearman.cor<-cor(result.no.rand,z$w,use="complete.obs",method="spearman")
        
        #adding values to db
        db<-rbind(db,c(N.inds.values[p],N.obs.values[o],
                       avalues[j],bvalues[i],
                       spearman.cor))
        
      }
    }
  }
}

#renaming variables in database
names(db) <- c("Ninds","Nobs","alevel","blevel","spearman")

proc.time() - ptm



###############################################################################
# Plotting Elo-rating ~ real rank: spearman correaltion
###############################################################################

for (p in 1:length(N.inds.values)){
  
  par(mfrow=c(3,2))
  
  db.2 <- db[db$Ninds==N.inds.values[p],]
  
  for (i in 1:length(bvalues)){
    
    db.3 <- db.2[db.2$blevel==bvalues[i],]
    
    palette(c("black","red","green3","blue","cyan","magenta","gray50"))
    
    plot(db.3$spearman~jitter(db.3$Nobs,0.6),col=as.factor(db.3$alevel),
         ylab="spearman correlation", xlab="number of interactions/individual",
         ylim=c(-0.6,1), pch=19)
    
    Nindtext <- paste("N.ind = ",N.inds.values[p])
    btext <- paste("\nb = ",bvalues[i])
    ttext <- paste0(Nindtext,btext,sep="\n")
    text(45,0.6,ttext)
    
    text(1,-0.55,"a:")
    
    par(xpd=TRUE)
    
    legend(-3.2,0.25,as.character(avalues),
           col=1:length(avalues),
           cex=0.95,bty='n',
           #y.intersp=0,
           x.intersp=0.1,
           horiz=TRUE,
           xjust=0,
           #yjust=0,
           inset=c(0,0),
           pch=rep(19,6),
           text.width=rep(0,7))
    
    lines(c(-1,52),c(-0.4,-0.4),col="black")
    
  }
  
  palette("default")
  par(mfrow=c(1,1))
  
}



###############################################################################
# SECTION 2: COMPARING THE 4 METHODS
###############################################################################

ptm <- proc.time()

db <- data.frame(Ninds=integer(),
                 Nobs=integer(),
                 alevel=integer(),
                 blevel=integer(),
                 Ndavid=numeric(),
                 elo.original=numeric(),
                 elo.no.rand=numeric(),
                 elo.rand=numeric(),
                 elochoice.no.rand=numeric(),
                 elochoice.rand=numeric(),
                 stringsAsFactors=FALSE)


# avalues <- c(0,5,10,15,20) # bvalues are the same as those are where elo-rating did to seem to do very well (see above plots)
# N.inds.values <- c(50)
# N.obs.values <- c(1,10,
#                   20,
#                   30,40,50)
# #to play around
avalues<-5
N.inds.values <- c(50)
N.obs.values <- c(1,5,9,
                  11,15,19,
                  21,25,29,
                  31,35,39,
                  41,45,49)
N.obs.values <- c(1,10,20,30,40,50)


for (j in 1:length(avalues)){
  
  for (p in 1:length(N.inds.values)){
    
    for (o in 1:length(N.obs.values)){
      
      for (sim in 1:10){
        
        output <- generate_interactions(N.inds.values[p],
                                        N.inds.values[p]*N.obs.values[o],
                                        a=avalues[j],
                                        b=avalues[j])
        
        winner <- output$interactions$Winner
        loser <- output$interactions$Loser
        date <- output$interactions$Date
        
        # generating elo.rating according to elo.seq from library(EloRating)
        x<-elo.seq(winner=as.factor(winner),
                   loser=as.factor(loser), 
                   Date=date,
                   k=200,
                   progressbar=FALSE)
        
        w<-extract.elo(x,standardize = FALSE)
        
        z <- data.frame(w,attributes(w),
                        row.names = as.character(seq(1,length(w),
                                                     1)))
        
        z$rank <- as.numeric(as.character(z$names))
        
        z$Elo.ranked <- rank(-z$w)
        
        spearman.original<-cor(z$rank,z$Elo.ranked,
                               use="complete.obs",method="spearman")
        
        # generating david's score
        dav<-DS(creatematrix(x, drawmethod="0.5"))
        
        dav$ID <- as.numeric(as.character(dav$ID))
        
        dav$normDSrank <- rank(-dav$normDS)
        
        Ndavid <- cor(dav$ID,dav$normDSrank,
                      use="complete.obs",method="spearman")
        
        # generating elo-rating according to elo.scores() from this script
        result.no.rand <- elo.scores(winner,loser,n.rands=1,
                                     init.score=1000,
                                     n.inds=N.inds.values[p])
        
        ranks.no.rand <- rank(-result.no.rand)
        
        spearman.cor.no.rand<-cor(output$hierarchy$Rank,
                                  ranks.no.rand,
                                  use="complete.obs",method="spearman")
        
        # generating elo-rating according to elo.scores() from this script and
        # randomizing the order of the interactions 1000 times
        result <- elo.scores(winner,loser,init.score=1000,
                             n.inds=N.inds.values[p])
        
        #mean.scores <- rowMeans(result)
        ranks <- apply(-result,2,function(x) rank(x, na.last="keep"))
        mean.ranks <- rowMeans(ranks)
        
        spearman.cor.rand<-cor(output$hierarchy$Rank,
                               mean.ranks,
                               use="complete.obs",method="spearman")
        
        #elochoice() no randomization
        
        eloc.1<-ratings(elochoice(winner,loser,
                                  kval=200,startvalue=1000,
                                  normprob=TRUE,runs=1),
                        drawplot=FALSE)
        
        z.eloc.1 <- data.frame(eloc.1,attributes(eloc.1),
                                       row.names = as.character(seq(1,length(eloc.1),
                                                             1)))
               
        z.eloc.1$rank <- as.numeric(as.character(z.eloc.1$names))
        
        z.eloc.1$Elo.ranked <- rank(-z.eloc.1$eloc.1)
        
        elochoice.no.rand<-cor(z.eloc.1$rank,z.eloc.1$Elo.ranked,
                               use="complete.obs",method="spearman")
        
        
        #elochoice() randomization
        
        eloc.2<-ratings(elochoice(winner,loser,
                                  kval=200,startvalue=1000,
                                  normprob=TRUE,runs=1000),
                        drawplot=FALSE)
        
        
        z.eloc.2 <- data.frame(eloc.2,attributes(eloc.2),
                                       row.names = as.character(seq(1,length(eloc.2),
                                                                    1)))
        
        z.eloc.2$rank <- as.numeric(as.character(z.eloc.2$names))
        
        z.eloc.2$Elo.ranked <- rank(-z.eloc.2$eloc.2)
        
        elochoice.rand<-cor(z.eloc.2$rank,z.eloc.2$Elo.ranked,
                            use="complete.obs",method="spearman")
        
      
        db<-rbind(db,c(N.inds.values[p],N.obs.values[o],
                       avalues[j],avalues[j],
                       Ndavid,spearman.original,
                       spearman.cor.no.rand,spearman.cor.rand,
                       elochoice.no.rand,elochoice.rand))
        
      }
    }
  }
}

names(db) <- c("Ninds","Nobs","alevel","blevel","Ndavid",
               "elo.original","elo.no.rand","elo.rand",
               "elochoice.no.rand","elochoice.rand")

proc.time() - ptm


# write.csv(db,
#           "db_1-10intandsoon.csv",row.names=FALSE)


for (p in 1:length(N.inds.values)){
  
  #par(mfrow=c(3,2))
  
  db.2 <- db[db$Ninds==N.inds.values[p],]
  
  for (i in 1:length(avalues)){
    
    db.3 <- db.2[db.2$alevel==avalues[i],]
    
    plot(db.3$elo.rand~db.3$Nobs,0.5,type="n",
         ylab="spearman correlation", 
         xlab="Number of interactions/individual",ylim=c(0,1), pch=19)
    
    points(db.3$Nobs,db.3$Ndavid,type="b",col="black")
    points(db.3$Nobs,db.3$elo.original,type="b",col="red",pch=19)
    points(db.3$Nobs,db.3$elo.no.rand,type="b",col="orange",pch=19)
    points(db.3$Nobs,db.3$elo.rand,type="b",col="blue")
    points(db.3$Nobs,db.3$elochoice.no.rand,type="b",col="pink")
    points(db.3$Nobs,db.3$elochoice.rand,type="b",col="green")
    
    Nindtext <- paste("N.ind = ",N.inds.values[p])
    atext <- paste("\na = ",avalues[i])
    btext <- paste("\nb = ",avalues[i])
    ttext <- paste0(Nindtext,atext,sep="\n")
    ttext2 <- paste0(ttext,btext,sep="\n")
    text(45,0.2,ttext2)
    
    par(xpd=TRUE)
    legend("bottom",c("Ndavid","Elo.original","Elo.no.rand",
                      "Elo.rand","elochoice.no.rand","elochoice.rand"),
           col=c("black","red","orange","blue","pink","green"),
             cex=1,bty='n',
             y.intersp=0.2,
             x.intersp=0.2,
             pch=rep(19,4))
#       legend("bottom",c("elo.seq()","our_elo_estimate"),
#            col=c("red","orange"),
#            cex=1,bty='n',
#            y.intersp=1,
#            x.intersp=0.2,
#            #horiz=TRUE,
#            pch=rep(19,2))
        
  }
  
  par(mfrow=c(1,1))
  
}


###############################################################################
# Plotting estimated rank ~ real rank: spearman correaltion for each method
# for the "tricky" scenarios. Adding 95% CI intervals 
###############################################################################
# db100sim <- read.table("db_100_simulations.csv",header=TRUE,sep=",")
# 
# avalues <- c(0,5,10,15,20) # bvalues are the same as those are where elo-rating did to seem to do very well (see above plots)
# N.inds.values <- c(50)
# N.obs.values <- c(1,10,
#                   20,
#                   30,40,50)

for (p in 1:length(N.inds.values)){
  
  #par(mfrow=c(3,2))
  
  db.2 <- db[db$Ninds==N.inds.values[p],]
  
  for (i in 1:length(avalues)){
    
    db.3 <- db.2[db.2$alevel==avalues[i],]
    
    db.4 <-summaryBy(Ndavid + elo.original + elo.no.rand +
                       elo.rand + elochoice.no.rand +
                       elochoice.rand ~ Nobs, 
                     data = db.3, 
                     FUN = function(x) { c(m = mean(x),
                                           q = quantile(x,probs=c(0.025,0.975))) })
    
    names(db.4) <- c("Nobs","Ndavid.m","Ndavid.lower","Ndavid.upper",
                     "elo.original.m","elo.original.lower","elo.original.upper",
                     "elo.no.rand.m","elo.no.rand.lower","elo.no.rand.upper",
                     "elo.rand.m","elo.rand.lower","elo.rand.upper",
                     "elochoice.no.rand.m","elochoice.no.rand.lower","elochoice.no.rand.upper",
                     "elochoice.rand.m","elochoice.rand.lower","elochoice.rand.upper")
    
    plot(db.4$elo.rand.m~db.4$Nobs,0.5,type="n",
         ylab="spearman correlation",
         xlab="number of interactions/individual",
         ylim=c(0,1))
        
    #adding points for the means and shadowed areas for the 95% CI
    points(db.4$Nobs,db.4$elo.original.m,type="b",col="red",pch=19)
    polygon(c(db.4$Nobs,rev(db.4$Nobs)),
            c(db.4$elo.original.lower,rev(db.4$elo.original.upper)),
            border=NA,col=rgb(1,0,0, 0.15))
    
    points(db.4$Nobs,db.4$elo.no.rand.m,type="b",col="orange",pch=19)
    polygon(c(db.4$Nobs,rev(db.4$Nobs)),
            c(db.4$elo.no.rand.lower,rev(db.4$elo.no.rand.upper)),
            border=NA,col=rgb(1,165/255,0,0.15))
   
    points(db.4$Nobs,db.4$Ndavid.m,type="b",col="black",pch=19)
    polygon(c(db.4$Nobs,rev(db.4$Nobs)),
           c(db.4$Ndavid.lower,rev(db.4$Ndavid.upper)),
           border=NA,col=rgb(120/255,120/255,120/255,0.15))
    
    points(db.4$Nobs,db.4$elo.rand.m,type="b",col="blue",pch=19)
    polygon(c(db.4$Nobs,rev(db.4$Nobs)),
            c(db.4$elo.rand.lower,rev(db.4$elo.rand.upper)),
            border=NA,col=rgb(0,0,1, 0.15))
    
    points(db.4$Nobs+0.05,db.4$elochoice.no.rand.m,type="b",col="pink",pch=19)
    polygon(c(db.4$Nobs+0.05,rev(db.4$Nobs+0.05)),
            c(db.4$elochoice.no.rand.lower,rev(db.4$elochoice.no.rand.upper)),
            border=NA,col=rgb(238/255,130/255,238/255, 0.15))
    
    points(db.4$Nobs,db.4$elochoice.rand.m,type="b",col="green",pch=19)
    polygon(c(db.4$Nobs,rev(db.4$Nobs)),
            c(db.4$elochoice.rand.lower,rev(db.4$elochoice.rand.upper)),
            border=NA,col=rgb(0,1,0, 0.15))
    
    Nindtext <- paste("N.ind = ",N.inds.values[p])
    atext <- paste("\na = ",avalues[i])
    btext <- paste("b = ",avalues[i])
    ttext <- paste0(Nindtext,atext,sep="\n")
    ttext2 <- paste0(ttext,btext,sep="\n")
    text(38,0.11,ttext2,adj = 0)
    
    par(xpd=TRUE)
    legend(#10,0.70,
           "bottom",
           c("David's score","Elo original","Elo no rand",
                     "Elo rand","elochoice.no.rand","elochoice.rand"),
           col=c("black","red","orange","blue","pink","green"),
           cex=1,bty='n',
           y.intersp=0.2,
           x.intersp=0.2,
           pch=rep(19,4),
           inset=c(0,0))
    
  }
  
  par(mfrow=c(1,1))
  
}



###############################################################################
# SECTION 3
###############################################################################
ptm <- proc.time()

db.split <- data.frame(Ninds=integer(),
                       Nobs=integer(),
                       alevel=integer(),
                       blevel=integer(),
                       #Ndavid.split=numeric(),
                       #elo.split=numeric(),
                       elo.rand.split=numeric(),
                       stringsAsFactors=FALSE)

# avalues <- seq(0,36,2)
# bvalues <- c(-5,0,5,10,15,20,25,30,35)
avalues <- c(0,5,10,15,20)
N.inds.values <- c(50)
N.obs.values <- c(1,2,3,4,5,6,7,8,9,10,20,30,40,50)


for (j in 1:length(avalues)){
  
  for (p in 1:length(N.inds.values)){
    
    for (o in 1:length(N.obs.values)){
      
      for (numsim in 1:100){
        
        
        output1 <- generate_interactions(N.inds.values[p],
                                         (N.inds.values[p])*(N.obs.values[o]/2),
                                         a=avalues[j],
                                         b=avalues[j])
        output2 <- generate_interactions(N.inds.values[p],
                                         (N.inds.values[p])*(N.obs.values[o]/2),
                                         a=avalues[j],
                                         b=avalues[j])
        
        winner1 <- output1$interactions$Winner
        loser1 <- output1$interactions$Loser
        winner2 <- output2$interactions$Winner
        loser2 <- output2$interactions$Loser
        date1 <- output1$interactions$Date
        date2 <- output2$interactions$Date
        
        #       result.no.rand1 <- elo.scores(winner1,loser1,n.rands=1)
        #       #result.no.rand1 <- as.data.table(cbind(result.no.rand1,c(1:nrow(result.no.rand1))))
        #       result.no.rand2 <- elo.scores(winner2,loser2,n.rands=1)
        #       #result.no.rand2 <- cbind(result.no.rand2,c(1:nrow(result.no.rand2)))
        #       
        #       #result.no.rand <- merge(result.no.rand1,result.no.rand2)
        #       
        #       elo.split<-cor(result.no.rand1,result.no.rand2,
        #                      use="complete.obs",method="spearman")
        #       
        
        result.rand1 <- elo.scores(winner1,loser1,init.score=1000,
                                   n.inds=N.inds.values[p])
        result.rand2 <- elo.scores(winner2,loser2,init.score=1000,
                                   n.inds=N.inds.values[p])
        
        ranks.rand1 <- apply(result.rand1,2,function(x) rank(x, na.last="keep"))
        mean.ranks.rand1 <- rowMeans(ranks.rand1)
        
        ranks.rand2 <- apply(result.rand2,2,function(x) rank(x, na.last="keep"))
        mean.ranks.rand2 <- rowMeans(ranks.rand2)
        
        elo.rand.split<-cor(mean.ranks.rand1,mean.ranks.rand2,
                            use="complete.obs",method="spearman")
        
        
        #       x.1<-elo.seq(winner=as.factor(winner1),
        #                    loser=as.factor(loser1), 
        #                    Date=date1,
        #                    k=200,
        #                    progressbar=FALSE)
        #       
        #       x.2<-elo.seq(winner=as.factor(winner2),
        #                    loser=as.factor(loser2), 
        #                    Date=date2,
        #                    k=200,
        #                    progressbar=FALSE)
        #                  
        #       
        #       dav.1<-DS(creatematrix(x.1, drawmethod="0.5"))
        #       dav.2<-DS(creatematrix(x.2, drawmethod="0.5"))
        #       
        #       dav.1$normDSrank.1 <- rank(dav.1$normDS)
        #       dav.2$normDSrank.2 <- rank(dav.2$normDS)
        #       
        #       dav <- merge(dav.1,dav.2,by="ID")
        #             
        #       dav$ID <- as.numeric(as.character(dav$ID))
        #       
        #           
        #       Ndavid <- cor(dav$normDSrank.1,dav$normDSrank.2,
        #                     use="complete.obs",method="spearman")
        
        db.split<-rbind(db.split,c(N.inds.values[p],N.obs.values[o],
                                   avalues[j],avalues[j],
                                   #Ndavid,elo.split,
                                   elo.rand.split))
        
        
      }
    }
  }
}

names(db.split) <- c("Ninds","Nobs","alevel","blevel",
                     #"Ndavid","elo",
                     "elo.rand")

proc.time() - ptm



for (p in 1:length(N.inds.values)){
  
  par(mfrow=c(3,2))
  
  db.2 <- db.split[db.split$Ninds==N.inds.values[p],]
  
  for (i in 1:length(avalues)){
    
    db.3 <- db.2[db.2$alevel==avalues[i],]
    
    plot(db.3$elo.rand~db.3$Nobs,0.5,type="n",
         ylab="spearman", xlab="Nint/individual",ylim=c(-1,1), pch=19)
    
    #points(db.3$Nobs,db.3$Ndavid,type="b",col="black")
    #points(db.3$Nobs,db.3$elo,type="b",col="orange")
    points(db.3$Nobs,db.3$elo.rand,type="b",col="blue")
    
    Nindtext <- paste("Nind = ",N.inds.values[p])
    atext <- paste("\na and b = ",avalues[i])
    ttext <- paste0(Nindtext,atext,sep="\n")
    text(45,-0.2,ttext)
    
#     par(xpd=TRUE)
#     legend("bottomright",c("Ndavid","Elo","Elo.rand"),
#            col=c("black","orange","blue"),
#            cex=1,bty='n',
#            y.intersp=0.2,
#            x.intersp=0.2,
#            #horiz=TRUE,
#            pch=rep(19,3))
    
  }
  
  par(mfrow=c(1,1))
  
}