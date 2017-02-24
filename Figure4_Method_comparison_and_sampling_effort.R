
# Authors: 
#         Alfredo Sanchez-Tojar, MPIO (Seewiesen) and ICL (Silwood Park), alfredo.tojar@gmail.com
#         Damien Farine, MPIO (Radolfzell) and EGI (Oxford), dfarine@orn.mpg.de

# Script first created on the 27th of October, 2016

###############################################################################
# Description of script and Instructions
###############################################################################

# This script aims to generate a figure that shows how different dominance
# methods perform inferring the latent hierachy, Figure 3 for:

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
library(EloRating)


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
# Simulating: COMPARING all METHODS
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
        
        for (sim in 1:100){
          
          output <- generate_interactions(N.inds.values[p],
                                          N.inds.values[p]*N.obs.values[o],
                                          a=avalues[j],
                                          b=bvalues[j],
                                          id.biased=poiss[typ],
                                          rank.biased=dombias[typ])
          
          
          filename <- paste(paste(paste(ifelse(poiss[typ]==TRUE,1,0),
                                        ifelse(dombias[typ]==TRUE,1,0),
                                        sep=""),
                                  paste0(bvalues[j],a=avalues[j]),
                                  sep="_"),
                            N.obs.values[o],
                            sep="_")
          
          #saving hierarchies
          # write.csv(output$hierarchy,
          #           paste0("databases_package/matrices_for_ADAGIO",
          #                  typ,"/hierarchies/hierarchy_",filename,"_sim",sim,".csv"),
          #           row.names=FALSE)
          
          winner <- output$interactions$Winner
          loser <- output$interactions$Loser
          date <- c("2010-01-25",rep("2010-01-26",length(loser)-1)) #fake date needed for elo.seq()
          
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
          
          z$Elo.ranked <- rank(-z$w,na.last="keep")
          
          spearman.original<-cor(z$rank,z$Elo.ranked,
                                 use="complete.obs",method="spearman")
          
          # generating David's score
          domatrix<-creatematrix(x, drawmethod="0.5")
          
          # saving matrices for running ADAGIO
          write.csv(as.data.frame(domatrix),
                    paste0("databases_package/matrices_10ind",
                           typ,"/matrices/matrix_",filename,"_sim",sim,".csv"),
                    row.names=TRUE,
                    quote = FALSE)
          
          dav<-DS(domatrix)
          
          dav$ID <- as.numeric(as.character(dav$ID))
          
          dav$normDSrank <- rank(-dav$normDS,na.last="keep")
          
          Ndavid <- cor(dav$ID,dav$normDSrank,
                        use="complete.obs",method="spearman")
          
          # generating elo-rating according to elo_scores() from aniDom
          result.no.rand <- elo_scores(winner,
                                       loser,
                                       identities=c(1:N.inds.values[p]),
                                       init.score=1000,
                                       randomise=FALSE,
                                       return.as.ranks=TRUE)
          
          #ranks.no.rand <- rank(-result.no.rand,na.last="keep")
          
          spearman.cor.no.rand<-cor(output$hierarchy$Rank,
                                    #ranks.no.rand,
                                    result.no.rand,
                                    use="complete.obs",method="spearman")
          
          # generating elo-rating according to elo_scores() from aniDom and
          # randomizing the order of the interactions 1000 times
          result <- elo_scores(winner,
                               loser,
                               identities=c(1:N.inds.values[p]),
                               init.score=1000,
                               randomise=TRUE,
                               return.as.ranks=TRUE)
          
          mean.scores <- rowMeans(result)
          #ranks <- apply(-result,2,function(x) rank(x, na.last="keep"))
          #mean.ranks <- rowMeans(ranks)
          
          spearman.cor.rand<-cor(output$hierarchy$Rank,
                                 mean.scores,
                                 use="complete.obs",method="spearman")
          
          db<-rbind(db,c(N.inds.values[p],
                         N.obs.values[o],
                         poiss[typ],
                         dombias[typ],
                         avalues[j],
                         bvalues[j],
                         Ndavid,
                         spearman.original,
                         spearman.cor.no.rand,
                         spearman.cor.rand))
          
        }
      }
    }
  }
  
}

names(db) <- c("Ninds","Nobs",
               "poiss","dombias",
               "alevel","blevel",
               "Ndavid",
               "elo.original",
               "elo.no.rand",
               "elo.rand")

proc.time() - ptm


# write.csv(db,
#           "databases_package/Fig3_db_methods_100sim_fixed_biases.csv",row.names=FALSE)

write.csv(db,
          "databases_package/Fig3_db_methods_100sim_fixed_biases_10ind.csv",row.names=FALSE)


###############################################################################
# ADAGIO: Adding the results from ADAGIO, which needed to be run at the terminal
###############################################################################

# importing all rank files (output from ADAGIO)

temp = list.files(pattern="*.csv.adagio.ranks") #check you are in the right folder - getwd() and setwd()

setwd("C:/Users/nb5093/Documents/allmatrices")

db <- data.frame(Ninds=integer(),
                 Nobs=integer(),
                 poiss=integer(),
                 dombias=integer(),
                 alevel=integer(),
                 blevel=integer(),
                 spearman=numeric(),
                 stringsAsFactors=FALSE)


#N.inds.values <- c(50)
N.inds.values <- c(10)


for (filename in 1:length(temp)){
  
  poiss <- substr(temp[filename], 8, 8)
  dombias <- substr(temp[filename], 9, 9)
  
  if(substr(temp[filename], 11, 12)=="-5"){
    
    blevel <- substr(temp[filename], 11, 12)
    
    if(substr(temp[filename], 14, 14)=="_"){
      
      alevel <- substr(temp[filename], 13, 13)
      
      if(substr(temp[filename], 16, 16)=="_"){
        
        Nobs <- substr(temp[filename], 15, 15)
        
      } else {
        
        Nobs <- substr(temp[filename], 15, 16)
        
      }
      
    } else {
      
      alevel <- substr(temp[filename], 13, 14)
      
      if(substr(temp[filename], 17, 17)=="_"){
        
        Nobs <- substr(temp[filename], 16, 16)
        
      } else {
        
        Nobs <- substr(temp[filename], 16, 17)
        
      }
      
    }
    
  } else {
    
    blevel <- substr(temp[filename], 11, 11)
    
    if(substr(temp[filename], 13, 13)=="_"){
      
      alevel <- substr(temp[filename], 12, 12)
      
      if(substr(temp[filename], 15, 15)=="_"){
        
        Nobs <- substr(temp[filename], 14, 14)
        
      } else {
        
        Nobs <- substr(temp[filename], 14, 15)
        
      }
      
    } else {
      
      alevel <- substr(temp[filename], 12, 13)
      
      if(substr(temp[filename], 16, 16)=="_"){
        
        Nobs <- substr(temp[filename], 15, 15)
        
      } else {
        
        Nobs <- substr(temp[filename], 15, 16)
        
      }
      
    }
    
  } 
  
  db_temp <- read.table(temp[filename],header=FALSE,sep="\t")
  
  spearman.cor<-cor(db_temp$V1,
                    db_temp$V2,
                    use="complete.obs",method="spearman")
  
  
  db <- rbind(db,as.numeric(c(N.inds.values,
                              Nobs,
                              poiss,
                              dombias,
                              alevel,
                              blevel,
                              spearman.cor)))
  
}

names(db) <- c("Ninds","Nobs",
               "poiss","dombias",
               "alevel","blevel","spearman")


#write.csv(db,
#          "databases_package/db_ADAGIO_100simulations_fixed_biases.csv",row.names=FALSE)

write.csv(db,
          "databases_package/final_data_for_Figures_backup/db_ADAGIO_100simulations_fixed_biases_10ind.csv",row.names=FALSE)


###############################################################################
# Plotting: MAIN TEXT: estimated rank ~ real rank: spearman correaltion for each 
# method. Adding 95% CI intervals 
###############################################################################

db5methods <- read.table("databases_package/final_data_for_Figures_backup/Fig4a_db_methods_100sim_fixed_biases.csv",
                                  header=TRUE,sep=",")


db5methods_sorted <- db5methods[order(db5methods$Ninds,
                                      db5methods$Nobs,
                                      db5methods$poiss,
                                      db5methods$dombias,
                                      db5methods$blevel,
                                      db5methods$alevel),]


dbADAGIO <- read.table("databases_package/final_data_for_Figures_backup/Fig4b_db_ADAGIO_100simulations_fixed_biases.csv",
                       header=TRUE,sep=",")


dbADAGIO_sorted <- dbADAGIO[order(dbADAGIO$Ninds,
                                  dbADAGIO$Nobs,
                                  dbADAGIO$poiss,
                                  dbADAGIO$dombias,
                                  dbADAGIO$blevel,
                                  dbADAGIO$alevel),]

# making sure everything is identical

iden <- as.factor(c("Nobs","poiss","dombias","alevel","blevel"))


for (column in levels(iden)){
  
  print(identical(db5methods_sorted[,c(column)],
            dbADAGIO_sorted[,c(column)],
            attrib.as.set=FALSE))
  
}


# cbinding both databases to add the ADAGIO measurements

dbADAGIO_cor <- dbADAGIO_sorted[,c("blevel", "spearman")]
names(dbADAGIO_cor) <- c("blevel","ADAGIO")


db.provisional <- cbind(db5methods_sorted,dbADAGIO_cor)


db.provisional.2 <- db.provisional[,c(1:11,13)]


# Database is ready, plotting starts!

avalues <- c(15,15,15,10,10,10,5,5,5)
bvalues <- c(5,5,5,5,5,5,5,5,5)
N.inds.values <- c(50)
N.obs.values <- c(1,4,7,10,15,20,30,40,50)


db<-db.provisional.2[db.provisional.2$poiss==1 & db.provisional.2$dombias==0,]


a <- c("(a)","x","x","(b)","x","x","(c)","x","x")


tiff(#"plots/Figure4_Method_comparison_and_sampling_effort_poisson.tiff",
     "plots/supplements/FigureS1_Comparing_original_Elo-rating_packages_poisson.tiff",
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

  
  db.2 <- db[db$Ninds==N.inds.values[p] & (db$Nobs %in% N.obs.values),]
  
  for (i in 1:length(avalues)){
    
    if(i %in% c(1,4,7)){
      
      db.3 <- db.2[db.2$alevel==avalues[i] & db.2$blevel==bvalues[i],]
      
      db.4 <-summaryBy(Ndavid + 
                         elo.original + 
                         elo.no.rand +
                         elo.rand + 
                         ADAGIO
                       ~ Nobs, 
                       data = db.3, 
                       FUN = function(x) { c(m = mean(x),
                                             q = quantile(x,probs=c(0.025,0.975))) })
      
      names(db.4) <- c("Nobs",
                       "Ndavid.m","Ndavid.lower","Ndavid.upper",
                       "elo.original.m","elo.original.lower","elo.original.upper",
                       "elo.no.rand.m","elo.no.rand.lower","elo.no.rand.upper",
                       "elo.rand.m","elo.rand.lower","elo.rand.upper",
                       "ADAGIO.m","ADAGIO.lower","ADAGIO.upper")
      
      plot(db.4$elo.rand.m~db.4$Nobs,0.5,type="n",
           ylab="",
           xlab="",
           xaxt="n",
           yaxt="n",
           ylim=c(-0.3,1))
      
    
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
    
    axis(2,at=round(seq(-0.3,1,0.1),1),cex.axis=1.2,las=2,tck=0.015)

    
    #adding points for the means and shadowed areas for the 95% CI
    points(db.4$Nobs,db.4$elo.original.m,type="b",col="green",pch=19)
    polygon(c(db.4$Nobs,rev(db.4$Nobs)),
            c(db.4$elo.original.lower,rev(db.4$elo.original.upper)),
            border=NA,col=rgb(0,1,0, 0.15))
    
    points(db.4$Nobs,db.4$elo.no.rand.m,type="b",col="red",pch=19)
    polygon(c(db.4$Nobs,rev(db.4$Nobs)),
            c(db.4$elo.no.rand.lower,rev(db.4$elo.no.rand.upper)),
            border=NA,col=rgb(1,0,0,0.15))
    
    # points(db.4$Nobs,db.4$Ndavid.m,type="b",col="black",pch=19)
    # polygon(c(db.4$Nobs,rev(db.4$Nobs)),
    #         c(db.4$Ndavid.lower,rev(db.4$Ndavid.upper)),
    #         border=NA,col=rgb(120/255,120/255,120/255,0.15))
    # 
    # points(db.4$Nobs,db.4$elo.rand.m,type="b",col="blue",pch=19)
    # polygon(c(db.4$Nobs,rev(db.4$Nobs)),
    #         c(db.4$elo.rand.lower,rev(db.4$elo.rand.upper)),
    #         border=NA,col=rgb(0,0,1, 0.15))
    # 
    # points(db.4$Nobs,db.4$ADAGIO.m,type="b",col="orange",pch=19)
    # polygon(c(db.4$Nobs,rev(db.4$Nobs)),
    #         c(db.4$ADAGIO.lower,rev(db.4$ADAGIO.upper)),
    #         border=NA,col=rgb(1,165/255,0,0.15))
    
    lines(c(0,51),c(0.7,0.7),col="red",lty=3,lwd=1.5)
    
    atext <- paste("\na = ",avalues[i])
    btext <- paste("\nb = ",bvalues[i])
    ttext <- paste0(atext,btext,sep="\n")
    text(48,-0.25,a[i],adj = 0 ,cex=1.5)
    text(31,-0.2,ttext,adj = 0,cex=2)

    
    } else {
      
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
        
        par(xpd=TRUE)
        # legend(0,0.8,
        #        c("David's score",
        #          "original Elo-rating",
        #          "randomized Elo-rating",
        #          "ADAGIO"),
        #        col=c("black",
        #              "red",
        #              "blue",
        #              "orange"),
        #        cex=1.35,bty='n',
        #        pch=rep(19,4),
        #        inset=c(0,0))
        
        legend(0,0.8,
               c("package:aniDom",
                 "package:EloRating"),
               col=c("red",
                     "green"),
               cex=1.35,bty='n',
               pch=rep(19,4),
               inset=c(0,0))
        
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
# Plotting: SUPPLEMENTARY MATERIAL, part 1: 
# intermediate scenarios: uniform, dominant bias, poisson*dominant bias
###############################################################################

avalues <- c(15,15,15,10,10,10,5,5,5)
bvalues <- c(5,5,5,5,5,5,5,5,5)
N.inds.values <- c(50)
N.obs.values <- c(1,4,7,10,15,20,30,40,50)


#db<-db.provisional.2[db.provisional.2$poiss==0 & db.provisional.2$dombias==0,]
#db<-db.provisional.2[db.provisional.2$poiss==0 & db.provisional.2$dombias==1,]
db<-db.provisional.2[db.provisional.2$poiss==1 & db.provisional.2$dombias==1,]


a <- c("(a)","x","x","(b)","x","x","(c)","x","x")


tiff(#"plots/supplements/FigureS2_Method_comparison_and_sampling_effort_uniform.tiff",
     #"plots/supplements/FigureS3_Comparing_original_Elo-rating_packages_uniform.tiff",
     #"plots/supplements/FigureS9_Method_comparison_and_sampling_effort_dombias.tiff",
     #"plots/supplements/FigureSX_Comparing_original_Elo-rating_packages_dombias.tiff",
     #"plots/supplements/FigureS16_Method_comparison_and_sampling_effort_poiss+dombias.tiff",
     "plots/supplements/FigureSX_Comparing_original_Elo-rating_packages_poiss+dombias.tiff",
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
  
  
  db.2 <- db[db$Ninds==N.inds.values[p] & (db$Nobs %in% N.obs.values),]
  
  for (i in 1:length(avalues)){
    
    if(i %in% c(1,4,7)){
      
      db.3 <- db.2[db.2$alevel==avalues[i] & db.2$blevel==bvalues[i],]
      
      db.4 <-summaryBy(Ndavid + 
                         elo.original + 
                         elo.no.rand +
                         elo.rand + 
                         ADAGIO
                       ~ Nobs, 
                       data = db.3, 
                       FUN = function(x) { c(m = mean(x),
                                             q = quantile(x,probs=c(0.025,0.975))) })
      
      names(db.4) <- c("Nobs",
                       "Ndavid.m","Ndavid.lower","Ndavid.upper",
                       "elo.original.m","elo.original.lower","elo.original.upper",
                       "elo.no.rand.m","elo.no.rand.lower","elo.no.rand.upper",
                       "elo.rand.m","elo.rand.lower","elo.rand.upper",
                       "ADAGIO.m","ADAGIO.lower","ADAGIO.upper")
      
      plot(db.4$elo.rand.m~db.4$Nobs,0.5,type="n",
           ylab="",
           xlab="",
           xaxt="n",
           yaxt="n",
           ylim=c(-0.3,1))
      
      
      if(i<7){
        
        axis(1,at=N.obs.values,
             cex.axis=1,tck=0.015,
             labels=FALSE)
        
        
      } else {
        
        axis(1,at=N.obs.values,
             labels=as.character(N.obs.values),
             cex.axis=1,tck=0.015)
        
        mtext("ratio of interactions to individuals",
              side=1, adj=0, line=4, cex=1.8); 
        
      }
      
      axis(2,at=round(seq(-0.3,1,0.1),1),cex.axis=1.2,las=2,tck=0.015)
      
      
      #adding points for the means and shadowed areas for the 95% CI
      points(db.4$Nobs,db.4$elo.original.m,type="b",col="green",pch=19)
      polygon(c(db.4$Nobs,rev(db.4$Nobs)),
              c(db.4$elo.original.lower,rev(db.4$elo.original.upper)),
              border=NA,col=rgb(0,1,0, 0.15))
      
      points(db.4$Nobs,db.4$elo.no.rand.m,type="b",col="red",pch=19)
      polygon(c(db.4$Nobs,rev(db.4$Nobs)),
              c(db.4$elo.no.rand.lower,rev(db.4$elo.no.rand.upper)),
              border=NA,col=rgb(1,0,0,0.15))
      
      # points(db.4$Nobs,db.4$Ndavid.m,type="b",col="black",pch=19)
      # polygon(c(db.4$Nobs,rev(db.4$Nobs)),
      #         c(db.4$Ndavid.lower,rev(db.4$Ndavid.upper)),
      #         border=NA,col=rgb(120/255,120/255,120/255,0.15))
      # 
      # points(db.4$Nobs,db.4$elo.rand.m,type="b",col="blue",pch=19)
      # polygon(c(db.4$Nobs,rev(db.4$Nobs)),
      #         c(db.4$elo.rand.lower,rev(db.4$elo.rand.upper)),
      #         border=NA,col=rgb(0,0,1, 0.15))
      # 
      # points(db.4$Nobs,db.4$ADAGIO.m,type="b",col="orange",pch=19)
      # polygon(c(db.4$Nobs,rev(db.4$Nobs)),
      #         c(db.4$ADAGIO.lower,rev(db.4$ADAGIO.upper)),
      #         border=NA,col=rgb(1,165/255,0,0.15))
      
      lines(c(0,51),c(0.7,0.7),col="red",lty=3,lwd=1.5)
      
      atext <- paste("\na = ",avalues[i])
      btext <- paste("\nb = ",bvalues[i])
      ttext <- paste0(atext,btext,sep="\n")
      text(48,-0.25,a[i],adj = 0 ,cex=1.5)
      text(31,-0.2,ttext,adj = 0,cex=2)
      
      
    } else {
      
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
        
        par(xpd=TRUE)
        # legend(0,0.8,
        #        c("David's score",
        #          "original Elo-rating",
        #          "randomized Elo-rating",
        #          "ADAGIO"),
        #        col=c("black",
        #              "red",
        #              "blue",
        #              "orange"),
        #        cex=1.35,bty='n',
        #        pch=rep(19,4),
        #        inset=c(0,0))
        
        legend(0,0.8,
               c("package:aniDom",
                 "package:EloRating"),
               col=c("red",
                     "green"),
               cex=1.35,bty='n',
               pch=rep(19,4),
               inset=c(0,0))
        
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
# Plotting: SUPPLEMENTARY MATERIAL, part 2: 
# steep and flat scenarios
###############################################################################

avalues <- c(10,10,10,15,15,15,30,30,30,0,0,0)
bvalues <- c(-5,-5,-5,0,0,0,5,5,5,5,5,5)
N.inds.values <- c(50)
N.obs.values <- c(1,4,7,10,15,20,30,40,50)


#db<-db.provisional.2[db.provisional.2$poiss==1 & db.provisional.2$dombias==0,]
db<-db.provisional.2[db.provisional.2$poiss==0 & db.provisional.2$dombias==0,]
#db<-db.provisional.2[db.provisional.2$poiss==0 & db.provisional.2$dombias==1,]
#db<-db.provisional.2[db.provisional.2$poiss==1 & db.provisional.2$dombias==1,]


a <- c("(a)","x","x","(b)","x","x","(c)","x","x","(d)","x","x")


tiff(#"plots/supplements/FigureS1_Method_comparison_and_sampling_effort_steep_and_flat_poisson.tiff",
     #"plots/supplements/FigureS2_Comparing_original_Elo-rating_packages_steep_and_flat_poisson.tiff",
     #"plots/supplements/FigureS3_Method_comparison_and_sampling_effort_steep_and_flat_uniform.tiff",
     "plots/supplements/FigureS4_Comparing_original_Elo-rating_packages_steep_and_flat_uniform.tiff",
     #"plots/supplements/FigureS10_Method_comparison_and_sampling_effort_steep_and_flat_dombias.tiff",
     #"plots/supplements/FigureSX_Comparing_original_Elo-rating_packages_steep_and_flat_dombias.tiff",
     #"plots/supplements/FigureS17_Method_comparison_and_sampling_effort_steep_and_flat_poiss+dombias.tiff",
     #"plots/supplements/FigureSX_Comparing_original_Elo-rating_packages_steep_and_flat_poiss+dombias.tiff",
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
  
  
  db.2 <- db[db$Ninds==N.inds.values[p] & (db$Nobs %in% N.obs.values),]
  
  for (i in 1:length(avalues)){
    
    if(i %in% c(1,4,7,10)){
      
      db.3 <- db.2[db.2$alevel==avalues[i] & db.2$blevel==bvalues[i],]
      
      db.4 <-summaryBy(Ndavid + 
                         elo.original + 
                         elo.no.rand +
                         elo.rand + 
                         ADAGIO
                       ~ Nobs, 
                       data = db.3, 
                       FUN = function(x) { c(m = mean(x),
                                             q = quantile(x,probs=c(0.025,0.975))) })
      
      names(db.4) <- c("Nobs",
                       "Ndavid.m","Ndavid.lower","Ndavid.upper",
                       "elo.original.m","elo.original.lower","elo.original.upper",
                       "elo.no.rand.m","elo.no.rand.lower","elo.no.rand.upper",
                       "elo.rand.m","elo.rand.lower","elo.rand.upper",
                       "ADAGIO.m","ADAGIO.lower","ADAGIO.upper")
      
      plot(db.4$elo.rand.m~db.4$Nobs,0.5,type="n",
           ylab="",
           xlab="",
           xaxt="n",
           yaxt="n",
           ylim=c(-0.6,1))
      
      
      if(i<10){
        
        axis(1,at=N.obs.values,
             cex.axis=1,tck=0.015,
             labels=FALSE)
        
        
      } else {
        
        axis(1,at=N.obs.values,
             labels=as.character(N.obs.values),
             cex.axis=1,tck=0.015)
        
        mtext("ratio of interactions to individuals",
              side=1, adj=0, line=4, cex=1.8); 
        
      }
      
      axis(2,at=round(seq(-0.6,1,0.1),1),cex.axis=1.2,las=2,tck=0.015)
      
      
      #adding points for the means and shadowed areas for the 95% CI
      points(db.4$Nobs,db.4$elo.original.m,type="b",col="green",pch=19)
      polygon(c(db.4$Nobs,rev(db.4$Nobs)),
              c(db.4$elo.original.lower,rev(db.4$elo.original.upper)),
              border=NA,col=rgb(0,1,0, 0.15))
      
      points(db.4$Nobs,db.4$elo.no.rand.m,type="b",col="red",pch=19)
      polygon(c(db.4$Nobs,rev(db.4$Nobs)),
              c(db.4$elo.no.rand.lower,rev(db.4$elo.no.rand.upper)),
              border=NA,col=rgb(1,0,0,0.15))
      
      # points(db.4$Nobs,db.4$Ndavid.m,type="b",col="black",pch=19)
      # polygon(c(db.4$Nobs,rev(db.4$Nobs)),
      #         c(db.4$Ndavid.lower,rev(db.4$Ndavid.upper)),
      #         border=NA,col=rgb(120/255,120/255,120/255,0.15))
      # 
      # points(db.4$Nobs,db.4$elo.rand.m,type="b",col="blue",pch=19)
      # polygon(c(db.4$Nobs,rev(db.4$Nobs)),
      #         c(db.4$elo.rand.lower,rev(db.4$elo.rand.upper)),
      #         border=NA,col=rgb(0,0,1, 0.15))
      # 
      # points(db.4$Nobs,db.4$ADAGIO.m,type="b",col="orange",pch=19)
      # polygon(c(db.4$Nobs,rev(db.4$Nobs)),
      #         c(db.4$ADAGIO.lower,rev(db.4$ADAGIO.upper)),
      #         border=NA,col=rgb(1,165/255,0,0.15))
      
      lines(c(0,51),c(0.7,0.7),col="red",lty=3,lwd=1.5)
      
      atext <- paste("\na = ",avalues[i])
      btext <- paste("\nb = ",bvalues[i])
      ttext <- paste0(atext,btext,sep="\n")
      text(48,-0.55,a[i],adj = 0 ,cex=1.5)
      text(31,-0.45,ttext,adj = 0,cex=2)
      
      
    } else {
      
      if(i %in% c(2,5,8,11)){
        
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
              side=2, adj=0, line=3, cex=0.85); 
        
      } else {
        
        plot(c(1,50),c(0,1),type="n",
             ylab="", 
             xlab="",
             xaxt="n",
             yaxt="n",
             frame.plot=FALSE)    
        
        par(xpd=TRUE)
        # legend(0,0.7,
        #        c("David's score",
        #          "original Elo-rating",
        #          "randomized Elo-rating",
        #          "ADAGIO"),
        #        col=c("black",
        #              "red",
        #              "blue",
        #              "orange"),
        #        cex=1.35,bty='n',
        #        pch=rep(19,4),
        #        inset=c(0,0))
        
        legend(0,0.8,
               c("package:aniDom",
                 "package:EloRating"),
               col=c("red",
                     "green"),
               cex=1.35,bty='n',
               pch=rep(19,4),
               inset=c(0,0))
        
        mtext("Difference in rank      ",
              side=3, adj=1, line=-2, cex=1); 
        
      }
      
    }
    
  }
  
  
  title(ylab = "Spearman rank correlation coefficient",
        outer = TRUE, line = 0)
  
  par(mfrow=c(1,1))
  
}

dev.off()