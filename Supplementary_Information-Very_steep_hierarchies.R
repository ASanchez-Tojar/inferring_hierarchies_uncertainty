
# Authors: 
#         Alfredo Sanchez-Tojar, MPIO (Seewiesen) and ICL (Silwood Park), alfredo.tojar@gmail.com
#         Damien Farine, MPIO (Radolfzell) and EGI (Oxford), dfarine@orn.mpg.de

# Script first created on the 8th of June, 2017

###############################################################################
# Description of script and Instructions
###############################################################################

# This script aims to generate a figure that shows how different dominance
# methods perform inferring the latent hierachy, Figure X for:

# Sanchez-Tojar, A., Schroeder, J., Farine, D.R. (In preparation) A practical 
# guide for inferring reliable dominance hierarchies and estimating their 
# uncertainty

# more info at the Open Science Framework: http://doi.org/10.17605/OSF.IO/9GYEK


###############################################################################
# Packages needed
###############################################################################

# packages needed for this script

library(aniDom)
library(doBy)
library(RColorBrewer)
library(EloRating)


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
# Plotting: Supplementary Information 3
###############################################################################

# db5methods <- read.table("databases_package/final_data_for_Figures_backup/Fig5a_db_methods_100sim_fixed_biases_ISIincluded_10ind.csv",
#                          header=TRUE,sep=",")

db5methods <- read.table("databases_package/final_data_for_Figures_backup/Fig5a_db_methods_100sim_fixed_biases_ISIincluded_25ind.csv",
                         header=TRUE,sep=",")

# db5methods <- read.table("databases_package/final_data_for_Figures_backup/Fig5a_db_methods_100sim_fixed_biases_ISIincluded_50ind.csv",
#                          header=TRUE,sep=",")


db5methods_sorted <- db5methods[order(db5methods$Ninds,
                                      db5methods$Nobs,
                                      db5methods$poiss,
                                      db5methods$dombias,
                                      db5methods$blevel,
                                      db5methods$alevel),]


# dbADAGIO <- read.table("databases_package/final_data_for_Figures_backup/Fig5b_db_ADAGIO_100simulations_fixed_biases_10ind_100int.csv",
#                        header=TRUE,sep=",")

dbADAGIO <- read.table("databases_package/final_data_for_Figures_backup/Fig5b_db_ADAGIO_100simulations_fixed_biases_25ind_100int.csv",
                       header=TRUE,sep=",")

# dbADAGIO <- read.table("databases_package/final_data_for_Figures_backup/Fig5b_db_ADAGIO_100simulations_fixed_biases_50ind_100int.csv",
#                        header=TRUE,sep=",")


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


db.provisional.2 <- db.provisional[db.provisional$poiss==1 & 
                                     #db.provisional$poiss==0 & 
                                     db.provisional$dombias==0 &
                                     db.provisional$blevel==-5,
                                   c(1:15,17)]



# # Database is ready, plotting starts!

N.obs.values <- c(1,4,7,10,15,20,30,40,50,100)


tiff("plots/supplements/after_revision/FigureS1_Very_steep_hierarchies_25ind_Poisson.tiff",
     #"plots/supplements/after_revision/FigureS5_Very_steep_hierarchies_25ind_uniform.tiff",
     #"plots/supplements/after_revision/FigureS07_Very_steep_hierarchies_10ind_Poisson.tiff",
     #"plots/supplements/after_revision/FigureS14_Very_steep_hierarchies_50ind_Poisson.tiff",
     height=29.7, width=21,
     units='cm', compression="lzw", res=600)

m <- rbind(c(1,1,1,2),c(1,1,1,3),
           c(4,4,4,5),c(4,4,4,6),
           c(7,7,7,8),c(7,7,7,9))

layout(m)

op <- par(oma = c(6,3,1,1) + 0.1,
          mar = c(0.5,5,1,0) + 0.1,
          cex.lab=2.5)


db.sum <-summaryBy(Ndavid + 
                   elo.original + 
                   elo.no.rand +
                   elo.rand + 
                   ADAGIO +
                   ISI98
                 ~ Nobs, 
                 data = db.provisional.2, 
                 FUN = function(x) { c(m = mean(x),
                                       q = quantile(x,probs=c(0.025,0.975))) })

names(db.sum) <- c("Nobs",
                 "Ndavid.m","Ndavid.lower","Ndavid.upper",
                 "elo.original.m","elo.original.lower","elo.original.upper",
                 "elo.no.rand.m","elo.no.rand.lower","elo.no.rand.upper",
                 "elo.rand.m","elo.rand.lower","elo.rand.upper",
                 "ADAGIO.m","ADAGIO.lower","ADAGIO.upper",
                 "ISI98.m","ISI98.lower","ISI98.upper")

plot(db.sum$elo.rand.m~db.sum$Nobs,0.5,type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(0,1))

axis(1,at=N.obs.values,
     cex.axis=1,tck=0.015,
     labels=FALSE)

axis(2,
     at=round(seq(0,1,0.1),1),
     cex.axis=1,las=2,tck=0.015)


points(db.sum$Nobs,db.sum$ISI98.m,type="b",col="green",pch=19)
polygon(c(db.sum$Nobs,rev(db.sum$Nobs)),
        c(db.sum$ISI98.lower,rev(db.sum$ISI98.upper)),
        border=NA,col=rgb(0,1,0, 0.15))

points(db.sum$Nobs,db.sum$elo.no.rand.m,type="b",col="red",pch=19)
polygon(c(db.sum$Nobs,rev(db.sum$Nobs)),
        c(db.sum$elo.no.rand.lower,rev(db.sum$elo.no.rand.upper)),
        border=NA,col=rgb(1,0,0,0.15))

points(db.sum$Nobs,db.sum$Ndavid.m,type="b",col="black",pch=19)
polygon(c(db.sum$Nobs,rev(db.sum$Nobs)),
        c(db.sum$Ndavid.lower,rev(db.sum$Ndavid.upper)),
        border=NA,col=rgb(120/255,120/255,120/255,0.15))

points(db.sum$Nobs,db.sum$elo.rand.m,type="b",col="blue",pch=19)
polygon(c(db.sum$Nobs,rev(db.sum$Nobs)),
        c(db.sum$elo.rand.lower,rev(db.sum$elo.rand.upper)),
        border=NA,col=rgb(0,0,1, 0.15))

points(db.sum$Nobs,db.sum$ADAGIO.m,type="b",col="orange",pch=19)
polygon(c(db.sum$Nobs,rev(db.sum$Nobs)),
        c(db.sum$ADAGIO.lower,rev(db.sum$ADAGIO.upper)),
        border=NA,col=rgb(1,165/255,0,0.15))

lines(c(0,101),c(0.7,0.7),col="red",lty=3,lwd=1.5)

text(98,0.05,"(a)",adj = 0 ,cex=1.5)

mtext("  correlation coefficient (rs)",
      side=2, adj=0, line=4, cex=1.6)

text(81,0.1,"a = 10\nb = -5",adj = 0,cex=2)


# plotting the scenario

plot(c(1,25),
     #c(1,10),
     #c(1,50),
     c(0.5,1),type="n",
     ylab="", 
     xlab="",
     xaxt="n",
     yaxt="n",
     cex=1)

axis(1,
     at=seq(1,25,2),
     #at=seq(1,10,1),
     #at=seq(1,50,6),
     cex.axis=0.75,tck=0.015)

axis(2,at=seq(0.5,1,0.1),cex.axis=0.75,las=2,tck=0.015) 

plot_winner_prob(1:25,a=10,b=-5,"black")
#plot_winner_prob(1:10,a=10,b=-5,"black")
#plot_winner_prob(1:50,a=10,b=-5,"black")


mtext("P (dominant wins)",
      side=2, adj=0, line=3, cex=1.10)


# ADDING LEGEND

plot(c(1,50),c(0,1),type="n",
     ylab="", 
     xlab="",
     xaxt="n",
     yaxt="n",
     frame.plot=FALSE)  

mtext("Difference in rank",
      side=3, adj=1, line=-2, cex=0.95)

par(xpd=TRUE)
legend(-20,0.8,
       c("David's score",
         "original Elo-rating",
         "randomized Elo-rating",
         "ADAGIO",
         "I&SI"),
       col=c("black",
             "red",
             "blue",
             "orange",
             "green"),
       cex=1.25,bty='n',
       pch=rep(19,4),
       inset=c(0,0),
       title.adj=0)




# REPEATABILITY PLOT

#db_rep <- read.table("databases_package/final_data_for_Figures_backup/Fig6_db_repeatability_100_simulations_fixed_biases_newaniDom_10ind.csv",header=TRUE,sep=",")
db_rep <- read.table("databases_package/final_data_for_Figures_backup/Fig6_db_repeatability_100_simulations_fixed_biases_newaniDom_25ind.csv",header=TRUE,sep=",")
#db_rep <- read.table("databases_package/final_data_for_Figures_backup/Fig6_db_repeatability_100_simulations_fixed_biases_newaniDom_50ind.csv",header=TRUE,sep=",")


#db_rep.100 <- read.table("databases_package/final_data_for_Figures_backup/Fig6_db_repeatability_100_simulations_fixed_biases_newaniDom_10ind_100int.csv",header=TRUE,sep=",")
db_rep.100 <- read.table("databases_package/final_data_for_Figures_backup/Fig6_db_repeatability_100_simulations_fixed_biases_newaniDom_25ind_100int.csv",header=TRUE,sep=",")
#db_rep.100 <- read.table("databases_package/final_data_for_Figures_backup/Fig6_db_repeatability_100_simulations_fixed_biases_newaniDom_50ind_100int.csv",header=TRUE,sep=",")


db_rep.full <- rbind(db_rep,db_rep.100)

db_rep.full.2 <- db_rep.full[db_rep.full$poiss==1 &
                             #  db_rep.full$poiss==0 &
                               db_rep.full$dombias==0 &
                               db_rep.full$blevel==-5,]

db_rep.sum <-summaryBy(rep ~ Nobs, 
                       data = db_rep.full.2, 
                       FUN = function(x) { c(m = mean(x),
                                             q = quantile(x,probs=c(0.025,0.975))) })


names(db_rep.sum) <- c("Nobs",
                       "rep.m","lower","upper")

print(summary(db_rep.sum))


plot(db_rep.sum$rep.m~db_rep.sum$Nobs,0.5,type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(0.2,1))

axis(1,at=N.obs.values,
     cex.axis=1,tck=0.015,
     labels=FALSE)

axis(2,
     at=round(seq(0.2,1,0.1),1),
     cex.axis=1,las=2,tck=0.015)

points(db_rep.sum$Nobs,db_rep.sum$rep.m,type="b",col="blue",pch=19)
polygon(c(db_rep.sum$Nobs,rev(db_rep.sum$Nobs)),
        c(db_rep.sum$lower,rev(db_rep.sum$upper)),
        border=NA,col=rgb(0,0,1, 0.15))


text(98,0.25,"(b)",adj = 0 ,cex=1.5)

mtext("     Elo-rating repeatability",
      side=2, adj=0, line=4, cex=1.6)


# plotting the scenario

plot(c(1,25),
     #c(1,10),
     #c(1,50),
     c(0.5,1),type="n",
     ylab="", 
     xlab="",
     xaxt="n",
     yaxt="n",
     cex=1)

axis(1,
     at=seq(1,25,2),
     #at=seq(1,10,1),
     #at=seq(1,50,6),
     cex.axis=0.75,tck=0.015)

axis(2,at=seq(0.5,1,0.1),cex.axis=0.75,las=2,tck=0.015) 

plot_winner_prob(1:25,a=10,b=-5,"black")
#plot_winner_prob(1:10,a=10,b=-5,"black")
#plot_winner_prob(1:50,a=10,b=-5,"black")


mtext("P (dominant wins)",
      side=2, adj=0, line=3, cex=1.10)


plot(c(1,50),c(0,1),type="n",
     ylab="", 
     xlab="",
     xaxt="n",
     yaxt="n",
     frame.plot=FALSE)  

mtext("Difference in rank",
      side=3, adj=1, line=-2, cex=0.95)

mtext("  a = 10\n  b = -5",
      side=3, adj=0, line=-9, cex=1.5)




# HALVE COMPARISON PLOT

# db_split <- 
#   read.table("databases_package/final_data_for_Figures_backup/Fig7_db_split_elorand_100sim_fixed_biases_10ind.csv",header=TRUE,sep=",")
db_split <-
  read.table("databases_package/final_data_for_Figures_backup/Fig7_db_split_elorand_100sim_fixed_biases_25ind.csv",header=TRUE,sep=",")
# db_split <- 
#   read.table("databases_package/final_data_for_Figures_backup/Fig7_db_split_elorand_100sim_fixed_biases_50ind.csv",header=TRUE,sep=",")


#db_split.100 <- read.table("databases_package/final_data_for_Figures_backup/Fig7_db_split_elorand_100sim_fixed_biases_10ind_100int.csv",header=TRUE,sep=",")
db_split.100 <- read.table("databases_package/final_data_for_Figures_backup/Fig7_db_split_elorand_100sim_fixed_biases_25ind_100int.csv",header=TRUE,sep=",")
#db_split.100 <- read.table("databases_package/final_data_for_Figures_backup/Fig7_db_split_elorand_100sim_fixed_biases_50ind_100int.csv",header=TRUE,sep=",")


db_split.full <- rbind(db_split,db_split.100)

db_split.full.2 <- db_split.full[db_split.full$poiss==1 &
                                 #  db_split.full$poiss==0 &
                                   db_split.full$dombias==0 &
                                   db_split.full$blevel==-5,]


db_split.sum <-summaryBy(elo.rand.split ~ Nobs, 
                         data = db_split.full.2, 
                         FUN = function(x) { c(m = mean(x),
                                               q = quantile(x,probs=c(0.025,0.975))) })

names(db_split.sum) <- c("Nobs",
                 "elo.rand.split.m","lower","upper")

print(summary(db_split.sum))

plot(db_split.sum$elo.rand.split.m~db_split.sum$Nobs,0.5,type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(0,1))

axis(1,at=N.obs.values,
     labels=as.character(N.obs.values),
     cex.axis=0.75,tck=0.015)

mtext("    ratio of interactions to individuals",
      side=1, adj=0, line=4, cex=1.8) 

axis(2,
     at=seq(0,1,0.1),
     cex.axis=0.9,las=2)


#adding points for the means and shadowed areas for the 95% CI
points(db_split.sum$Nobs,db_split.sum$elo.rand.split.m,type="b",col="blue",pch=19)
polygon(c(db_split.sum$Nobs,rev(db_split.sum$Nobs)),
        c(db_split.sum$lower,rev(db_split.sum$upper)),
        border=NA,col=rgb(0,0,1, 0.15))

text(98,0.05,"(c)",adj = 0 ,cex=1.5)

mtext("  correlation coefficient (rs)",
      side=2, adj=0, line=4, cex=1.6)

# plotting the scenario

plot(c(1,25),
     #c(1,10),
     #c(1,50),
     c(0.5,1),type="n",
     ylab="", 
     xlab="",
     xaxt="n",
     yaxt="n",
     cex=1)

axis(1,
     at=seq(1,25,2),
     #at=seq(1,10,1),
     #at=seq(1,50,6),
     cex.axis=0.75,tck=0.015)

axis(2,at=seq(0.5,1,0.1),cex.axis=0.75,las=2,tck=0.015) 

plot_winner_prob(1:25,a=10,b=-5,"black")
#plot_winner_prob(1:10,a=10,b=-5,"black")
#plot_winner_prob(1:50,a=10,b=-5,"black")


mtext("P (dominant wins)",
      side=2, adj=0, line=3, cex=1.10)


plot(c(1,50),c(0,1),type="n",
     ylab="", 
     xlab="",
     xaxt="n",
     yaxt="n",
     frame.plot=FALSE)  

mtext("Difference in rank",
      side=3, adj=1, line=-2, cex=0.95)

mtext("  a = 10\n  b = -5",
      side=3, adj=0, line=-9, cex=1.5)



dev.off()
