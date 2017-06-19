
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

# more info at the Open Science Framework: http://doi.org/10.17605/OSF.IO/9GYEK


###############################################################################
# Packages needed
###############################################################################

# packages needed for this script

library(aniDom)
library(doBy)
library(RColorBrewer)
library(igraph)
library(EloRating)


# Clear memory
rm(list=ls())


avalues <- c(0)
bvalues <- c(-5)
N.inds.values <- c(5)
N.obs.values <- c(10)
poiss <- c(TRUE)
dombias <- c(FALSE)

output <- generate_interactions(N.inds.values,
                                N.inds.values*N.obs.values,
                                a=avalues,
                                b=bvalues,
                                id.biased=poiss,
                                rank.biased=dombias)

int.matrix <- output$interactions

int.matrix$Winner <- ifelse(int.matrix$Winner==1,"A",
                            ifelse(int.matrix$Winner==2,"B",
                                   ifelse(int.matrix$Winner==3,"C",
                                          ifelse(int.matrix$Winner==4,"D","E"
                                          ))))

int.matrix$Loser <- ifelse(int.matrix$Loser==1,"A",
                            ifelse(int.matrix$Loser==2,"B",
                                   ifelse(int.matrix$Loser==3,"C",
                                          ifelse(int.matrix$Loser==4,"D","E"
                                          ))))

int.matrix$Loser[50]<-"E"

write.csv(int.matrix,
          "plots/Fig1_Fake_interaction_network.csv",row.names=FALSE)



# Plotting

data <- read.csv("plots/Fig1_Fake_interaction_network.csv", stringsAsFactors=FALSE)

ids <- unique(c(data[,1],data[,2]))

network <- matrix(0,nrow=length(ids),ncol=length(ids))

for (i in 1:length(ids)) {
  for (j in 1:length(ids)) {
    network[i,j] <- sum(data[,1]==ids[i] & data[,2]==ids[j])
  }
}
rownames(network) <- colnames(network) <- ids

write.csv(network,
          "plots/Fig1_Fake_interaction_network_matrix.csv",row.names=FALSE)


net <- graph.adjacency(network, weighted=TRUE)


tiff("plots/Figure1_Fake_network.tiff",
     height=21, width=29.7,
     units='cm', compression="lzw", res=600)


#ll <- layout_(net, as_star())

V(net)$label.cex <- c(3)

plot(net, vertex.color="chartreuse4", 
     vertex.label.color="black",
     edge.width=E(net)$weight/4,
     edge.arrow.width=1.75,
     vertex.size=30,
     vertex.label.dist=0,
     edge.arrow.size=1.25,
     edge.curved=FALSE,
     edge.color="grey25",
     #layout=ll
     vertex.label = V(net)$name)

dev.off()



result.no.rand <- elo_scores(data$Winner,
                             data$Loser,
                             identities=c("A","B","C","D","E"),
                             init.score=1000,
                             randomise=FALSE,
                             return.as.ranks=FALSE)

w<-scale.elo(result.no.rand)
