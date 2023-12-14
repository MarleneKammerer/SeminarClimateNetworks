##### ERGM Script, CH
##### Authors: Marlene Kammerer

##### Date: 12.12.2023

setwd("C:/Users/kammerer/OneDrive - Universitaet Bern/Lehre/UniBE_HS23/EnvPol1_networks/SeminarClimateNetworks")

#set up packages

# install.packages("ergm")
# install.packages("statnet")
# install.packages("texreg")
# install.packages("network")
# install.packages("latticeExtra") # needed for GOF

library(ergm)
library(statnet)
library(texreg)
library(network)

seed <- 12345
set.seed (seed)

#Import data

OrgTypeCH <- read.csv2("Data/Actor_CH.csv", header=TRUE, fileEncoding = "latin1")
BeliefCH <- read.csv2("Data/Beliefs_CH.csv", header=TRUE, fileEncoding = "latin1")
BeliefCH$Name <- as.character(BeliefCH$Name)
RepCH <- read.csv2("Data/REP_CH.csv", header=TRUE, fileEncoding = "latin1")

coopCH <- as.matrix(read.table(file = "Data/COLLAB_CH.csv", 
                               header = TRUE, row.names = 1, sep = ";"))
RepCH[is.na(RepCH)] <- 0 
rownames(RepCH) <- RepCH$Name

RepCH <- subset (RepCH, RepCH$Name %in% BeliefCH$Name)
RepCH$Name <- NULL

list <- BeliefCH$Name
RepCH = RepCH[, (names(RepCH) %in% list)]

#Convert to matrices

RepCH <- as.matrix(RepCH)

#prepare orgtype data, so that it can be set as an attribute
BeliefCH$ActorType = as.character(BeliefCH$ActorType)
BeliefCH$ActorType[BeliefCH$ActorType== "1" ]<-"GOV"
BeliefCH$ActorType[BeliefCH$ActorType== "2" ]<-"SCI"
BeliefCH$ActorType[BeliefCH$ActorType== "3" ]<-"BUS"
BeliefCH$ActorType[BeliefCH$ActorType== "15" ]<-"GOV"
BeliefCH$ActorType[BeliefCH$ActorType== "4" ]<-"CIV"
BeliefCH$ActorType[BeliefCH$ActorType== "5" ]<-"CIV"

nw.coopCH <- network(coopCH, directed = TRUE, loops = FALSE)

set.vertex.attribute (nw.coopCH, "OrgType", BeliefCH$ActorType)

#create matrix of ties between bus actors and NGOs (expect a negative paramater estimate)

ch.bus.ngo <- matrix (0 , nrow = nrow (coopCH) , ncol = nrow (coopCH))
for (i in 1: nrow (ch.bus.ngo )) {
  for (j in 1: ncol ( ch.bus.ngo)) {
    if (( BeliefCH$ActorType [i] == "BUS" && BeliefCH$ActorType [j] == "CIV") ||
        (BeliefCH$ActorType [i] == "CIV" && BeliefCH$ActorType [j] == "BUS")) {
      ch.bus.ngo [i, j] <- 1
      ch.bus.ngo [j, i] <- 1
    }
  }
}

rownames(ch.bus.ngo) = rownames(coopCH)
colnames (ch.bus.ngo) = colnames (coopCH)
# 
# # #order of actor types for the nodefactor terms WHY??

u <- sort(unique(BeliefCH$ActorType)) 
u
nodecov <- match(BeliefCH$ActorType,u)
nodecov

# # #"BUS" "CIV" "GOV" "SCI"


#import and calculate beliefs distance matrix. 
BeliefCH$Name <- NULL
BeliefCH$ActorType <- NULL

Pol.dist.CH <- dist(BeliefCH, method = "manhattan")
PrefSimMatCH <- max(Pol.dist.CH) - Pol.dist.CH
PBsCH <-as.matrix(PrefSimMatCH)
rownames(PBsCH) <- list
colnames(PBsCH) <- list


#set Influence attribute

CHInf <- network(RepCH, directed = TRUE, loops = TRUE )
set.vertex.attribute (nw.coopCH, "Influence", degree (RepCH, cmode = "indegree"))

plot(nw.coopCH)

delete.vertices(nw.coopCH, which(degree(nw.coopCH )==0))

# Run ERGMs

# Bernoulli model
m1.CH <- ergm(nw.coopCH ~ edges,
              eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 2500 , MCMC.interval = 2500))
summary(m1.CH)

# Bernoulli model with mutuality
m2.CH <- ergm(coopCH ~ edges + mutual ,
              eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 2500 , MCMC.interval = 2500))
summary(m2.CH)

# Dependence model
m3.CH <- ergm(nw.coopCH ~ edges + mutual + gwidegree(1, fixed = FALSE ) + gwodegree(0, fixed = FALSE ) + twopath +
                gwesp(0 , fixed = FALSE ) ,
              eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 2500 , MCMC.interval = 2500))
summary(m3.CH)
plot(mcmc.diagnostics(m3.CH))

# Including beliefs
m4.CH <- ergm(nw.coopCH ~ edges + mutual + gwidegree(1, fixed = FALSE ) + gwodegree(0, fixed = FALSE ) + twopath +
                gwesp(0 , fixed = FALSE ) +  edgecov(PBsCH),
              eval.loglik = TRUE, check.degeneracy = TRUE, control = control.ergm ( seed = seed , MCMC.samplesize = 2500 , MCMC.interval = 2500))
summary(m4.CH)
plot(mcmc.diagnostics(m4.CH))


m6.CH <- ergm(nw.coopCH ~ edges + mutual + gwidegree(2, fixed = TRUE ) + gwodegree(0.5, fixed = TRUE ) + twopath +
                gwesp(1 , fixed = TRUE ) +  edgecov(PBsCH) + nodeifactor ("OrgType", base=c(1,2,4)) + 
                nodeofactor ("OrgType", base=c(1,3,4)) + nodematch ("OrgType") +
                edgecov (ch.bus.ngo) + edgecov (RepCH) + nodeicov("Influence")  + 
                absdiff ("Influence") 
                , eval.loglik = TRUE, check.degeneracy = TRUE, 
              control = control.ergm ( seed = seed , MCMC.samplesize = 5000 , MCMC.interval = 5000))
summary(m6.CH)
plot(mcmc.diagnostics(m6.CH))


#---------------------------------------------
screenreg(list(m1.CH, m2.CH, m3.CH, m4.CH, m5.CH, m6.CH))

####### Plot Network

library(sna)
# install.packages("GGally")
library(GGally)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("RColorBrewer")
library(RColorBrewer)
# install.packages("intergraph")
library(intergraph)
# install.packages("scales")
library(scales)


ggnet2(nw.coopCH, alpha = 0.5, label = FALSE, label.size = 3, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.5),  node.size = 6, 
       color = "OrgType",
       color.legend = "Actor Type", edge.color = "black", 
       palette = c("CIV" = "green", "GOV" = "tomato", "BUS" = "grey", "SCI" = "blue"), 
       size = "indegree", size.min = 1, size.cut = 10, size.legend = "Indegree centrality",
       legend.position = "bottom")
ggsave("CH59_graph.png", width = 29, height = 18,  device = NULL, dpi = 300)    

# Vertex Centrality

# Function for centrality stats

detach(package:sna)

install.packages("igraph")
library(igraph)

myCentrality <- function(net) {
  
  if (!is.igraph(net)) stop ("Input is not an igraph object")
  ideg <- degree (net, mode="in", loop=F, normalized = T)
  odeg <- degree (net, mode="out", loop=F, normalized = T)
  btw <- betweenness (net, normalized=TRUE)
  clo <- closeness (net, normalized=TRUE) # consider removing
  evc <- evcent(net)$vector
  pgr <- page.rank(net)$vector
  id <- V(net)
  
  ret <-data.frame (cbind(ID=id, IDegree=ideg, ODegree=odeg,
                           Betweenness=btw, Closeness=clo,
                          Eigenvector=evc, PageRank=pgr))
  return(ret)
  
}

# make network graph for igraph

nw.coopFin.igraph <- graph.adjacency(coopFin)
Results1 <- myCentrality(nw.coopFin.igraph)

# Rounding 

Results1$Betweenness <- round(Results1$Betweenness, digits =3)
Results1$Closeness <- round(Results1$Closeness, digits =3)
Results1$Eigenvector <- round(Results1$Eigenvector, digits =3)
Results1$PageRank <- round(Results1$PageRank, digits =3)
Results1$IDegree <- round(Results1$IDegree, digits =3)
Results1$ODegree <- round(Results1$ODegree, digits =3)

Results1

write.csv2(Results1, file = "centrality.scores.Fin.csv")

Results1$Actor.type <- OrgTypeFin

## Calculate means per actor type
Results1$ID <- NULL
Results1$Actor.type <- as.character(Results1$Actor.type)

mean_by_actor <- aggregate(Results1, by=list(Results1$Actor.type), FUN=mean, na.rm=TRUE )
mean_by_actor

mean_by_actor$Actor.type <- NULL

mean_by_actor$Betweenness <- round(mean_by_actor$Betweenness, digits =3)
mean_by_actor$Closeness <- round(mean_by_actor$Closeness, digits =3)
mean_by_actor$Eigenvector <- round(mean_by_actor$Eigenvector, digits =3)
mean_by_actor$PageRank <- round(mean_by_actor$PageRank, digits =3)
mean_by_actor$IDegree <- round(mean_by_actor$IDegree, digits =3)
mean_by_actor$ODegree <- round(mean_by_actor$ODegree, digits =3)

write.csv2(mean_by_actor, file = "centrality.scores-mean-by-actor.csv")

