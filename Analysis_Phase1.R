#===================================================================================================
### BASIC ANALYSIS PHASE I
### Author: Marlene Kammerer
### Date: 06.07.2017, revised 06.09.2017, updated 19.11.2023
#===================================================================================================

# Set up

# detach(package:sna)

library(igraph)

igraph.options(sparsematrices=FALSE)

# Load data

Net <- read.csv2("Data/Network_Phase1.csv", header = TRUE)
Covs <- read.csv2("Data/Covariates_Phase1.csv", header = TRUE)

rownames(Covs) <- Covs$Abbreviation
rownames(Net) <- Net$ID
Net$ID <- NULL
Net_mat <- as.matrix(Net)

# Dichotomized network (done in UCINET, cut-off point 6, correlation around 0.7)

Net_binar <- graph.adjacency(Net_mat, mode="directed", diag=FALSE)

str(Net_binar)
plot(Net_binar)

# Calculate basic statistics
# Vertex degree

dev.off()

hist (degree (Net_binar, mode = "out"), col="blue", xlim=c(0,40),
      xlab = "Vertex Outdegree", ylab="Frequency", main= "Degree distribution Implementation")

# Vertex Centrality

# Function for centrality stats

myCentrality <- function(net) {
  
  if (!is.igraph(net)) stop ("Input is not an igraph object")
  ideg <- degree (net, mode="in", loop=F, normalized = T)
  odeg <- degree (net, mode="out", loop=F, normalized = T)
  istr <- strength(net, mode = "in", loop=F)
  ostr <- strength(net, mode = "out", loop=F)
  btw <- betweenness (net, normalized=TRUE)
  clo <- closeness (net, normalized=TRUE) # consider removing
  evc <- evcent(net)$vector
  pgr <- page.rank(net)$vector
  id <- V(net)
  
  ret <-data.frame (cbind(ID=id, IDegree=ideg, ODegree=odeg, IStrength=istr,
                          OStrength=ostr, Betweenness=btw, Closeness=clo,
                          Eigenvector=evc, PageRank=pgr))
  return(ret)
  
}

Results1 <- myCentrality(Net_binar)
Results1


# Rounding 

Results1$Betweenness <- round(Results1$Betweenness, digits =3)
Results1$Closeness <- round(Results1$Closeness, digits =3)
Results1$Eigenvector <- round(Results1$Eigenvector, digits =3)
Results1$PageRank <- round(Results1$PageRank, digits =3)
Results1$IDegree <- round(Results1$IDegree, digits =3)
Results1$ODegree <- round(Results1$ODegree, digits =3)

Results1

write.csv2(Results1, file = "Results/centrality.scores-adjusted.dichotomized.csv")

Results1$Actor.type <- Covs$ActorType

## Calculate means per actor type
Results1$ID <- NULL
Results1$Actor.type <- as.character(Results1$Actor.type)

mean_by_actor <- aggregate(Results1, by=list(Results1$Actor.type), FUN=mean, na.rm=TRUE )
mean_by_actor

mean_by_actor$Actor.type <- NULL

write.csv2(mean_by_actor, file = "centrality.scores-mean-by-actor.csv")

# Plots

detach(package:igraph)
library(sna)
library(GGally)
library(network)
library(ggplot2)
library(RColorBrewer)
library(intergraph)
library(scales)


net <- network(Net_mat, directed = TRUE) # create network w/ network package
net

# attach attributes

## it's important that attributes are transformed into characters first
Covs$ActorType <- as.character(Covs$ActorType)
Covs$voluntaryMeasures <- as.character(Covs$voluntaryMeasures)
Covs$CO2Levy <- as.character(Covs$CO2Levy)
Covs$Permits <- as.character(Covs$Permits)
Covs$ClimateCent <- as.character(Covs$ClimateCent)

# now add to network; could also be done with set.vertex.attribute()
net%v%"actorType" <-Covs$ActorType
net%v%"voluntary measures" <-Covs$voluntaryMeasures
net%v%"CO2Levy" <- Covs$CO2Levy
net%v%"permits" <-Covs$Permits
net%v%"climatecent" <- Covs$ClimateCent

list.vertex.attributes(net) # check if attributes are loaded 

# I create my plots

dev.off()

ggnet2(net, mode = "fruchtermanreingold",  node.size = 6, node.color = "steelblue", edge.size = 0.5,
       edge.color = "black")

# I test different layouts and decide for fruchtermanreingold
ggnet2(net, "kamadakawai" )
ggnet2(net, "target")
ggnet2(net, "fruchtermanreingold")

# my graphs

palette = c("Administration" = "#252525", "Legislative Branch" = "#d9d9d9", 
            "Parties" = "#bdbdbd", "Citizen Group" ="#969696", "Private Sector" = "#737373", "Science" = "#525252", 
            "International" = "#f7f7f7" )

ggnet2(net, alpha = 0.75, label = TRUE, label.size = 3, mode = "kamadakawai",  node.size = 6, color = "actorType", edge.size = 0.4, 
       color.legend = "Actor Type", edge.color = "grey", 
       palette = palette, 
       size = "degree", size.min = 1, size.cut = 5, size.legend = "Degree centrality",
       shape = "CO2Levy", shape.legend = "CO2 Tax",
       shape.palette = c("support" = 16, "oppose" = 17, "no" = 3),
       legend.size = 10, legend.position = "bottom")
ggsave("Results/graph_degree_actor_tax_2005.png", width = 29, height = 18,  device = NULL, dpi = 300)

# Run simple regression

## Recoding belief variables

library(texreg)
data = cbind(Covs, Results1)


fit.1a <- lm(data$IDegree ~ 
               data$ClimateCent + data$ActorType)
summary(fit.1a)


fit.1b <- lm(data$Betweenness ~ 
               data$ClimateCent + data$ActorType)
summary(fit.1b)


fit.1c <- lm(data$Closeness ~ 
               data$ClimateCent + Covs$ActorType)
summary(fit.1c)

screenreg(list(fit.1a, fit.1b, fit.1c))

####################################################################


fit.4a <- lm(data$IDegre ~ 
               data$Permits + data$ActorType)
summary(fit.4a)


fit.4b <- lm(data$Betweenness ~ 
               data$Permits + data$ActorType)
summary(fit.4b)


fit.4c <- lm(data$Closeness ~ 
               data$Permits + data$ActorType)
summary(fit.4c)

screenreg(list(fit.4a, fit.4b, fit.4c))
