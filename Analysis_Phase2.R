#===================================================================================================
### BASIC ANALYSIS  PHASE II
### Author: Marlene Kammerer
### Date: 06.07.2017, revised 06.09.2017
#===================================================================================================

# Set up

detach(package:sna)
library(igraph)
igraph.options(sparsematrices=FALSE)

# Load data

Net <- read.csv2("final-mat_2013-for-bookchapter.csv", header = TRUE)
Covs <- read.csv2("covariates_2013-for-bookchapter.csv", header = TRUE)
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

# Rounding 

Results1$Betweenness <- round(Results1$Betweenness, digits =3)
Results1$Closeness <- round(Results1$Closeness, digits =3)
Results1$Eigenvector <- round(Results1$Eigenvector, digits =3)
Results1$PageRank <- round(Results1$PageRank, digits =3)
Results1$IDegree <- round(Results1$IDegree, digits =3)
Results1$ODegree <- round(Results1$ODegree, digits =3)

Results1

write.csv2(Results1, file = "centrality.scores-adjusted.dichotomized.csv")

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


net <- network(Net_mat, directed = TRUE)
net

# attach attributes

## it's important that attributes are transformed into charachters first

Covs$ActorType <- as.character(Covs$ActorType)
Covs$targets <- as.character(Covs$targets)
Covs$tax <- as.character(Covs$tax)
Covs$scope <- as.character(Covs$scope)


# now add to network; could also be done with set.vertex.attribute()

net%v%"actorType" <-Covs$ActorType
net%v%"targets" <-Covs$targets
net%v%"scope" <- Covs$scope
net%v%"tax" <-Covs$tax

list.vertex.attributes(net) # check if attributes are loaded 

# I create my plots

dev.off()

ggnet2(net, mode = "fruchtermanreingold",  node.size = 6, node.color = "steelblue", edge.size = 0.5,
       edge.color = "black")

# I test different layouts and decide for fruchtermanreingold
ggnet2(net_w, "kamadakawai" )
ggnet2(net_w, "target")
ggnet2(net_w, "fruchtermanreingold")

# my graphs

palette = c("Administration" = "#252525", "Legislative Branch" = "#d9d9d9", 
            "Parties" = "#bdbdbd", "Citizen Group" ="#969696", "Private Sector" = "#737373", "Science" = "#525252", 
            "International" = "#f7f7f7" )

ggnet2(net, alpha = 0.75, label = TRUE, label.size = 8, mode = "kamadakawai",  node.size = 6, color = "actorType", edge.size = 0.4, 
       color.legend = "Actor Type", edge.color = "grey", 
       palette = palette, 
       size = "degree", size.min = 1, size.cut = 5, size.legend = "Degree centrality",
       shape = "tax", shape.legend = "CO2 Tax",
       shape.palette = c("support" = 16, "oppose" = 17, "no" = 3),
       legend.size = 20, legend.position = "bottom")
ggsave("new_13.2.jpg", width = 25, height = 18,  device = NULL, dpi = 400)


# ggnet2(net, alpha = 0.75, label = TRUE, label.size = 5, mode = "kamadakawai",  node.size = 15, color = "actorType", edge.size = 0.4, 
#        color.legend = "Actor Type", edge.color = "grey", 
#        palette = c("Administration" = "grey", "Legislative Branch" = "tomato", 
#                    "Parties" = "orange", "Citizen Group" ="green", "Private Sector" = "lightblue", "Science" = "yellow", "International" = "pink" ), 
#        size = "degree", size.min = 1, size.cut = 5, size.legend = "Degree centrality",
#        shape = "CO2Levy", shape.legend = "CO2 Levy",
#        shape.palette = c("support" = 16, "oppose" = 17, "no" = 3),
#        legend.size = 10, legend.position = "bottom")
# ggsave("graph_degree_actor_CO2levy_final.png", width = 29, height = 18,  device = NULL, dpi = 300)
# 
# 
# ggnet2(net, alpha = 0.75, label = TRUE, label.size = 3,
#        mode = "kamadakawai",  node.size = 6, 
#        color = "targets", edge.size = 0.4, 
#        color.legend = "Reduction targets", edge.color = "grey", 
#        palette = c("high" = "green", "medium" = "yellow", "soft" = "tomato", "no" = "grey"), 
#        size = "degree", size.min = 1, size.cut = 5, size.legend = "Degree centrality",
#        legend.position = "bottom")
# 
# ggnet2(net, alpha = 0.75, label = TRUE, label.size = 3,
#        mode = "kamadakawai",  node.size = 6, 
#        color = "tax", edge.size = 0.4, 
#        color.legend = "Carbon levy on fuels", edge.color = "grey", 
#        palette = c("support" = "green", "oppose" = "tomato", "no" = "grey"), 
#        size = "degree", size.min = 1, size.cut = 5, size.legend = "Degree centrality",
#        legend.position = "bottom")
# 
# ggnet2(net, alpha = 0.75, label = TRUE, label.size = 3,
#        mode = "kamadakawai",  node.size = 6, 
#        color = "scope", edge.size = 0.4, 
#        color.legend = "Flexibiity of reduction", edge.color = "grey", 
#        palette = c("domestic" = "green", "flex" = "tomato", "no" = "grey"), 
#        size = "degree",size.min = 1, size.cut = 5, size.legend = "Degree centrality",
#        legend.position = "bottom")
# 
# # Now I combine Carbon tax graph with positins on instruments
# 
# ggnet2(net, alpha = 0.75, label = TRUE, label.size = 5, mode = "kamadakawai",  node.size = 15, color = "actorType", edge.size = 0.4, 
#        color.legend = "Actor Type", edge.color = "grey", 
#        palette = c("Administration" = "grey", "Legislative Branch" = "tomato", 
#                    "Parties" = "orange", "Citizen Group" ="green", "Private Sector" = "lightblue", "Science" = "yellow", "International" = "pink" ), 
#        size = "degree", size.min = 1, size.cut = 5, size.legend = "Degree centrality",
#        shape = "voluntary measures", shape.legend = "Voluntary measures",
#        shape.palette = c("support" = 16, "oppose" = 17, "no" = 3),
#        legend.size = 10, legend.position = "bottom")
# ggsave("graph_degree_actor_va_final.png", width = 29, height = 18,  device = NULL, dpi = 300)
# 
# ggnet2(net, alpha = 0.75, label = TRUE, label.size = 5, mode = "kamadakawai",  node.size = 15, color = "actorType", edge.size = 0.4, 
#        color.legend = "Actor Type", edge.color = "grey", 
#        palette = c("Administration" = "grey", "Legislative Branch" = "tomato", 
#                    "Parties" = "orange", "Citizen Group" ="green", "Private Sector" = "lightblue", "Science" = "yellow", "International" = "pink" ), 
#        size = "degree", size.min = 1, size.cut = 5, size.legend = "Degree centrality",
#        shape = "CO2Levy", shape.legend = "CO2 Levy",
#        shape.palette = c("support" = 16, "oppose" = 17, "no" = 3),
#        legend.size = 10, legend.position = "bottom")
# ggsave("graph_degree_actor_CO2levy_final.png", width = 29, height = 18,  device = NULL, dpi = 300)
# 
# ggnet2(net, alpha = 0.75, label = TRUE, label.size = 6, mode = "kamadakawai",  node.size = 20, color = "actorType", edge.size = 0.4, 
#        color.legend = "Actor Type", edge.color = "grey", 
#        palette = c("Administration" = "grey", "Legislative Branch" = "tomato", 
#                    "Parties" = "orange", "Citizen Group" ="green", "Private Sector" = "lightblue", "Science" = "yellow", "International" = "pink" ), 
#        size = "degree", size.min = 1, size.cut = 5, size.legend = "Degree centrality",
#        shape = "permits", shape.legend = "Permits",
#        shape.palette = c("support" = 16, "oppose" = 17, "no" = 3),
#        legend.size = 15, legend.position = "bottom")
# ggsave("graph_degree_actor_permits_final.png", width = 29, height = 18,  device = NULL, dpi = 300)
# 
# 
# ggnet2(net, alpha = 0.75, label = TRUE, label.size = 6, mode = "kamadakawai",  node.size = 20, color = "actorType", edge.size = 0.4, 
#        color.legend = "Actor Type", edge.color = "grey", 
#        palette = c("Administration" = "grey", "Legislative Branch" = "tomato", 
#                    "Parties" = "orange", "Citizen Group" ="green", "Private Sector" = "lightblue", "Science" = "yellow", "International" = "pink" ), 
#        size = "degree", size.min = 1, size.cut = 5, size.legend = "Degree centrality",
#        shape = "climatecent", shape.legend = "Climate Cent",
#        shape.palette = c("support" = 16, "oppose" = 17, "no" = 3),
#        legend.size = 15, legend.position = "bottom")
# ggsave("graph_degree_actor_climatecent_final.png", width = 29, height = 18,  device = NULL, dpi = 300)


# Average degree

density <- graph.density(graph = Net_binar)

avdegree <- density*((vcount(Net_binar)-1))

# Centralization

centralize(degree(Net_binar), normalized = TRUE) # Need to check this

# Run simple regression


## Recoding belief variables

levels(Covs$CO2Levy)[levels(Covs$CO2Levy)=="no"] <- 3
levels(Covs$CO2Levy)[levels(Covs$CO2Levy)=="oppose"] <- 2
levels(Covs$CO2Levy)[levels(Covs$CO2Levy)=="support"] <- 1

levels(Covs$Permits)[levels(Covs$Permits)=="no"] <- 3
levels(Covs$Permits)[levels(Covs$Permits)=="oppose"] <- 2
levels(Covs$Permits)[levels(Covs$Permits)=="support"] <- 1

levels(Covs$voluntaryMeasures)[levels(Covs$voluntaryMeasures)=="no"] <- 0
levels(Covs$voluntaryMeasures)[levels(Covs$voluntaryMeasures)=="oppose"] <- -1
levels(Covs$voluntaryMeasures)[levels(Covs$voluntaryMeasures)=="support"] <- 1

levels(Covs$ClimateCent)[levels(Covs$ClimateCent)=="no"] <- 0
levels(Covs$ClimateCent)[levels(Covs$ClimateCent)=="oppose"] <- -1
levels(Covs$ClimateCent)[levels(Covs$ClimateCent)=="support"] <- 1

Covs$Permits = as.numeric(Covs$Permits)
Covs$voluntaryMeasures = as.numeric(Covs$voluntaryMeasures)
Covs$ClimateCent = as.numeric(Covs$ClimateCent)
Covs$CO2Levy = as.numeric(Covs$CO2Levy)


data = cbind(Covs, Results1)

# I control for different centrality measures

library(texreg)
library(car)

levels(data$ActorType)

fit.1a <- lm (CO2Levy ~ IDegree + C(ActorType, contr.treatment(6, base=1)) 
                               , data = data)
summary(fit.1a) 

fit.1b <- lm (CO2Levy~ Betweenness + C(ActorType, contr.treatment(6, base=1)) 
               , data = data)
summary(fit.1b)  

fit.1c <- lm (CO2Levy ~ Closeness+ C(ActorType, contr.treatment(6, base=1)) 
               , data = data)
summary(fit.1c)                    
                   

screenreg(list(fit.1a, fit.1b, fit.1c))

# The different centrality measures, do not really make a difference

####################################################################

# I control for other instruments


fit.2a <- lm (CO2Levy ~ IDegree + C(ActorType, contr.treatment(6, base=1)) + 
                 
               Permits + ClimateCent + voluntaryMeasures , data = data)
summary(fit.2a) 
vif(fit.2a)


screenreg(list(fit.1a, fit.2a))

####################################################################


fit.3a <- lm(cent.scores$IDegree ~ 
               Covs$ClimateCent + Covs$ActorType)
summary(fit.3a)


fit.3b <- lm(cent.scores$Betweenness ~ 
               Covs$ClimateCent + Covs$ActorType)
summary(fit.3b)


fit.3c <- lm(cent.scores$Closeness ~ 
               Covs$ClimateCent + Covs$ActorType)
summary(fit.3c)

screenreg(list(fit.3a, fit.3b, fit.3c))

####################################################################


fit.4a <- lm(cent.scores$IDegree ~ 
               Covs$Permits + Covs$ActorType)
summary(fit.4a)


fit.4b <- lm(cent.scores$Betweenness ~ 
               Covs$Permits + Covs$ActorType)
summary(fit.4b)


fit.4c <- lm(cent.scores$Closeness ~ 
               Covs$Permits + Covs$ActorType)
summary(fit.4c)

screenreg(list(fit.4a, fit.4b, fit.4c))