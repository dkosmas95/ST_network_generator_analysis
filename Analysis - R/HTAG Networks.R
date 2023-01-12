library(readr)
library(stringr)
library(dplyr)
library(tidyverse)
library(prodlim)
library(pracma)
library(igraph)
library(graphlayouts)
library(doBy)
library(partitions)

gAG1 <- make_empty_graph(n=4, directed=FALSE)
V(gAG1)$role <- c("trafficker","bottom","victim","victim")
gAG1 <- add.edges(gAG1, c(1,2,1,4,2,3,2,4,3,4))

V(gAG1)$color <- c("lightgoldenrod1","lightgreen","lightskyblue","lightskyblue")

gAG2 <- make_empty_graph(n=17, directed=FALSE)
V(gAG2)$role <- c("trafficker",rep("victim",5), "trafficker", "bottom", "bottom", rep("victim",3), "trafficker", "bottom", rep("victim",3))
gAG2 <- add.edges(gAG2, c(1,2,1,3,1,4,1,5,1,6)) # operation 1
gAG2 <- add.edges(gAG2, c(7,8,7,9,8,10,8,11,8,12,9,12)) #operation 2
gAG2 <- add.edges(gAG2, c(13,14,13,15,13,16,14,16,14,17)) #operation 3
gAG2 <- add.edges(gAG2, c(1,7,7,13)) #trafficker social arcs
gAG2 <- add.edges(gAG2, c(6,12,10,11,11,12,10,14,10,17,16,17)) #victim social arcs

V(gAG2)$color <- c("lightgoldenrod1",rep("lightskyblue",5), "lightgoldenrod1", "lightgreen", "lightgreen", rep("lightskyblue",3),"lightgoldenrod1", "lightgreen", rep("lightskyblue",3))

vicAG1 <- induced.subgraph(gAG1, which(V(gAG1)$role=="victim"))
vicAG2 <- induced.subgraph(gAG2, which(V(gAG2)$role=="victim")) # full second network

subAG2_1 <- induced.subgraph(gAG2, 1:6) # second network, operation 1
subAG2_2 <- induced.subgraph(gAG2, 7:12) # second network, operation 2
subAG2_3 <- induced.subgraph(gAG2, c(10,13:17)) # second network, operation 3 (includes overlapping victim)

subvicAG2_1 <- induced.subgraph(subAG2_1, which(V(subAG2_1)$role=="victim"))
subvicAG2_2 <- induced.subgraph(subAG2_2, which(V(subAG2_2)$role=="victim")) 
subvicAG2_3 <- induced.subgraph(subAG2_3, which(V(subAG2_3)$role=="victim"))


#### Victim social network centrality
#AG1
ag1_Soc_ed <- edge_density(vicAG1)
ag1_Soc_degc <- centr_degree(vicAG1)$centralization
ag1_Soc_betc <- centr_betw(vicAG1)$centralization 
ag1_Soc_cloc <- centr_clo(vicAG1)$centralization 
ag1_Soc_eigc <- centr_eigen(vicAG1)$centralization

#AG2
ag2_Soc_ed <- edge_density(vicAG2)
ag2_Soc_degc <- centr_degree(vicAG2)$centralization
ag2_Soc_betc <- centr_betw(vicAG2)$centralization 
ag2_Soc_cloc <- centr_clo(vicAG2)$centralization 
ag2_Soc_eigc <- centr_eigen(vicAG2)$centralization

#sub AG2-1
sag21_Soc_ed <- edge_density(subvicAG2_1)
sag21_Soc_degc <- centr_degree(subvicAG2_1)$centralization
sag21_Soc_betc <- centr_betw(subvicAG2_1)$centralization 
sag21_Soc_cloc <- centr_clo(subvicAG2_1)$centralization 
sag21_Soc_eigc <- centr_eigen(subvicAG2_1)$centralization

#sub AG2-2
sag22_Soc_ed <- edge_density(subvicAG2_2)
sag22_Soc_degc <- centr_degree(subvicAG2_2)$centralization
sag22_Soc_betc <- centr_betw(subvicAG2_2)$centralization 
sag22_Soc_cloc <- centr_clo(subvicAG2_2)$centralization 
sag22_Soc_eigc <- centr_eigen(subvicAG2_2)$centralization

#sub AG2-3
sag23_Soc_ed <- edge_density(subvicAG2_3)
sag23_Soc_degc <- centr_degree(subvicAG2_3)$centralization
sag23_Soc_betc <- centr_betw(subvicAG2_3)$centralization 
sag23_Soc_cloc <- centr_clo(subvicAG2_3)$centralization 
sag23_Soc_eigc <- centr_eigen(subvicAG2_3)$centralization


##### eigen comparison
## centrality
ag1_Soc_eig <- eigen_centrality(vicAG1)
sag21_Soc_eig <- eigen_centrality(subvicAG2_1)
sag22_Soc_eig <- eigen_centrality(subvicAG2_2)
sag23_Soc_eig <- eigen_centrality(subvicAG2_3)

## spectrum
ag1mat <- as_adjacency_matrix(vicAG1,sparse = FALSE)
ag1_spec <- eigen(ag1mat)
sag21mat <- as_adjacency_matrix(subvicAG2_1,sparse = FALSE)
sag21_spec <- eigen(sag21mat)
sag22mat <- as_adjacency_matrix(subvicAG2_2,sparse = FALSE)
sag22_spec <- eigen(sag22mat)
sag23mat <- as_adjacency_matrix(subvicAG2_3,sparse = FALSE)
sag23_spec <- eigen(sag23mat)

ag1lapmat <- laplacian_matrix(vicAG1, sparse=FALSE)
ag1lap_spec <- eigen(ag1lapmat)
sag21lapmat <- laplacian_matrix(subvicAG2_1, sparse=FALSE)
sag21lap_spec <- eigen(sag21lapmat)
sag22lapmat <- laplacian_matrix(subvicAG2_2, sparse=FALSE)
sag22lap_spec <- eigen(sag22lapmat)
sag23lapmat <- laplacian_matrix(subvicAG2_3, sparse=FALSE)
sag23lap_spec <- eigen(sag23lapmat)


#### comparison with generated operations
gGen <- read.graph(file="/HT5_1_40.graphml", format="graphml")
vAug <- which(grepl("victim", V(gGen)$role))  # identifies all victim nodes
tList <- which(grepl("trafficker", V(gGen)$role))
bList <- tList <- which(grepl("bottom", V(gGen)$role))
vByTraff <- list()
vByTraff[[1]] <- 2:4 # correct for appropriate data set
vByTraff[[2]] <- 6:8 # correct for appropriate data set
vByTraff[[3]] <- 11:16 # correct for appropriate data set
vByTraff[[4]] <- 19:21 # correct for appropriate data set
vByTraff[[5]] <- 24:28 # correct for appropriate data set

gOp <- list()
gOpMat <- list()
gOpLapMat <- list()
gOp_spec <- list()
gOpLap_spec <- list()
for (i in 1:5){
  gOp[[i]] <- induced.subgraph(gGen, vByTraff[[i]])
  gOpMat[[i]] <- as_adjacency_matrix(gOp[[i]])
  gOp_spec[[i]] <- eigen(gOpMat[[i]])
  gOpLapMat[[i]] <- laplacian_matrix(gOp[[i]], sparse=FALSE)
  gOpLap_spec[[i]] <- eigen(gOpLapMat[[i]])
}


## distance between spectra 
spec1 <- gOpLap_spec[[1]]$values # operation 1 for comparision, update as needed
spec2 <- gOpLap_spec[[2]]$values # operation 2 for comparision, update as needed

longer <- max(c(length(spec1),length(spec2)))

spec1 <- sort(c(spec1, rep(0, longer-length(spec1))))
spec2 <- sort(c(spec2, rep(0, longer-length(spec2))))
Norm(spec1-spec2)

