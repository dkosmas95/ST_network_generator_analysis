library(readr)
library(stringr)
library(dplyr)
library(tidyverse)
library(prodlim)
library(pracma)
library(igraph)
library(doBy)

#### Operations Central
centralV <- make_empty_graph(n=4, directed = FALSE)
centralV_edges <- c(1,2,1,3,1,4,2,3,2,4,3,4)
centralV <- add.edges(centralV, centralV_edges)
plot(centralV)
cen_ed <- edge_density(centralV) # edge density 1
cen_degc <- centr_degree(centralV)$centralization # degree centrality 0
cen_betc <- centr_betw(centralV)$centralization # betweenness centrality 0
cen_cloc <- centr_clo(centralV)$centralization # closeness centrality 0 
cen_eigc <- centr_eigen(centralV)$centralization # eigen centraility 0

centralmat <- as_adjacency_matrix(centralV,sparse = FALSE)
central_spec <- eigen(centralmat)
centrallapmat <- laplacian_matrix(centralV, sparse=FALSE)
centrallap_spec <- eigen(centrallapmat)

#### Operation Chalice
chaliceV <- make_empty_graph(n=6, directed = FALSE)
chaliceV_edges <- c(1,2,1,3,1,5,1,6,2,3,2,4)
chaliceV <- add.edges(chaliceV, chaliceV_edges)
plot(chaliceV)
chal_ed <- edge_density(chaliceV) # edge density 0.4
chal_degc <- centr_degree(chaliceV)$centralization # degree centrality 0.4
chal_betc <- centr_betw(chaliceV)$centralization # betweenness centrality 0.62
chal_cloc <- centr_clo(chaliceV)$centralization # closeness centrality 0.617776
chal_eigc <- centr_eigen(chaliceV)$centralization # eigen centraility 0.5500983

chalicemat <- as_adjacency_matrix(chaliceV,sparse = FALSE)
chalice_spec <- eigen(chalicemat)
chalicelapmat <- laplacian_matrix(chaliceV, sparse=FALSE)
chalicelap_spec <- eigen(chalicelapmat)

#### Operation Retriever
retrieverV <- make_empty_graph(n=25, directed = FALSE)
retrieverV_edges <- c(1,11,1,4,1,2,1,9,1,3,1,6,1,24,1,7,1,13,2,9,2,7,2,8,2,14,3,9,3,10,3,24,3,6,4,11,4,5,4,16,4,8,4,9,5,16,5,12,5,17,5,11,6,24,7,9,7,8,9,10,10,15,11,25,12,17,18,19)
retrieverV <- add.edges(retrieverV, retrieverV_edges)
plot(retrieverV)
r_ed <- edge_density(retrieverV) # edge density 0.1133333
r_degc <- centr_degree(retrieverV)$centralization # degree centrality 0.2616667
r_betc <- centr_betw(retrieverV)$centralization # betweenness centrality 0.20003573 
r_cloc <- centr_clo(retrieverV)$centralization # closeness centrality 0.06308256  (not well defined, disconnected)
r_eigc <- centr_eigen(retrieverV)$centralization # eigen centraility 0.7647236

retrievermat <- as_adjacency_matrix(retrieverV,sparse = FALSE)
retriever_spec <- eigen(retrievermat)
retrieverlapmat <- laplacian_matrix(retrieverV, sparse=FALSE)
retrieverlap_spec <- eigen(retrieverlapmat)


retrieverV_connect <- decompose(retrieverV)[[1]] 
plot(retrieverV_connect)
rc_ed <- edge_density(retrieverV_connect) # edge density 0.4
rc_degc <- centr_degree(retrieverV_connect)$centralization # degree centrality 0.1929825
rc_betc <- centr_betw(retrieverV_connect)$centralization # betweenness centrality 0.3070175
rc_cloc <- centr_clo(retrieverV_connect)$centralization # closeness centrality 0.3888632
rc_eigc <- centr_eigen(retrieverV_connect)$centralization # eigen centraility 0.6816848



#### Operation Span
spanV <- make_empty_graph(n=5, directed = FALSE)
spanV_edges <- c(1,3,1,4,1,5,2,4,3,4)
spanV <- add.edges(spanV, spanV_edges)
plot(spanV)
s_ed <- edge_density(spanV) # edge density 0.5
s_degc <- centr_degree(spanV)$centralization # degree centrality 0.25
s_betc <- centr_betw(spanV)$centralization # betweenness centrality 0.375
s_cloc <- centr_clo(spanV)$centralization # closeness centrality 0.4277778
s_eigc <- centr_eigen(spanV)$centralization # eigen centraility 0.4209886

spanmat <- as_adjacency_matrix(spanV,sparse = FALSE)
span_spec <- eigen(spanmat)
spanlapmat <- laplacian_matrix(spanV, sparse=FALSE)
spanlap_spec <- eigen(spanlapmat)
