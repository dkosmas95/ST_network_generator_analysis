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
library(readxl)
library(lubridate)

# read excel tabs
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)


#### read in graph
graphFile <- "/HT5_1_40.graphml" # insert file location here
gSoc <- read.graph(file=graphFile, format="graphml")
traf <- which(V(gSoc)$role=="trafficker")
bot <- which(V(gSoc)$role=="bottom")

# change to grayscale
V(gSoc)$color[which(V(gSoc)$color=="lightgoldenrod1")] <- "gray30"
V(gSoc)$color[which(V(gSoc)$color=="lightgreen")] <- "gray50"
V(gSoc)$color[which(V(gSoc)$color=="royalblue")] <- "gray70"
V(gSoc)$color[which(V(gSoc)$color=="lightskyblue")] <- "gray90"
E(gSoc)$color <- "darkgray"
E(gSoc)$lty[which(E(gSoc)$lty==2)] <- 3

V(gSoc)$shape <- "circle"
V(gSoc)$shape[traf] <- "square"
V(gSoc)$shape[bot] <- "triangle"

# layout traffickers in ring, with their operations radially off of them
gSoc <- delete.edges(gSoc, which(E(gSoc)$color=="lavenderblush3"))
tNodes <- traf
temp <- add.vertices(gSoc, 1)
for (i in 1:(length(tNodes)-1)){
  for (j in (i+1):length(tNodes)){
    if (get.edge.ids(temp,c(tNodes[i],tNodes[j]))==0){
      temp <- add.edges(temp, c(tNodes[i],tNodes[j]))
    }
  }
}
for (i in 1:length(tNodes)){
  temp <- add.edges(temp, c(tNodes[i],vcount(temp)))
}
layT <- layout_with_stress(temp)
lay <- layT[-vcount(temp),]

## adjust layout
lay[2,1] <- 2.25
lay[6,] <- c(1.6,-1.6)
lay[10,2] <- 2.1
lay[11,] <- c(0.01,1.85)
lay[12,1] <- 1.3
lay[13,] <- c(1.1,1.9)
lay[14,] <- c(1,0.9)
lay[15,1] <- -0.6
lay[16,] <- c(-0.07,1.25)
lay[23,] <- c(-1.45,-0.6)
lay[27,] <- c(-1.7,0)
lay[28,2] <- 0.55
V(gSoc)$label <- ""

plot(gSoc, layout=lay, vertex.size=15, edge.width=3)
roleList <- c("trafficker","bottom","victim, minor","victim, adult")
#colList <- c("lightgoldenrod1","lightgreen","lightskyblue","royalblue")
colList <- c("gray30","gray50","gray70","gray90")
legend("topleft",bty = "n",
       legend=roleList,
       fill=colList, border=NA, cex = 1.6)

#### isolate operations
vOpList <- list()
for (i in 1:(length(traf)-1)){
  vOpList[[i]] <- traf[i]:(traf[i+1]-1)
}
vOpList[[length(traf)]] <- traf[length(traf)]:vcount(gSoc)

gOpList <- list()
for (i in 1:length(traf)){
  gOpList[[i]] <- induced.subgraph(gSoc,vOpList[[i]])
}

layOp <- list()


for (i in 1:length(traf)){
  layOp[[i]] <- layout_nicely(gOpList[[i]])
  plot(gOpList[[i]], layout=layOp[[i]], vertex.size=18, edge.width=4)
  if (i==1){
    legend("topleft",bty = "n",
         legend=roleList,
         fill=colList, border=NA, cex = 2)
  }
}

## convert spreadsheet of results into plot
# column 1: attacker budget, column 2: MFNIP flow, column 3: flow restructuring in response to MFNIP interdictions, column 4: MFNIP-R flow
loc <- "/Computational Results.xlsx" # insert file location here
data <- read_excel_allsheets(loc)
tabNames <- excel_sheets(loc)

results <- data[[4]]
budgetA <- results[,1]
initFlow <- results[,2]
restFlow <- results[,3]
optFlow <- results[,4]
maxFlow <- 23 # correct for appropriate data set

plot(budgetA, initFlow, col="blue", pch=1, xlab="Attacker Budget", ylab="Flow", ylim = c(0,45), cex=1.8, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
title(main="Budget vs. Flow on Test Network",cex=1.8)
points(budgetA, restFlow, col="red", pch=8,cex=1.8)
points(budgetA, optFlow, col="green3", pch=9,cex=1.8)
lines(budgetA, rep(maxFlow,length(budgetA)), pch=19, col="black", type="b",cex=1.8)
legend(23,47,legend=c("Max Flow", "MFNIP Flow Before Restructuring","MFNIP Flow After Restructuring","MFNIP-R Flow"), col=c("black", "blue","red","green3"),
       pch = c(19,1,8,9),ncol=1, box.lty = 0, bty="n",cex=1.5)

