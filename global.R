# setwd("C:/Users/Yered/Copy/Shiny/PathCoexp06")
# setwd("/home/yeredh/ShinyApps/PathCoexp06")

load("data/GPL570.gs.RData")

# mouse.gs <- readRDS("data/LSC.mouse.gs.RDS")

path.net.dframe <- readRDS("data/net.dframe.RDS")

path.net.dframe$PathCor <- round(path.net.dframe$PathCor,4)
path.net.dframe$Overlap.Coeff <- round(path.net.dframe$Overlap.Coeff,4)

# # Get the names for all pathways
# path.names <- sort(names(mouse.gs))
# path.list <- as.list(path.names)
# names(path.list) <- path.names

# Get the names for all pathways
path.names <- sort(names(GPL570.gs))
path.list <- as.list(path.names)
names(path.list) <- path.names


coeffSelString <- "The shaded area corresponds to the pathway correlation coefficients included in the network."

# Node color by pathway databases
library(RColorBrewer)
node.col <- brewer.pal(n=8,name="Set2")
node.col <- node.col[-1*5:7]



