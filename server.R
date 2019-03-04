# helper function to apply D.BIC and correlation cut-offs
ApplyCutOffs <- function(path.dframe,D.BIC.cut,cor.cut){
  ind <-  intersect(which(path.dframe$D.BIC > D.BIC.cut),
                    which(abs(path.dframe$PathCor) > cor.cut))
  return(path.dframe[ind,])
}


ExtractTopConnected <- function(path.dframe,path.target,top.n=25,topMethod="abs",BIC.cut=0.75,cor.cut=0.25){
  # Function to return the top n pathways connected with the target pathway
  # 
  # Args:
  #   path.dframe: a data frame with the names of the pathways
  #                (Pathway.A and Pathway.B) and their corresponding
  #                value for the pathway correlation (PathCor)
  #   path.target: string with the name of the pathway of interest
  #   top.n: the cut-off for the top n pathways connected to the 
  #          the pathway of interest.
  #   topMethod: method to sort the edge weights (PathCor). There are
  #              3 methods available
  #              1. "abs" sort by the absolute value of the edge weight
  #                 (|PathCor|)
  #              2. "decreasing" sort the edge weight in decreasing order
  #              2. "increasing" sort the edge weight in increasing order
  #
  # Returns:
  #   A data frame with the same format as the input data frame path.dframe
  #   with the top n pathways connected to the pathway of interest
    
  
  path.dframe <- ApplyCutOffs(path.dframe,BIC.cut,cor.cut)
  
  # get all edges for target pathway
  edges.vec <- paste(path.dframe$Pathway.A,path.dframe$Pathway.B,sep="_")
  path.ind <- grep(path.target,edges.vec,fixed=T)
  path.subnet <- path.dframe[path.ind,]
  # get top edges based on type of ordering
  selectMethod <- c("abs","decreasing","increasing") %in% topMethod
  if(selectMethod[1]){
    # sort by |PathCor|
    path.subnet <- path.subnet[order(abs(path.subnet$PathCor),decreasing=T),]
  }
  if(selectMethod[2]){
    path.subnet <- path.subnet[order(path.subnet$PathCor,decreasing=T),]
  }
  if(selectMethod[3]){
    path.subnet <- path.subnet[order(path.subnet$PathCor,decreasing=F),]
  }
  if(sum(selectMethod)==0){
    print("Not a valid choice to sort edges")
    return()
  }
  # return top edges
  if(length(path.ind) >= top.n){
    return(path.subnet[1:top.n,])
  }else{
    return(path.subnet)
  }
  
}

SwapPathwayNames <- function(path.dframe,path.target,top.n){
  # Auxiliary function to swap the pathways names in a data frame so that
  # the target pathway name is in the Pathway.A and the top n connected
  # pathway names are in the Pathway.B column
  #
  # Args:
  #   path.dframe: a data frame with the names of the pathways
  #                (Pathway.A and Pathway.B) 
  #   path.target: string with the name of the pathway of interest
  #   top.n: the cut-off for the top n pathways connected to the 
  #          the pathway of interest.
  #
  # Returns:
  #   A data frame with the same format as the input data frame path.dframe
  
  swap.ind <- which(path.dframe$Pathway.B[1:top.n] == path.target)
  path.dframe$Pathway.B[swap.ind] <- path.dframe$Pathway.A[swap.ind]
  path.dframe$Pathway.A[swap.ind] <- path.target
  return(path.dframe)
}

ExtractSubNet <- function(path.dframe,path.target,top.n=25,top.k=NULL,topMethod="abs",BIC.cut=10){
  # Function to return the top n pathways connected with the target pathway and the
  # edges among the top n pathways connected with the pathway of interest
  #
  # Args:
  #   path.dframe: a data frame with the names of the pathways
  #                (Pathway.A and Pathway.B) and their corresponding
  #                value for the pathway correlation (PathCor)
  #   path.target: string with the name of the pathway of interest
  #   top.n: the cut-off for the top n pathways connected to the 
  #          the pathway of interest.
  #   top.k: the cut-off for the number of edges to include among the
  #          top connected pathways. If top.n=NULL, then all connections
  #          between the top pathways are considered
  #   topMethod: method to sort the edge weights (PathCor). There are
  #              3 methods available
  #              1. "abs" sort by the absolute value of the edge weight
  #                 (|PathCor|)
  #              2. "decreasing" sort the edge weight in decreasing order
  #              2. "increasing" sort the edge weight in increasing order
  #
  # Returns:
  #   A data frame with the same format as the input data frame path.dframe
  #   with the top n pathways connected to the pathway of interest and the 
  #   connections among the top connected pathways
  
  
  subnet.dframe <- ExtractTopConnected(path.dframe,path.target,top.n,topMethod,BIC.cut)
  # swap pathway names
  subnet.dframe <- SwapPathwayNames(subnet.dframe,path.target,top.n)
  path.list <- subnet.dframe$Pathway.B[1:top.n]
  tmp.dframe <- as.data.frame(c())
  for(i in 1:length(path.list)){
    # find all edges 
    edges.vec <- paste(path.dframe$Pathway.A,path.dframe$Pathway.B,sep="_")
    path.ind <- grep(path.list[i],edges.vec,fixed=T)
    #     tmp.ind <- path.dframe$Pathway.A %in% path.list[i]
    #     tmp.ind <- cbind(tmp.ind,path.dframe$Pathway.B %in% path.list[i])
    #     path.ind <- apply(tmp.ind,1,function(x){x[1]||x[2]})
    path.edges <- path.dframe[path.ind,]
    # swap pathway names
    path.edges <- SwapPathwayNames(path.edges,path.list[i],dim(path.edges)[1])
    path.match <- path.edges$Pathway.B %in% path.list[-1*1:i]
    # add edges
    tmp.dframe <- rbind(tmp.dframe,path.edges[path.match,])
  }
  # sort weighted edges
  if(is.numeric(top.k)){
    selectMethod <- c("abs","decreasing","increasing") %in% topMethod
    
    if(selectMethod[1]){
      # sort by |PathCor|
      tmp.dframe <- tmp.dframe[order(abs(tmp.dframe$PathCor),decreasing=T),]
    }
    if(selectMethod[2]){
      tmp.dframe <- tmp.dframe[order(tmp.dframe$PathCor,decreasing=T),]
    }
    if(selectMethod[3]){
      tmp.dframe <- tmp.dframe[order(tmp.dframe$PathCor,decreasing=F),]
    }
    if(sum(selectMethod) == 0){
      print("Not a valid choice to sort edges")
      return()
    }
    tmp.dframe <- tmp.dframe[1:top.k,]
  }
  subnet.dframe <- rbind(subnet.dframe,tmp.dframe)
  return(subnet.dframe)
}
# 
# ExtractTopConnected2 <- function(path.dframe,path.target,top.n=25,topMethod="abs",
#                                  BIC.cut=10,edge.weight="PathCor"){
#   # Function to return the top n pathways connected with the target pathway
#   # 
#   # Args:
#   #   path.dframe: a data frame with the names of the pathways
#   #                (Pathway.A and Pathway.B) and their corresponding
#   #                value for the pathway correlation (PathCor)
#   #   path.target: string with the name of the pathway of interest
#   #   top.n: the cut-off for the top n pathways connected to the 
#   #          the pathway of interest.
#   #   topMethod: method to sort the edge weights (PathCor). There are
#   #              3 methods available
#   #              1. "abs" sort by the absolute value of the edge weight
#   #                 (|PathCor|)
#   #              2. "decreasing" sort the edge weight in decreasing order
#   #              2. "increasing" sort the edge weight in increasing order
#   #
#   # Returns:
#   #   A data frame with the same format as the input data frame path.dframe
#   #   with the top n pathways connected to the pathway of interest
#   
#   path.dframe <- path.dframe[which(path.dframe$D.BIC > BIC.cut),]
#   # get all edges for target pathway
#   edges.vec <- paste(path.dframe$Pathway.A,path.dframe$Pathway.B,sep="_")
#   path.ind <- grep(path.target,edges.vec,fixed=T)
#   path.subnet <- path.dframe[path.ind,]
#   # get top edges based on type of ordering
#   selectMethod <- c("abs","decreasing","increasing") %in% topMethod
#   if(selectMethod[1]){
#     # sort by |PathCor|
#     tmp.exp <- paste("path.subnet <- path.subnet[order(abs(path.subnet$",
#                      edge.weight,"),decreasing=T),]",sep="")
#     
#   }
#   if(selectMethod[2]){
#     tmp.exp <- paste("path.subnet <- path.subnet[order(path.subnet$P",
#                      edge.weight,",decreasing=T),]",sep="")
#   }
#   if(selectMethod[3]){
#     tmp.exp <- paste("path.subnet <- path.subnet[order(path.subnet$",
#                      edge.weight,",decreasing=F),]",sep="")
#   }
#   if(sum(selectMethod)==0){
#     print("Not a valid choice to sort edges")
#     return()
#   }
#   eval(parse(text=tmp.exp))
#   # return top edges
#   return(path.subnet[1:top.n,])
# }



ExtractSubNet2 <- function(path.dframe,path.target,top.n=25,top.k=NULL,topMethod="abs",
                           BIC.cut=10,cor.cut=0.25,edge.weight="PathCor"){
  # Function to return the top n pathways connected with the target pathway and the
  # edges among the top n pathways connected with the pathway of interest
  #
  # Args:
  #   path.dframe: a data frame with the names of the pathways
  #                (Pathway.A and Pathway.B) and their corresponding
  #                value for the pathway correlation (PathCor)
  #   path.target: string with the name of the pathway of interest
  #   top.n: the cut-off for the top n pathways connected to the 
  #          the pathway of interest.
  #   top.k: the cut-off for the number of edges to include among the
  #          top connected pathways. If top.n=NULL, then all connections
  #          between the top pathways are considered
  #   topMethod: method to sort the edge weights (PathCor). There are
  #              3 methods available
  #              1. "abs" sort by the absolute value of the edge weight
  #                 (|PathCor|)
  #              2. "decreasing" sort the edge weight in decreasing order
  #              2. "increasing" sort the edge weight in increasing order
  #
  # Returns:
  #   A data frame with the same format as the input data frame path.dframe
  #   with the top n pathways connected to the pathway of interest and the 
  #   connections among the top connected pathways
  
  
  subnet.dframe <- ExtractTopConnected(path.dframe,path.target,top.n,topMethod,BIC.cut,cor.cut)
  if(dim(subnet.dframe)[1] < top.n){
    top.n <- dim(subnet.dframe)[1]
  }
  # swap pathway names
  subnet.dframe <- SwapPathwayNames(subnet.dframe,path.target,top.n)
  path.list <- subnet.dframe$Pathway.B[1:top.n]
  tmp.dframe <- as.data.frame(c())
  for(i in 1:length(path.list)){
    # find all edges 
    edges.vec <- paste(path.dframe$Pathway.A,path.dframe$Pathway.B,sep="_")
    path.ind <- grep(path.list[i],edges.vec,fixed=T)
    #     tmp.ind <- path.dframe$Pathway.A %in% path.list[i]
    #     tmp.ind <- cbind(tmp.ind,path.dframe$Pathway.B %in% path.list[i])
    #     path.ind <- apply(tmp.ind,1,function(x){x[1]||x[2]})
    path.edges <- path.dframe[path.ind,]
    # swap pathway names
    path.edges <- SwapPathwayNames(path.edges,path.list[i],dim(path.edges)[1])
    path.match <- path.edges$Pathway.B %in% path.list[-1*1:i]
    # add edges
    tmp.dframe <- rbind(tmp.dframe,path.edges[path.match,])
  }
  # sort weighted edges
  if(is.numeric(top.k)){
    selectMethod <- c("abs","decreasing","increasing") %in% topMethod
    
    if(selectMethod[1]){
      # sort by |PathCor|
      tmp.exp <- paste("tmp.dframe <- tmp.dframe[order(abs(tmp.dframe$",
                       edge.weight,"),decreasing=T),]",sep="")
    }
    if(selectMethod[2]){
      tmp.exp <- paste("tmp.dframe <- tmp.dframe[order(tmp.dframe$",
                       edge.weight,",decreasing=T),]",sep="")
    }
    if(selectMethod[3]){
      tmp.exp <- paste("tmp.dframe <- tmp.dframe[order(tmp.dframe$",
                       edge.weight,",decreasing=F),]",sep="")
    }
    if(sum(selectMethod) == 0){
      print("Not a valid choice to sort edges")
      return()
    }
    eval(parse(text=tmp.exp))
    tmp.dframe <- tmp.dframe[1:top.k,]
  }
  subnet.dframe <- rbind(subnet.dframe,tmp.dframe)
  return(subnet.dframe)
}


TopTable <- function(path.dframe,top.n,edge.type,path.target){
  path.dframe <- path.dframe[1:top.n,]
  path.dframe <- SwapPathwayNames(path.dframe,path.target,top.n)
  tab <- path.dframe$Pathway.B
  if(edge.type == "PathCor"){
    tab <- cbind(tab,path.dframe$PathCor)
    colnames(tab) <- c("Pathway","PathCor")
    return(tab)
  }
  if(edge.type == "Overlap.Coeff"){
    tab <- cbind(tab,path.dframe$Overlap.Coeff)
    colnames(tab) <- c("Pathway","Overlap.Coeff")
    return(tab)
  }
  if(edge.type == "both"){
    tab <- cbind(tab,path.dframe$PathCor)
    tab <- cbind(tab,path.dframe$Overlap.Coeff)
    colnames(tab) <- c("Pathway","PathCor","Overlap.Coeff")
    return(tab)
  }
}

TopTable2 <- function(path.dframe,top.n,path.target){
  path.ind <- c(grep(path.target,path.dframe$Pathway.A,fixed=T),
                grep(path.target,path.dframe$Pathway.B,fixed=T))
  if(length(path.ind) < top.n){
    top.n <- length(path.ind)
  }
  path.dframe <- path.dframe[1:top.n,]

  
  path.dframe <- SwapPathwayNames(path.dframe,path.target,top.n)
  tab <- data.frame(Pathway=rep(NA,top.n),
                    PathCor=rep(NA,top.n),
                    Cor=rep(NA,top.n),
                    Overlap.Coeff=rep(NA,top.n))
  tab$Pathway <- path.dframe$Pathway.B
  tab$PathCor <- path.dframe$PathCor
  tab$Cor     <- path.dframe$Cor
  tab$Overlap.Coeff <-path.dframe$Overlap.Coeff
  return(tab)
}


GetPathCol <- function(subnet.dframe,db.col){
  # Creates a tab delimited text file with the pathway names and
  # their database. The text file can be imported into Cytoscape
  # as node attributes
  #
  # Args:
  #   subnet.dframe: a data frame with three columns, the first two
  #                  contain the pathway names and the third the edge weight
  #   file.name: file name for the text file
  # Returns:
  #   a text file [file.name].txt
  path.names <- unique(c(subnet.dframe$Pathway.A,subnet.dframe$Pathway.B))
  #tag <- c("kegg","wikipathways","reactome","netpath","static modules")
  path.db <- c("KEGG","Wikipathways","Reactome","Netpath","Static Module")
  
  if(length(path.db) != length(db.col)){
    warning("Color palette doesn't match tags.")
    return(0)
  }
  
  # path.tag <- rep(NA,length(path.names))
  path.col <- rep(NA,length(path.names))
  for(i in 1:length(path.db)){
    mtch <- grep(path.db[i],path.names,fix=T)
    if(is.vector(mtch)){
      # path.tag[mtch] <- tag[i]
      path.col[mtch] <- db.col[i]
    }
  }
  res <- data.frame(cbind(path.names,path.col))
  res[,2] <- as.character(res[,2])
  #res <- cbind(res,path.col)
  colnames(res) <- c("ID","Node.Color")
  return(res)
}

GetEdgeCol <- function(edge.weight){
  edge.col <- rep(NA,length(edge.weight))
  edge.col[which(edge.weight >= 0)] <- rgb(0,0.4,0,0.65)
  edge.col[which(edge.weight < 0)] <- rgb(0.8,0.1,0.1,0.75)
  return(edge.col)
}

PlotSubnet <- function(subnet.dframe,node.col){
#   # Node color by pathway databases
#   node.col <- brewer.pal(n=8,name="Set2")
#   node.col <- node.col[-1*5:7]
  
  subnet.col <- GetPathCol(subnet.dframe=subnet.dframe,db.col=node.col)
  # color edges by positive/negative correlations
  edge.col <- GetEdgeCol(subnet.dframe$PathCor)
  
  g <- graph.data.frame(subnet.dframe,directed=F)

  plot(g,vertex.color=subnet.col$Node.Color,vertex.size=17,
       vertex.label.font=11,vertex.label.color="black",
       vertex.label.cex=1,vertex.frame.color="black",
       edge.width=abs(E(g)$PathCor)*7,edge.color=edge.col)
  
}

FormatPathwayList <- function(path.dframe,node.col){
  tmp.names <- unique(c(path.dframe$Pathway.A,path.dframe$Pathway.B))
  tmp.links <- matrix(NA,nrow=dim(path.dframe)[1],ncol=3)
  colnames(tmp.links) <- c("source","target","weight")
  for(i in 1:length(tmp.names)){
    tmp.ind00 <- path.dframe$Pathway.A %in% tmp.names[i]
    tmp.links[tmp.ind00,1] <- i-1
    tmp.ind01 <- path.dframe$Pathway.B %in% tmp.names[i]
    tmp.links[tmp.ind01,2] <- i-1
  }
  tmp.links[,3] <- abs(path.dframe[,3])
  tmp.col <- GetPathCol(subnet.dframe=path.dframe,db.col=node.col)
  return(list(names=tmp.names,nodeCol=tmp.col$Node.Color,links=tmp.links))
}

BICplot <- function(path.dframe,D.BIC.cutoff,cor.cutoff){
  # function to approximate BIC curve linearly for a smoother plot
  BIC.fun <- approxfun(x=c(path.dframe$PathCor),y=c(path.dframe$D.BIC))
  x <- seq(-1,1,length.out=1000)
  plot(x,BIC.fun(x),type="l",xlim=c(-.8,0.8),ylim=c(0,59730),lwd=2,
       ylab=expression(paste(Delta,"BIC")),xlab="PathCor",cex.lab=1.5)
  abline(h=D.BIC.cutoff )
  abline(v=cor.cutoff)
  abline(v=-cor.cutoff)
  polygon(x=c(-0.86,-cor.cutoff,-cor.cutoff,-0.86),
          y=c(D.BIC.cutoff,D.BIC.cutoff,BIC.fun(.96),BIC.fun(.96)),
          border=NA,density=10,angle=45)
  polygon(x=c(0.86,cor.cutoff,cor.cutoff,0.86),
          y=c(D.BIC.cutoff,D.BIC.cutoff,BIC.fun(.96),BIC.fun(.96)),
          border=NA,density=10,angle=45)
  
}

library(shiny)
library(igraph)
# library(RColorBrewer)

shinyServer(function(input, output){
  output$caption <- renderText({
    input$tool1
  })
  
  datasetInput <- reactive({
    sMethod <- switch(input$method,
                       "Absolute Value" = "abs",
                       "Decreasing" = "decreasing",
                       "Increasing" = "increasing")
    na.omit(ExtractSubNet2(path.dframe=path.net.dframe,
                  path.target=input$tool1,
                  top.n=input$top.n,
                  top.k=5,
                  topMethod=sMethod,
                  BIC.cut=quantile(path.net.dframe$D.BIC,input$BIC.cut,na.rm=T),#BIC.cut=input$BIC.cut, 
                  cor.cut=input$cor.cut,
                  edge.weight=input$edge.weight))
  })
  # for d3.js network
  output$net.list <- reactive({
    FormatPathwayList(datasetInput(),node.col)
  })
  
  output$view <- renderTable({
    TopTable2(datasetInput(),input$top.n,input$tool1)
  },digits=4)
  
  # plot graph with pathways of interest
  output$graphplot <- renderPlot({
    par(mar = c(0.25,0.25,0.25,0.25), oma= c(0.05,0.05,0.05,0.05))
    PlotSubnet(datasetInput(),node.col)
    
  })
  
  output$BICplot <- renderPlot({
    # par(mar=c(bottom,left,top,right))
    par(pin=c(3.75,3.5))
    BICplot(path.dframe=path.net.dframe,
            D.BIC.cutoff=quantile(path.net.dframe$D.BIC,input$BIC.cut,na.rm=T),
            cor.cutoff=input$cor.cut)
    })

})