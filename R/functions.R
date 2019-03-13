# startPrint <- function(){
#   output$out <- renderText({ 
#     
#     
#     curr.txt  <- print('I am here<br>') 
#     startTxt  <<- paste(startTxt, curr.txt, collapse = '')
#     
#     # Call custom javascript to scroll window
#     session$sendCustomMessage(type = "scrollCallback", 1)
#     
#     return(startTxt)
#   })
# }

getClusterNodeLabels <- function(rvalues, explode.net){
  
  if(!is.null(rvalues$clusterNum[1]) && rvalues$clusterNum!=0){
    pointsInCluster <- explode.net[explode.net[,5]==rvalues$clusterNum[1], ]
    
    varLabels <- pointsInCluster[pointsInCluster[,6]==2,1]
    entityLabels <- pointsInCluster[pointsInCluster[,6]==1,1]
    labels.str <- ""
    labels.str <- paste(labels.str, "\n                    VARIABLES (", length(varLabels), ")", sep = "")
    
    
    for(i in 1:length(varLabels)){
      if(i>0 && i<10){
        varWithNum <- paste(i, ".       ", varLabels[i], sep = "")  
      }
      else if(i>=10 && i<100){
        varWithNum <- paste(i, ".     ", varLabels[i], sep = "")  
      }
      else if(i>=100){
        varWithNum <- paste(i, ".   ", varLabels[i], sep = "")  
      }
      labels.str <- paste(labels.str , varWithNum  , sep = "\n")
    }
    entityTitle <- paste("\n                     ENTITIES (",length(entityLabels), ")", sep = "")
    labels.str <- paste(labels.str, entityTitle , sep = "\n")
    
    for(i in 1:length(entityLabels)){
      if(i>0 && i<=9){
        entityWithNum <- paste(i,  ".       ", entityLabels[i], sep = "")
      }
      else if(i>=10 && i<=99){
        entityWithNum <- paste(i, ".     ", entityLabels[i], sep = "")
      }
      else if(i>=100 && i<=999){
        entityWithNum <- paste(i, ".   ", entityLabels[i], sep = "")
      }
      else if(i>=1000){
        entityWithNum <- paste(i, ".  ", entityLabels[i], sep = "")
      }
      #entityWithNum <- paste(i, ". ", entityLabels[i], sep = "")
      labels.str <- paste(labels.str ,  entityWithNum, sep = "\n")
    }
  }
  else{
    labels.str <- "Click on Entity/Variable to view the labels"
  }  
  
 return(labels.str) 
}

renderlegend <- function(input){
  outfile <- tempfile(fileext='.png')
  # Generate the PNG
  png(outfile, width = 150, height = 30, bg = "#272b30")
  par(mar = c(0,0,0,0), bg = "white")
  plot(c(0, 1), c(0,2 ), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  
  
  if(input$visualDirectorypath != "Folder not selected" && isStatusOK(input$visualDirectorypath)){
    pathToData <- gsub('Results','Data/', input$visualDirectorypath)    
    prefFilePath <- paste(pathToData, "/Preferences.txt", sep = "")
    pref <- as.matrix(read.table(prefFilePath, sep = "\t", header = FALSE)) #store user preferences in pref variable.  
    
    
    if(input$nodeColor == 3) # entities and variables
    {
      text(  x = 0.57,   y = 1.05,pref[1,2], cex = 0.7, col = "black") # ptn
      points(x = 0.1, y = 1.05, pch = 16, cex = 1.0, col = "red")
      text(  x = 0.28,   y = 1.65, pref[2,2], cex = 0.7, col = "black")  #snp
      points(x = 0.1, y = 1.55, pch = 17, cex = 0.8, col = "black")
      segments(0.07, 0.35, 0.18, 0.35, col= 'gray60') #edge
      text(  x = 0.39,   y = 0.4, "Association ", cex = 0.7, col = "black")  #edge
      
    }
    else if(input$nodeColor == 1) # outcomes
    {
      text(  x = 0.75,   y = 1.0, pref[3,2], cex = 0.7, col = "black")  #control
      points(x = 0.57, y = 0.95, pch = 16, cex = 1.3, col = "sky blue")
      text(  x = 0.35,   y = 1.0, pref[4,2], cex = 0.7, col = "black")  #case
      
      points(x = 0.2, y = 0.95, pch = 16, cex = 1.3, col = "magenta")
      
      text( x = 0.5,   y = 1.75, pref[2,2], cex = 0.7, col = "black")  #snp
      points(x = 0.2, y = 1.65, pch = 15, cex = 1.2, col = "black")
      
      segments(0.2, 0.35, 0.3, 0.35, col= 'gray60') #edge
      text(  x = 0.6,   y = 0.4, "Association ", cex = 0.7, col = "black")  #edge
      
    }
    else # cluster
    {
      text(x = 0.5, y = 1.5, paste("Node Colors = Cluster Membership"), cex = 0.7, col = "black")
      segments(0.07, 0.6, 0.28, 0.6, col= 'gray60') #edge
      text(  x = 0.6,   y = 0.7, "Association ", cex = 0.7, col = "black")  #edge
      
    }
  }
  else{
    text(x = 0.5, y = 1, paste("Legend"), cex = 0.75, col = "black")
  }
  dev.off()
  
  # Return a list containing the filename
  res <- list(src = outfile,
          contentType = 'image/png',
            width = 150, height = 30,
              alt = "Select Folder to view the legend")
  return(res)
}
getNetwork <- function(session, resultsFolder, displace, rvalues){
  
      path <- paste(resultsFolder, '/coordinates.dat', sep="")
      coordinates <- as.data.frame(fread(path , sep ="\t"))
      names <- colnames(coordinates)
      coordX <- paste('EL.',names[2],displace, sep = "" )
      coordY <- paste('EL.',names[3],displace, sep = "" )
      explodeNetwork <- cbind.data.frame(coordinates[,1], coordinates[,colnames(coordinates) == coordX],
                                         coordinates[,colnames(coordinates) == coordY], coordinates[,4],
                                         coordinates[,5], coordinates[,6])
      colnames(explodeNetwork) <- c(names[1], coordX, coordY, names[4], names[5], names[6] )
      return (explodeNetwork)
}

renderDynamicPlot <- function(input){
  return(shiny::renderImage({ 
    outfile <- tempfile(fileext='.png')
    if(input$visualDirectorypath != "Folder not selected" && isStatusOK(input$visualDirectorypath))
    { 
      # A temp file to save the output.
      # This file will be removed later by renderImage
      outfile <- tempfile(fileext='.png')
      
      # Generate the PNG
      png(outfile, width = 300, height = 280, bg = "#272b30")
      filePath <- paste(input$visualDirectorypath,"/statList.dat",sep = "")
      createELGraph(filePath, input$displace)
      
    }
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 300, height = 280,
         alt = "Folder not selected")
  }, deleteFile = TRUE))
}

createELGraph <- function(filePath, radiusSelected){
  dataframe <- read.delim(file = filePath, header = TRUE, sep = "\t") 
  
  radius <- dataframe[,1]
  ccs <- dataframe[,2]
  
  xmin <- min(dataframe[,1]) ##min(radius)
  xmax <- max(dataframe[,1])#2.0
  
  #minimum <- signif(min(ccs), digits = 2)
  #minimum <- signif(min(pmin(modim, coCluster, na.rm = TRUE)), digits = 2)
  #ymin <- minimum - 0.1*minimum
  ymin <- 0.0 
  #maximum <- signif(max(pmax(modim, coCluster, na.rm = TRUE)), digits = 2)
  #maximum <- signif(max(ccs), digits = 2)
  #ymax <- maximum + 0.1*maximum
  ymax <- 1.0
  
  #Creates a 3x1 grid where the first plot covers the first 2 rows 
  #pp <- layout(matrix(c(1,1,1), 1, 1, byrow=F))
  
  #Set plot margins
  par(mar = c(2,2,1,1)) 
  par( bg = "#272b30", fg="lightgrey")
  
  #plot ccs vs explode radisu graph. Change name of plot manually
  plot(radius, ccs, type="o", col="grey",
       xlab="", 
       ylab="", 
       lwd = 3, lty = 4,
       axes = "false",
       ylim = c(ymin,ymax),
       frame=T,
       col.lab="grey",
       col.main="grey"
  )
  
  cord.x <- c(-0.135, radiusSelected, radiusSelected, -0.135)
  cord.y <- c(-0.035, -0.035, 1.035, 1.035)
  polygon(cord.x,cord.y, density = 200, col = "#428bca" , fillOddEven = TRUE)
   
  #find optimal radius
  dataWithoutNA  <- na.omit(dataframe)
  m <- max(dataWithoutNA[,2])
  radiusOptimal <- as.numeric(dataWithoutNA[dataWithoutNA[,2]==m,1][1])
  abline(v=radiusOptimal, col='grey', lty=2, lwd=2)
  
  #adjust axes parameters, create major and minor ticks, add labels
  #axis(1, at=seq(xmin,xmax,by=0.1), lwd=1, lwd.ticks=1, col.ticks = "grey",col.axis = "grey")
  axis(2, at=seq(0,1.0, by = 0.2),lwd=1, lwd.ticks=1, col.ticks = "grey",col.axis = "grey")
  
  #redraw points
  points(radius, ccs,col="grey",lwd=3)
  lines(radius, ccs,col="grey",lwd=2)
  
  legend(x = "topright", c("CCS Score"), lty =c(1), col=c("grey"), border = "grey", lwd = c(2))
  dev.off()
}


getVertexSize <- function(vertexSize, session, input){
  
  if ((vertexSize < 1L ) || is.na(vertexSize)) {
    #updateNumericInput(session, "vertexSize", value = 1)
    return(0)
  } else {
    if (vertexSize > 10L ) {
      #updateNumericInput(session, "vertexSize", value = 10)
      return(10L)
    } else {
      return (vertexSize)
    }
  }
}

getThickness <- function(input){
  edgeThickness <- input$edgeThickness
  
  if((edgeThickness < 0.01) || is.na(edgeThickness)) {
    #updateNumericInput(session, "labelSize", value = 0)
    return(0)
  }
  else 
  {
    if(edgeThickness > 1 ) 
    {
      #updateNumericInput(session, "labelSize", value = 8)
      return(1L)
    }
    else 
    {
      return(as.numeric(edgeThickness))
    }
  }
}


isStatusOK <- function(path){
  result = tryCatch({ 
    fileName <- paste(path, "/status.dat", sep = "/")
    return(read.table(fileName) == "OK")
  }, warning = function(w) {
    return(FALSE)
  }, error = function(e) {
    return(FALSE)
  })
}

updateReactiveValue <- function(session, rvalues, msg){
  #updateTextInput(session,'progressText', value="Status : IDLE")
  
  rvalues$i <- 0
  rvalues$CS_list <- 0
  rvalues$best <- 0
  rvalues$attrs <- 0
  #suppressWarnings(library(sna,  warn.conflicts = FALSE, quietly = TRUE))
  
  #debug:initialize status.dat file
  #status <- paste(rvalues$pathToSaveResults,"status.dat", sep = "")
  #cat(paste("hiii"), file=status, append = TRUE)
  
  fileConn <- c(wd=paste(rvalues$resultsFolder,'/status.dat', sep = ''))
  #open(fileConn, open ='rw')
  #if(fileConn != "" ){
  writeLines(msg, fileConn)
  #close(fileConn)
  #}
  
}


getExplodeNetwork <- function(resultsFolder){
  
  path <- paste(resultsFolder, '/coordinates.dat', sep="")
  if(file.exists(path)){
      explodeNetwork <- data.table::fread(path , sep ="\t")
      explodeNetwork <- data.frame(explodeNetwork)
      return (explodeNetwork)
  }
}

formatExplodeNetMatrix <- function(explode.net, rvalues, input){
  #this method changes colnames; adds entity and outcome to the explode.net matrix
  
  if(input$layoutAlgoType == 'fr'){
    colnames(explode.net) <- c("Label", "Cluster", "FRX", "FRY")
  }
  else{
    colnames(explode.net) <- c("Label", "Cluster", "KKX", "KKY")
  }
  
  colnames(rvalues$entity) <- c("Label", "Entity")
  #Modified: AA
  #add entity column to the explode.net
  
  explode.net <- as.matrix(merge(explode.net, rvalues$entity, by = "Label", sort= FALSE))
  if(input$layoutAlgoType == 'fr'){
    colnames(explode.net) <- c("Label", "Cluster", "FRX", "FRY", "Entity")
  }
  else{
    colnames(explode.net) <- c("Label", "Cluster", "KKX", "KKY", "Entity")
  }
  
  #add outcome file if available
  if(!is.null(rvalues$outcome)){
    explode.net <- cbind(explode.net, rvalues$outcome[,2])
    colnames(explode.net)[6] <- "Outcome"
  }
  
  return(explode.net)
}

updateTextArea <- function(session, msg, input){
  val <- paste(input$searchUpdate, msg, sep = "\n")
  updateTextInput(session,'searchUpdate', value=val)
  session$sendCustomMessage(type = "scrollCallback", 1)
}


getArrayofEdgelist <- function(resultsFolder) {
  
  # Read edgelist file and create graphList 
  edgeList <- read.delim(paste(resultsFolder,"/", i, sep=""), sep="\t", 	header=TRUE)
  edgeList <- data.matrix(edgeList)
  edgeList <- edgeList[(edgeList[,3]>0),]	
  
  graphList[[n]] <- graph.data.frame(edgeList, directed=FALSE)
  
  
  ##################
  # #create node list files and read it
  # listEdgeList <- list()
  # listNodeInfo <- list()
  # n <- 2
  # j <- 1
  # detach("package:igraph")
  # for(i in edgeListFiles)
  # {
  #   nodeInfo <- read.delim(paste(resultsFolder,"/", paste(nodeListFiles[j]), 	sep=""), sep="\t")
  #   #plotcord <- plotVal[,2:3]
  #   
  #   graph <-  graphList[[n]]
  #   net <- asNetwork(graph)
  #   edglist <- as.matrix.network.edgelist(net)
  #   listEdgeList[[n]] <- edglist
  #   listNodeInfo[[n]] <- nodeInfo
  #   n <- n + 1
  #   j <- j + 1
  # }
  # library(igraph,  warn.conflicts = FALSE, quietly = TRUE)
  # 
  #returnList <- list("listEdgeList" = listEdgeList, "listNodeInfo" = listNodeInfo)
  
  #return(returnList)
  return()
}