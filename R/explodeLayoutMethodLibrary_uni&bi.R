toggleValidity <- function(updateSliderValues)
{
  result <- tryCatch(
     if(updateSliderValues) { return(FALSE) }
     else { return(TRUE) },
     warning = function(w){return(FALSE)},
     error   = function(e){return(FALSE)})
  return(result)
}

#Define Min-Max normalization function
minmax <- function(vec){
  vec <- as.numeric(vec)
  mi <- min(vec)
  ma <- max(vec)

  output <- sapply(vec, function(x) ((x-mi)/(ma-mi)))
  return(as.numeric(output))
  #return(as.numeric(vec))
}

#Find the co-ordinates from file/FR/KK layout algorithm
getOriginalNetworkCoordinates <-function(input, g, userCoords){

  if(!input$coordfileCheckBox)
  {
    algoType <- input$layoutAlgoType
    vac <- data.frame()
    coord <- getCoordinatesForTheNetwork(g, vac, algoType)
  }
  else
  {
    #coordFile <- paste(basePath, '/', input$coordFile, sep = "")
    coord <- getCoordinatesForTheNetwork(g, userCoords, NULL)
  }

  return(coord)
}

#Returns a cleaned up data frame read from the file.
#If a row is all zero or has NA in it then its removed
#If column sum is greater than number of rows then it removed as well
getDataFrameFromFile <- function(networkfile_name, rvalues){

  df  <- read.delim(networkfile_name, sep="\t", row.names=1, check.names = FALSE )
  #df  <- read.table(networkfile_name, header = F, sep="\t",check.names = FALSE )

  #set all na values to zero
  df[is.na(df)] <- 0

  #in the data frame first column is phenotype, we are not using this information in the code
  phenotype <- df[,1, drop=FALSE]

  #clean up data frame to look like the one we need
  # All rows and columns with all zeros are removed and data frame is converted to matrix
  df <- df[,-1]
  df <- df[sapply(df, function(x) !any(is.na(x)))]

  #df <- df[rowSums(df[,]!=0, na.rm=TRUE)!=0,]
  #df <- df[,colSums(df==0, na.rm=TRUE)<nrow(df)]
  df <- data.matrix(df)

  return(df)
}

getGraphForLayout <- function(df, rvalues){
  #Convert data from matrix to edgelist
  edgelist <- data.matrix(df)
  if (rvalues$netType == 1) {
    n <- dim(df)[1]
    df2 <- df
    for (i in 1:(n-1)) {
      for (j in (1+i):n) {
        df2[i,j] <- 0
      }
    }
    edgelist <- data.matrix(df2)
    rownames(edgelist) <- 1:dim(edgelist)[1]
    colnames(edgelist) <- 1:dim(edgelist)[1]
    g <- igraph::graph_from_adjacency_matrix(data.matrix(df), mode = 'undirected', weighted = T)
  } else {
    elist <- data.table::melt(df)
    elist <- elist[(elist[,3]>0),]
    rownames(edgelist) <- 1:dim(edgelist)[1]
    colnames(edgelist) <- (dim(edgelist)[1]+1):(dim(edgelist)[1]+dim(edgelist)[2])
    g <- igraph::graph.data.frame(elist, directed=FALSE)
  }

  #generate edgelist file
  edgelist <- data.table::melt(edgelist)
  edgelist2 <- edgelist[edgelist[,3]!=0,]


  #filename2 <- paste(rvalues$resultsFolder,"/edgeList.dat", sep="")
  #filename2 <- "edgeList.dat"
  #colnames(edgelist2) <- c("Patient", "Variable", "Weight")
  #write.table(edgelist2, filename2  , sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  return(list(g, edgelist2))
}


#Extract co-ordinates from the the coordinates file, if its supplied by the user
#else generate coordinates using the Fruchterman Reingold Algorithm
getCoordinatesForTheNetwork <- function(g, coordinates = data.frame, algoType = NULL){
  #read co-ordinates from the text file, file should be supplied in the same folder as the other files

  if(length(coordinates) == 0)
  {
    if(algoType == 'fr'){
      coord <- suppressWarnings(igraph::layout_with_fr(g, niter = 1000, dim=2))
      #coord <- gplot.layout.fruchtermanreingold(as.matrix(g, matrix.type = "edgelist"))
    }
    else{
      coord <- suppressWarnings(igraph::layout_with_kk(g, niter = 1000, dim = 2))
    }
  }
  else
  {
    #use the user specified co-ordinates file
    #coord <- read.delim(coordinatesFile,  sep="\t" , header = FALSE)
    coord <- coordinates
    coord <- coord[match(igraph::V(g)$name, coord[,1]),2:3]

  }

  coord <- data.matrix(coord)

  return(coord)
}

#Calculates and returns angles on a circle at equal distance. Starting point is specified when the function is called
calculateEquidistantAngles <- function(nclust, rings){

  #Parameter to place centroid of first cluster
  degreeoffset <- 360/(nclust)

  #Find number of clusters on each circle starting from the outside circle
  clustersOnEachRing <- rep(0, rings)
  index <- as.numeric(rings)
  for (i in 1: nclust)
  {
    clustersOnEachRing[index] <- clustersOnEachRing[index] + 1
    index <- index - 1
    if(index < 1){ index <- as.numeric(rings) }
  }

  #Each circle needs to have a starting degree at which first cluster of the circle will go to
  startDegreeForEachRing <- c()
  for (i in 1: rings){
    startDegreeForEachRing <- c(startDegreeForEachRing, (i-1)*degreeoffset)
  }
  startDegreeForEachRing <- rev(startDegreeForEachRing)

  #Find the offset for each circle, offet is the angle between two clusters from the center of the circle
  offsetForEachRing <- c()
  for ( i in 1:rings){
    offsetForEachRing <- c(offsetForEachRing, 360/clustersOnEachRing[i])
  }

  index <- as.numeric(rings)
  count <- 0
  equidistantDegrees <- c()
  for( i in 1:nclust){
    equidistantDegrees <- c(equidistantDegrees, startDegreeForEachRing[index] + count*offsetForEachRing[index])
    index <- index - 1
    if(index < 1)
    {
      index <- as.numeric(rings)
      count <- count + 1
    }
  }

  #Calculate correct angle around rings to place clusters
  #equidistantDegrees <- c(startdegree, sapply(1:(nclust-1), function(x) startdegree+x*degreeoffset))
  #Make sure its converted from degrees to radians
  equidistantDegrees <- sapply(equidistantDegrees, function(x) x*pi/180)

  return(equidistantDegrees)
}

#Calculates and returns the angle extended by the center of the cluster to the center of the network
#We use this angle in Radial Explode Layout Algorithm to calculate new location of centroids after explosion
calculateRadialAngles <- function(nclust, networkCentroidX, networkCentroidY, centroids){
  #Start off with an empty vector
  radialDegrees <- c()

  for(i in 1:nclust)
  {
    #Angle extended by the centroid of the cluster to the centroid of the network, calculated in radians
    angle <- atan2(centroids[i,3] - networkCentroidY, centroids[i,2] - networkCentroidX)

    #create a vector of radial angles
    radialDegrees <- c(radialDegrees, angle)
  }

  return(radialDegrees)
}

#Centroid of each cluster is seen as roughly the center point of the cluster.
#Its calculated by calculating the median of x and y locations
getCentroidForEachCluster <- function(nclust,net, selectAlgoForCentroids){
  centroids <- mat.or.vec(nclust, 3)
  centroids[,1] <- 1:nclust

  #Determine current centroid location of all clusters. Centroids are calculated using medians of X/Y coordinates of members of each cluster.
  for(i in 1:nclust)
  {
    clusterCentroid <- getCentroid(net[net[,2]==i,3], net[net[,2]==i,4], selectAlgoForCentroids)
    centroids[i,2] <- clusterCentroid$X
    centroids[i,3] <- clusterCentroid$Y
  }
  return(centroids)
}

#Extracts coordinates into a matrix
extractCoordinatesIntoMatrix <- function(g, groups, coord){
  net <- cbind(igraph::V(g)$name, groups[match(igraph::V(g)$name, row.names(groups)),1], coord)

  #Min-Max data
  net[,3] <- minmax(net[,3])
  net[,4] <- minmax(net[,4])

  #hard code
  #net[,2] <- 1
  return(net)
}

#Calculates Centroid of the whole network
getCentroid <- function(pointsX, pointsY, selectAlgoForCentroids){

  if(selectAlgoForCentroids == 'median')
  {
    centroidX <- median(as.numeric(pointsX))
    centroidY <- median(as.numeric(pointsY))
  }
  else if(selectAlgoForCentroids == 'mean')
  {
    centroidX <- mean(as.numeric(pointsX))
    centroidY <- mean(as.numeric(pointsY))
  }
  else
  {
    minX = min(as.numeric(pointsX))
    maxX = max(as.numeric(pointsX))

    minY = min(as.numeric(pointsY))
    maxY = max(as.numeric(pointsY))

    centroidX <- (minX+maxX)/2
    centroidY <- (minY+maxY)/2
  }
  returnList <- list("X" = centroidX, "Y" = centroidY)
  return(returnList)
}

###Extrat coords, outcomes, cluster and entity from nodelist file
convertNodeListFile <- function(nodelistFile_name, input){

  nodeListData <-   data.table::fread(nodelistFile_name, sep = "\t", header = TRUE)
  nodeListData <- as.data.frame(nodeListData)

  #if use user coords checked, use FRX and FRY.
  if(input$coordfileCheckBox){
    coordFile <- nodeListData[,1:3]  #coordinates.txt
  }
  else{
    coordFile <- data.frame()
    # if(input$layoutAlgoType == "fr"){
    #   coordFile <- nodeListData[,1:3]  #coordinates.txt
    # }
    # else{
    #   coordFile <- nodeListData[,c(1,4,5)]  #coordinates.txt
    # }
  }


  entityFile <- nodeListData[,c(1,8)] #entity.txt
  moduleFile <- nodeListData[,c(1,7)] #modules
  outcomeFile <- nodeListData[,c(1,6)] #outcomes

  inputList <- list("coords" = coordFile, "outcomes" = outcomeFile, "modules" = moduleFile,  "entity" = entityFile )
  return(inputList)
}


#This method gets user specified project folder and reads all the data from the files that reside within it.
#Method expects to get the network file and module file , coordinates file is optional
getDataFromFile <- function(networkfile_name, modules, rvalues){

  df     <- getDataFrameFromFile(networkfile_name, rvalues)

  # groups <- read.delim(modulefile_name,  sep="\t", row.names=1, header=FALSE)
  groups <- as.matrix(modules[,2])
  row.names(groups) <- modules[,1]
  #rownames(df) <- modules[,1]
  #colnames(df) <- modules[,1]
  net <- getGraphForLayout(df, rvalues)
  g <- net[[1]]
  edgelist <- net[[2]]

  nclust <- length(unique(groups[,1]))

  returnList <- list("df" = df, "nclust" = nclust, "groups" = groups, "g" = g, "edgelist" = edgelist)

  return(returnList)
}


removeOutliers <- function(x) {
  outliers <- NULL
  test <- round(x, digits = 5)
  #print(paste("Length of test before grubbs = ", length(test)))

  result <- outliers::grubbs.test(test, type = 10, two.sided = TRUE)
  pv <- result$p.value
  counter <- 0
  #0.01
  #10^-5
  #0.1
  #
  while(pv > (0.01)) {
    outliers <- as.numeric(strsplit(result$alternative," ")[[1]][3])
    test <- test[!test %in% outliers]
    result <- tryCatch(outliers::grubbs.test(test, type = 10, two.sided = TRUE),
                       warning = function(w){return("NA")},
                       error   = function(e){return("NA")})

    counter <- counter + 1
    if(class(result)== "character" || counter > 1000)
    {
      break
    }
    pv <- result$p.value
  }

  #print(paste("Length of test before grubbs = ", length(test)))
  return(test)
}

# clusterX and clusterY are X and Y co-ordinates of a cluster. Based on these numbers we will find
# the centroid of a cluster and also the distance of each point from the centroid. Centroid is simply the
# mean of clusterX and clusterY
removeOutliersFromEachCluster <- function(clusterX, clusterY){

  #find mean of the cluster
  meanX <- (min(clusterX) + max(clusterX))/2
  meanY <- (min(clusterY) + max(clusterY))/2

  #initialize a vector that holds distance of each point from the centroid
  distFromCentroid <- c()

  for( i in 1:length(clusterX))
  {
    #calculate distance of the point from the centroid
    distance <- dist(rbind(c(meanX, meanY), c(clusterX[i], clusterY[i])))

    #accumulate the value of distance in distFromCentroid vector
    distFromCentroid <- c(distFromCentroid, distance)
  }

  #create a data frame to associate each distance value to the corresponding x,y value
  dataFrame <- data.frame(distFromCentroid, clusterX, clusterY)

  # remove outliers from the distance vector, remain distance values will help us to filter out
  # the X,Y co-ordinates from the
  nonOutliers <- removeOutliers(distFromCentroid)

  # returns a TRUE/FALSE vector denoting if the value is in nonoutlier vector or not
  isAvailable <- distFromCentroid %in% nonOutliers

  #this new dataframe will only contain X,Y values that are non outliers
  newDataFrame <- dataFrame[isAvailable, ]

  newClusterX <- newDataFrame[,2]
  newClusterY <- newDataFrame[,3]

  returnList <- list("clusterX" = clusterX, "clusterY" = clusterY)
  return(returnList)

}

calculateScore <- function(input, net, nclust){

  networkMinX = networkMaxX = networkMinY = networkMaxY <- c()
  clusterList <- list()

  #capture all the clusters without outliers as a list
  for( i in 1: nclust)
  {
    clusterX <- as.numeric(net[net[,2]==i,3])
    clusterY <- as.numeric(net[net[,2]==i,4])

    sizeX <- length(clusterX)
    sizeY <- length(clusterY)

    # cluster <- removeOutliersFromEachCluster(clusterX, clusterY)
    # clusterX <- cluster$clusterX
    # clusterY <- cluster$clusterY

    clusterList <- c(clusterList, list(clusterX, clusterY))
  }


  #Find area of intersection of clusters
  intersectAllClusters <- spatstat::owin(c(0,0), c(0,0))
  unionOfAllClusters <- spatstat::owin(c(0,0), c(0,0))
  for( i in 1: nclust){
    clusterX <- as.numeric(clusterList[[(i-1)*2 + 1]])
    clusterY <- as.numeric(clusterList[[(i-1)*2 + 2]])

    clusteri.minX <- min(clusterX)
    clusteri.maxX <- max(clusterX)
    clusteri.minY <- min(clusterY)
    clusteri.maxY <- max(clusterY)

    networkMinX <- c(networkMinX, clusteri.minX)
    networkMaxX <- c(networkMaxX, clusteri.maxX)
    networkMinY <- c(networkMinY, clusteri.minY)
    networkMaxY <- c(networkMaxY, clusteri.maxY)

    unionOfAllClusters <- spatstat::union.owin(unionOfAllClusters,
                                     spatstat::owin(xrange = c(clusteri.minX, clusteri.maxX),
                                          yrange = c(clusteri.minY, clusteri.maxY)))

    for(j in i:nclust)
    {
      if(i != j)
      {

        clusterX <- as.numeric(clusterList[[(j-1)*2 + 1]])
        clusterY <- as.numeric(clusterList[[(j-1)*2 + 2]])

        clusterj.minX <- min(clusterX)
        clusterj.maxX <- max(clusterX)
        clusterj.minY <- min(clusterY)
        clusterj.maxY <- max(clusterY)

        clusteriWindow <- spatstat::owin(xrange = c(clusteri.minX, clusteri.maxX), yrange = c(clusteri.minY, clusteri.maxY))
        clusterjWindow <- spatstat::owin(xrange = c(clusterj.minX, clusterj.maxX), yrange = c(clusterj.minY, clusterj.maxY))

        intersection <- spatstat::intersect.owin(clusteriWindow, clusterjWindow, fatal = FALSE)
        if(!is.null(intersection))
        {
          intersectAllClusters <- spatstat::union.owin(intersectAllClusters,intersection)
        }
      }
    }
  }

  #Area of union of clusters
  unionArea <- spatstat::area.owin(unionOfAllClusters)
  intersectionArea <- spatstat::area.owin(intersectAllClusters)

  nonOverlapArea <- unionArea - intersectionArea

  rectangle <- spatstat::owin(c(min(networkMinX), max(networkMaxX)), c(min(networkMinY), max(networkMaxY)))
  networkArea <- as.numeric(spatstat::area.owin(w = rectangle))

  score <- nonOverlapArea/networkArea

  # cat(paste("Union area", "Intersect area", "Non-overlap area",  "Score", sep = "\t"), "\n")
  # cat(paste(unionArea, intersectionArea, nonOverlapArea, networkArea, score, sep = "\t"), "\n")

  # print(paste(input$tupleSize, " ", score, sep =""))

  return(score)

}

getAreaofCluster <- function(boundingBox){

  cluster.minX <- boundingBox[1]
  cluster.minY <- boundingBox[2]
  cluster.maxX <- boundingBox[3]
  cluster.maxY <- boundingBox[4]

  # width  <- abs(boundingBox[3] - boundingBox[1])
  # height <- abs(boundingBox[4] - boundingBox[2])
  #
  # return(as.numeric(width*height))

  rectangle <- spatstat::owin(c(cluster.minX, cluster.maxX), c(cluster.minY, cluster.maxY))
  return(spatstat::area.owin(rectangle))
}

getAreaOfOverlap <- function(boundingBox1, boundingBox2){

  cluster1.minX <- boundingBox1[1]
  cluster1.minY <- boundingBox1[2]
  cluster1.maxX <- boundingBox1[3]
  cluster1.maxY <- boundingBox1[4]

  cluster2.minX <- boundingBox2[1]
  cluster2.minY <- boundingBox2[2]
  cluster2.maxX <- boundingBox2[3]
  cluster2.maxY <- boundingBox2[4]

  x_overlap <- max(0, (min(cluster1.maxX, cluster2.maxX) - max(cluster1.minX, cluster2.minX)))
  y_overlap <- max(0, (min(cluster1.maxY, cluster2.maxY) - max(cluster1.minY, cluster2.minY)))
  return(x_overlap*y_overlap)

  # if(cluster1.maxX < cluster2.minX || cluster1.maxY < cluster2.minY || cluster1.minX > cluster2.maxX || cluster1.minY > cluster2.maxY)
  # {
  #   return(0)
  # }
  # else
  # {
  #   #find the area
  #   axis.X <- c(cluster1.minX, cluster1.maxX, cluster2.minX, cluster2.maxX)
  #   axis.Y <- c(cluster1.minY, cluster1.maxY, cluster2.minY, cluster2.maxY)
  #
  #   axis.X <- sort(axis.X)
  #   axis.Y <- sort(axis.Y)
  #
  #   return((axis.X[2] - axis.X[1]) * (axis.Y[2] - axis.Y[1]))
  # }
}

applyRadialExplosion <- function(rings, displace, net, radialDegrees, nclust, centroids, selectAlgoForCentroids) {

  ringassign <- as.numeric(rings)

  networkCentroid <- getCentroid(net[,3], net[,4], selectAlgoForCentroids)
  netX <- networkCentroid$X
  netY <- networkCentroid$Y

  idX <- c()
  idY <- c()

  #Determine final centroid location for each cluster
  for(i in 1:nclust)
  {
    x <- radialDegrees[i]

    idealx <- ringassign*displace*cos(as.numeric(x))
    idealy <- ringassign*displace*sin(as.numeric(x))

    ringassign <- ringassign - 1
    if (ringassign < 1){
      ringassign <- as.numeric(rings)
    }

    shiftx <- idealx- centroids[i,2] - netX
    shifty <- idealy- centroids[i,3] - netY

    net[net[,2]==i,3] <- (as.numeric(net[net[,2]==i,3]) + shiftx)
    net[net[,2]==i,4] <- (as.numeric(net[net[,2]==i,4]) + shifty)

  }

  return(net)
}

applyEquidistantExplosion <- function(net, radialDegrees, nclust, shiftOrRotate, selectAlgoForCentroids, rings, displace) {

  #Map the order of the clusters such that correct cluster goes to the right point on the circle

  mapping <- orderClustersVisually(radialDegrees)

  #Calculate equidistant angles on the circle
  degrees <- calculateEquidistantAngles(nclust, rings)

  networkCentroid <- getCentroid(net[,3], net[,4], selectAlgoForCentroids)
  netX <- networkCentroid$X
  netY <- networkCentroid$Y

  ringassign <- as.numeric(rings)

  #degrees is the equidistant angle that each cluster should make with the circle
  #radialDegree is the angle that cluster currently makes with the centroid
  for(i in 1:nclust)
  {
    diffAngle <- degrees[i] - radialDegrees[mapping[i]]

    clusterCentroid <- getCentroid(net[net[,2]==mapping[i],3], net[net[,2]==mapping[i],4], selectAlgoForCentroids)
    centroidclusterX <- clusterCentroid$X
    centroidclusterY <- clusterCentroid$Y

    #Shift/Rotate whole of the cluster to the new location
    clusterX <- as.numeric(net[net[,2]==mapping[i],3])
    clusterY <- as.numeric(net[net[,2]==mapping[i],4])

    newcentroidX <- ringassign*displace*cos(as.numeric(degrees[i]))
    newcentroidY <- ringassign*displace*sin(as.numeric(degrees[i]))

    shiftx <- newcentroidX - centroidclusterX - netX
    shifty <- newcentroidY - centroidclusterY - netY

    net[net[,2]==mapping[i],3] <- (clusterX + shiftx)
    net[net[,2]==mapping[i],4] <- (clusterY + shifty)

    ringassign <- ringassign -1
    if(ringassign < 1) { ringassign <- as.numeric(rings)}


    if(shiftOrRotate == 'rotate')
    {
      clusterCentroid <- getCentroid(net[net[,2]==mapping[i],3], net[net[,2]==mapping[i],4], selectAlgoForCentroids)

      #Rotation Matrix
      rotm <- matrix(c(cos(diffAngle),sin(diffAngle),-sin(diffAngle),cos(diffAngle)),ncol=2)
      centroidLocationMatrix <- cbind(as.numeric(net[net[,2]==mapping[i],3]),as.numeric(net[net[,2]==mapping[i],4]))
      newClusterLocation <- t( (rotm %*% (t(centroidLocationMatrix) - c(clusterCentroid$X, clusterCentroid$Y))) + c(clusterCentroid$X, clusterCentroid$Y))

      net[net[,2]==mapping[i],3] <- newClusterLocation[,1]
      net[net[,2]==mapping[i],4] <- newClusterLocation[,2]
    }

  }

  return(net)
}

#Clusters are not ordered as they appear visually on the plot, this method maps the visual ordering to cluster ordering in the file
orderClustersVisually <- function(radialAngles){

  len <- length(radialAngles)

  locationArray <- radialAngles

  if(length(locationArray) >1){
  locationArray <- sort(locationArray) #Sort the array in ascending order

  indexNearZero <- which(locationArray >0)[1] #Closest +ve value is the angle that should go to Zero degree on the circle


  # rearrange the array in the neew order
  reOrderedArray <- c(locationArray[indexNearZero:len], locationArray[1:(indexNearZero-1)])

  #Mapping the new order with respect to the radial angles that we have calculated for each cluster
  #we have to map because the cluster numbering does not match with what we have visually avaialble
  for(i in 1:len)
  {
    index <- which(radialAngles == reOrderedArray[i])[1]
    locationArray[i] = index
  }

  return(locationArray)
  }
  else
  {return (1)}
}

applyExplodeLayout <- function(rings, displace, net, selectAlgoForCentroids, shiftOrRotate, layoutAlgoType, nclust, input){

  net[,3] <- minmax(net[,3])
  net[,4] <- minmax(net[,4])


  centroids <- getCentroidForEachCluster(nclust,net, selectAlgoForCentroids)

  # #Calculate the centroid of the whole network
  networkCentroid <- getCentroid(net[,3], net[,4], selectAlgoForCentroids)


  #calculate the angle of each cluster wrt the center of the network for Radial Explosion Layout
  radialDegrees <- calculateRadialAngles(nclust, networkCentroid$X, networkCentroid$Y, centroids)

  net <- applyEquidistantExplosion(net, radialDegrees, nclust, shiftOrRotate, selectAlgoForCentroids, rings, displace)

  if(input$radialOrEquiDist == "radial")
  {
    net <- applyRadialExplosion(rings, displace, net, radialDegrees, nclust, centroids, selectAlgoForCentroids)
  }
  return(net)
}

findRadiusWithMaxGOE  <- function(input, dataFromFile, coord){
  groups <- dataFromFile$groups
  g <- dataFromFile$g
  nclust <- dataFromFile$nclust
  net <- cbind(igraph::V(g)$name, groups[match(igraph::V(g)$name, row.names(groups)),1], coord)

  net <- extractCoordinatesIntoMatrix(g, groups, coord)

  net[,3] <- minmax(net[,3])
  net[,4] <- minmax(net[,4])

  buffer   <- 5  # buffer/vicinity for the location of maxima
  maxScore  =  maxPos = startPos = radius = clusteriScore <- 0
  step     <- 0.1

  while(TRUE){
    tempScore <- maxScore
    endPos <- startPos + buffer*step

    for(i in seq(startPos, endPos, by = step)){
      explode.net <- applyExplodeLayout(input$rings, i, net, input$selectAlgoForCentroids,
                                        input$shiftOrRotate, input$layoutAlgoType, nclust, input)
      clusteriScore <- calculateScore(input, explode.net, nclust)

      if((clusteriScore > maxScore)){
          radius <- i
          maxScore <- clusteriScore
          startPos <- i + step
          break
      }
    }

    if(tempScore == maxScore) { break }
  }
  result <- c(radius, clusteriScore)
  return(result)
}

findClosestRadius <- function(input, dataFromFile, coord){
  groups <- dataFromFile$groups
  g <- dataFromFile$g
  nclust <- dataFromFile$nclust

  net <- cbind(igraph::V(g)$name, groups[match(igraph::V(g)$name, row.names(groups)),1], coord)
 # net[is.na(net)] <- 0

  explode.net <- extractCoordinatesIntoMatrix(g, groups, coord)

  #Get Original Network Score
  originalScore <- calculateScore(input, explode.net, nclust)

  print(paste("Original Score of Layout = ", originalScore, sep = ""))
  #Now calculate score for all the possible circle radius
  minVal <- 0.1
  maxVal <- 2.0
  step   <- 0.1
  pos <- 0
  maxclusteriScore <- 0

  for(i in seq(0.1, 5.0, by = 0.1)){
    explode.net <- applyExplodeLayout(input$rings, i, explode.net, input$selectAlgoForCentroids,
                                      input$shiftOrRotate, input$layoutAlgoType, nclust, input)
    clusteriScore <- calculateScore(input, explode.net, nclust)

        if((clusteriScore - originalScore) > 0){
                radius <- i
                break
        }
   }
  res <- list("radius"= radius, "originalScore" = originalScore)
  return(res)
}

#This method will use the data extracted from the Project folder files and generate the co-ordinates for the explode layout
#Its called each time user changes the sliderinput value and/or textinput value.
getCoordinatesForExplodeLayout <- function(dataFromFile, input, coord, rvalues){
  #Extract graph coordinates into a matrix

  explode.net <- extractCoordinatesIntoMatrix(dataFromFile$g, dataFromFile$groups, coord)

  #to display best network when user navigates from search tab to Visualize tab
  nclust <- length(unique(dataFromFile$groups[,1]))
  #if((!input$showOriginalPlot))
  if(input$viewTypes != "2")
   {   explode.net <- applyExplodeLayout(input$rings, input$displace, explode.net, input$selectAlgoForCentroids,
           input$shiftOrRotate, input$layoutAlgoType, dataFromFile$nclust, input)
   }
  explode.net[,3] <- minmax(explode.net[,3])
  explode.net[,4] <- minmax(explode.net[,4])

  return(explode.net)
}

#This method will save coordinates at different circle radius and no. of circles.
storeCoordinatesForExplodeLayout <- function(dataFromFile, input, coord, rvalues){
  #Extract graph coordinates into a matrix

  explode.net.initial <- extractCoordinatesIntoMatrix(dataFromFile$g, dataFromFile$groups, coord)

  # Store coordinates; these coordinates will be used when the user just loads the data into visualize tab
  #calculate coordinates for original network
  explode.net<- matrix(0,  nrow=nrow(explode.net.initial),  ncol=ncol(explode.net.initial))
  explode.net <- explode.net.initial
  explode.net[,3] <- minmax(explode.net[,3])
  explode.net[,4] <- minmax(explode.net[,4])
  explode.net <- formatExplodeNetMatrix(explode.net, rvalues, input)

  #coordsELFile contains original coords, outcomes, cluster, entity and EL coords from 0.1 to 2.0
  #re-order explode.net to match nodelist format in modim. #Idea is to maintain same file format in both modim and EL.
  if(ncol(explode.net) == 6){
    coordsELFile <- explode.net[,c(1,3,4,6,2,5)] #data with outcomes
  }
  else{
    coordsELFile <- explode.net[,c(1,3,4,2,5)]  #data without outcomes
  }

  #create statList file
  #statFilename <- paste(rvalues$resultsFolder, "/statList.dat", sep= "")
  statFilename <- "statList.dat"
  statFile <- matrix(0,20 , 2)
  colnames(statFile) <- c("Circle-Radius", "CCS")

  #modularityFile <- paste(rvalues$dataFolder,'/', "zscore.txt", sep = "")
  #if(file.exists(modularityFile)){
   # mod <- read.delim(modularityFile, header = TRUE)
  #}

  j <- 1 # number of rings set to by default
  nclust <- length(unique(dataFromFile$groups[,1]))

  if(nclust != 1)
  {
  #calculate coordinates at each circle radius (min=0.1; max=2)
  for(i in seq(0.1, 2, by = 0.1)){

      #initialise variables and matrix
      explode.net<- matrix(0,  nrow=nrow(explode.net.initial),  ncol=ncol(explode.net.initial))
      explode.net <- explode.net.initial

      #ApplyExplode Layout
      explode.net <- applyExplodeLayout(j, i, explode.net.initial, input$selectAlgoForCentroids,
                                        input$shiftOrRotate, input$layoutAlgoType, dataFromFile$nclust, input)
      clusteriScore <- calculateScore(input, explode.net, dataFromFile$nclust)


      if(i >= rvalues$ClosestRadius){

        statFile[(i*10),1] <- as.numeric(i)
        statFile[(i*10),2] <- as.numeric(sprintf("%.4f",as.numeric(clusteriScore)))
        #statFile[(i*10),3] <- 2.5
        #statFile[(i*10),4] <- mod[1,1]
        #statFile[(i*10),5] <- mod[1,2]
      }
      else{
        statFile[(i*10),2] <- ""
        #statFile[(i*10),3] <- ""
        #statFile[(i*10),4] <- ""
        #statFile[(i*10),5] <- ""
      }

      explode.net[,3] <- minmax(explode.net[,3])
      explode.net[,4] <- minmax(explode.net[,4])
      explode.net <- formatExplodeNetMatrix(explode.net, rvalues, input)
      ELcoords <- explode.net[, 3:4]
      if(input$layoutAlgoType == "fr"){
          colnames(ELcoords) <- c(paste("EL.FRX", i, sep =""), paste("EL.FRY", i, sep =""))
      }
      else{
          colnames(ELcoords) <- c(paste("EL.KKX", i, sep =""), paste("EL.KKY", i, sep =""))
      }

      coordsELFile <- cbind(coordsELFile, ELcoords)
  }

    statFile <- statFile[statFile[,1]!=0,]
    write.table(statFile[,c(1,2)],statFilename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )

  }
  #filename <- paste(rvalues$resultsFolder, "/coordinates.dat", sep= "")
  #filename <- "coordinates.dat"
  #write.table(coordsELFile, filename, sep = "\t", quote = FALSE, row.names = FALSE )
  coordsELFile <- data.frame(coordsELFile)
  coordsELFile[,1] <- as.character(coordsELFile[,1])
  coordsELFile[,2:length(coordsELFile)] <- lapply(coordsELFile[,2 : length(coordsELFile)], function(x) as.numeric(as.character(x)))
  return (coordsELFile)
}

#Updates the TextInput field on the File Popup after user selects the project folder
showFilesToUser <- function(session, userSelectedDataFolder, userSelectedDataFolderFullPath, rvalues){
  updateTextInput(session, 'projectFolder', value=paste(userSelectedDataFolder, '/Data/'))
  updateTextInput(session, 'networkFile', value=paste(userSelectedDataFolder, '/Data/', userSelectedDataFolder, '_subset.txt',sep = ""))
  updateTextInput(session, 'nodeListFile', value=paste(userSelectedDataFolder, '/Data/', userSelectedDataFolder, '_nodelist.dat',sep = ""))

  nodelist <- read.delim(paste(getwd(),'/Project/',userSelectedDataFolder, '/Data/', userSelectedDataFolder, '_nodelist.dat',sep = ""), sep = "\t", header = TRUE)
  if(is.na(nodelist[,2:3])){
    updateCheckboxInput(session, "coordfileCheckBox", value = FALSE)
    shinyjs::disable("coordfileCheckBox")
  }
  else{
    shinyjs::enable("coordfileCheckBox")
  }
  #coordFile <- paste(userSelectedDataFolderFullPath,'/Data/', userSelectedDataFolder, '_subset_coords.txt',sep = "")

  # if(file.exists(coordFile)){
  #   updateTextInput(session, 'coordFile', value=paste(userSelectedDataFolder, '/Data/',userSelectedDataFolder, '_subset_coords.txt',sep = ""))
  #   updateCheckboxInput(session, "coordfileCheckBox", value = TRUE)
  #   shinyjs::enable("coordfileCheckBox")
  # }
  # else{
  #   updateTextInput(session, 'coordFile', value="No File")
  #   updateCheckboxInput(session, "coordfileCheckBox", value = FALSE)
  #   shinyjs::disable("coordfileCheckBox")
  # }
  #
  # outcomeFile <-   paste(userSelectedDataFolderFullPath,'/Data/', userSelectedDataFolder, '_subset_outcomes.txt',sep = "")
  # if(file.exists(outcomeFile)){
  #   updateTextInput(session, 'outcomeFile', value=paste(userSelectedDataFolder, '/Data/',userSelectedDataFolder, '_subset_outcomes.txt',sep = ""))
  #   }
  # else{
  #   updateTextInput(session, 'outcomeFile', value="No File")
  #   }
}

#This method is the wrapper method that calls getDataFromFile. server.R calls this method to get data extracted from the files
getData <- function(basePath, input, rvalues, modules){

   networkFile <- paste(basePath, '/',input$networkFile, sep = "")
   #moduleFile <- paste(basePath, '/',input$moduleFile, sep = "")

   dataFromFile <- getDataFromFile(networkFile, modules, rvalues)

  return(dataFromFile)
}

getEdgeList <- function (path){

  readEdgeList <- data.table::fread(path , sep ="\t")
  readEdgeList <- as.data.frame(readEdgeList)
  return(readEdgeList)
}

plotnetwork <- function(input, sliderValues, rvalues, libUsedForPlot, labelsData= NULL){

  if(libUsedForPlot == 3)
  {
    #IGRAPH PLOT
    p <- plot(rvalues$dataFromFile$g, layout=matrix(as.numeric(sliderValues[,3:4]), nrow(sliderValues), 2, byrow=F),
         vertex.label = c(rep(NA, rvalues$nrows), rvalues$colnames),
         vertex.label.color = "black",
         vertex.size=input$vertexSize,
         vertex.label.font =2,
         vertex.color = rvalues$vcolors,
         vertex.label.cex=1.5,
         vertex.label.dist = 0.25,
         vertex.label.degree = pi/4,
         vertex.label.family = "TT Ariel",
         vertex.shape = rvalues$vertexShape, xlim=c(-1, 1), ylim=c(-1, 1)
    )
  }
  else
  {
    #library(igraph,  warn.conflicts = FALSE, quietly = TRUE)

    #plotcord <- as.data.frame(matrix(as.numeric(sliderValues[,2:3]), nrow(sliderValues), 2, byrow=F))
    plotcord <- as.data.frame(sliderValues[,2:3])
    colnames(plotcord) = c("X1", "X2")

    #edgeListpath <- paste(rvalues$resultsFolder, '/edgeList.dat', sep="")
    #elist <- getEdgeList(edgeListpath)
    elist <- rvalues$dataFromFile$edgelist

    #net <- asNetwork(rvalues$dataFromFile$g)
    net <- intergraph::asNetwork(igraph::graph.data.frame(elist, directed=FALSE))

    #detach("package:igraph")

    edglist <- network::as.matrix.network.edgelist(net, attrname = 'Weight')
    if (rvalues$netType == 1) {
    	edges <- data.frame(plotcord[elist[, 1], ], plotcord[elist[, 2], ])

    } else {
    	edges <- data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
    }


    #plotcord$elements <- as.factor(get.vertex.attribute(net, "elements"))

    #set vertex color, shape  and stroke
    rvalues$vcolors <- getVertexColor(sliderValues, input, rvalues)
    rvalues$vertexShape <- getVertexShape(sliderValues, input, rvalues)

    colnames(edges) <- c("X1", "Y1", "X2", "Y2")

    rvalues$vertexSize <- setVertexSize(sliderValues, input, rvalues, edglist)
    rvalues$edgeThickness <- setEdgeThickness(sliderValues, input, rvalues, edglist)
    rvalues$stroke <- setStrokeSize(sliderValues, input, rvalues)

    if(libUsedForPlot == 1)
    {
      p <- plotUsingGGplot(plotcord, edges, rvalues, input, labelsData)
    }
    else if(libUsedForPlot == 2)
    {
      p <- plotUsingGraphics(plotcord, edges, rvalues, input)
    }
    else
    {
      p <- forceNetwork(Links = edges, Nodes = plotcord, NodeID = "name")
    }

    #library(igraph,  warn.conflicts = FALSE, quietly = TRUE)
  }

  return(p)

}
setVertexSize <-function(sliderValues, input, rvalues, edgeList) {

  nodesize <- getVertexSize(input$vertexSize, session, input)

  ##calculate weight of edges each all the nodes.
  if((input$nodeDegree) && (nodesize !=0)){

    nodeInfo <- sliderValues
    edgeInfo <- edgeList

    newEdges <- as.matrix(data.frame(nodeInfo[edgeInfo[, 1], c(1,6) ], nodeInfo[edgeInfo[, 2], c(1,4,6) ]))
    colnames(newEdges) <- c('Label', 'Entity', 'Label2', 'EdgeWeight', 'Entity2')
    newEdges[,4] <- edgeInfo[,3]

    weightedNodes <- nodeInfo[, c(1,5,6,6)] #1,2,5,5
    colnames(weightedNodes) <- c('Label', 'WeightedNodeDegree', 'Entity', 'VertexSize')
    weightedNodes[,c(2,4)] <- 0
    weightedNodes[weightedNodes[,3]==2,4]<-nodesize

    for(i in 1:nrow(weightedNodes)){
      if(weightedNodes[i,3] == 1){
          weightedNodes[i,2]  <- sum(as.numeric(newEdges[weightedNodes[i,1] == newEdges[,1],4]))
      }
    }

    minweightedNodeDegree <- min(weightedNodes[weightedNodes[,3]==1,2])

    for(i in 1:nrow(weightedNodes)){
      if(weightedNodes[i,2] == minweightedNodeDegree){#weighted node degree equals min weighted node degree
        weightedNodes[i,4]  <- nodesize
      }
      else{
        weightedNodes[i,4]  <- nodesize + weightedNodes[i,2]
      }
    }

    vertexSize <- weightedNodes[,4]
  }
  else{
    vertexSize <- nodesize
  }
  return(vertexSize)

}

setEdgeThickness <- function(sliderValues, input, rvalues, edgeList) {

  edgeThickness <- getThickness(input)

  if((input$edgeWeight) && (edgeThickness!=0)){

    nodeInfo <- sliderValues
    edgeInfo <- edgeList

    newEdges <- as.matrix(data.frame(nodeInfo[edgeInfo[, 1], c(1,6) ], nodeInfo[edgeInfo[, 2], c(1,6,4,6) ]))
    #newEdges <- as.matrix(data.frame(nodeInfo[edgeInfo[, 1], c(1,8)], nodeInfo[edgeInfo[, 2], c(1,8,6,7)]))

    colnames(newEdges) <- c('Label', 'Entity', 'Label2', 'Entity2', 'EdgeWeight', 'EdgeThickness')
    newEdges[,6] <- 0
    newEdges[,5] <- edgeInfo[,3] #*edgeWeight

    minEdgeWeight <- min(as.numeric(newEdges[,5]))

    for(i in 1:nrow(newEdges)){
      if(newEdges[i,5] == minEdgeWeight){#weighted node degree equals min weighted node degree
        newEdges[i,6]  <- edgeThickness
      }
      else{
        newEdges[i,6]  <- edgeThickness * as.numeric( newEdges[i,5])
      }
    }

    edgeThickness <- as.numeric(newEdges[,6])
  }


  return(edgeThickness)
}

setStrokeSize <- function(sliderValues, input, rvalues){

  if((!is.null(rvalues$clusterNum)) && (rvalues$clusterNum !=0)){

    stroke <- 0
    strokeMat <- cbind(sliderValues, stroke)
    #for(i in 1:length(rvalues$clusterNum)){
     strokeMat[strokeMat[,5] == rvalues$clusterNum[1], 'stroke'] <- 1
    #}

    strokeMat[strokeMat[,'stroke']==0, 'stroke'] <- 0.3
    stroke <- strokeMat[,'stroke']
  }
  else{
    stroke = 0.3
  }
  return(stroke)
}

plotUsingGGplot <- function(plotcord, edges, rvalues, input, labelsData = NULL){

  p <- ggplot2::ggplot()

  if(rvalues$edgeThickness!=0){
      p <- p + ggplot2::geom_segment(ggplot2::aes(x = X1, y = Y1, xend = X2, yend = Y2), data = edges, size = rvalues$edgeThickness, colour = "#888888")
  }

  if(rvalues$vertexSize!=0){
    p <- p + ggplot2::geom_point(ggplot2::aes(X1, X2),colour = "black",  fill = rvalues$vcolors, shape = rvalues$vertexShape, size = rvalues$vertexSize, data = plotcord, stroke = rvalues$stroke)  #add shape, vertexcolor: original
  }

  p <-p + ggplot2::scale_colour_brewer(palette = "Set1")  + ggplot2::scale_x_continuous(breaks = NULL) + ggplot2::scale_y_continuous(breaks = NULL)
  p <-p + ggplot2::theme(panel.background = ggplot2::element_blank()) + ggplot2::theme(legend.position = "none")
  p <-p + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank())
  # p <-p + theme(legend.background = element_rect(colour = NA))+ theme(panel.background = element_rect(fill = "white", colour = NA))
  # p <-p + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  #p <- p + geom_text(aes(plotcord$X1, plotcord$X2,label=c(rep(NA, rvalues$nrows), rvalues$colnames)),
  #                   hjust=0,vjust=0, size = 6, fontface = "bold", check_overlap = TRUE )

  return(p)
}

getVertexColor <- function(explode.net, input, rvalues){
  #explodeNet is a matrix with colnames: "PID", "Cluster", "X1", "Y1", "Entity", "Outcome".
  entityColor.index  <- c("red", "black")
  #838383
  #C9197A - dark pink ; #CE9A89 - pale brown; #0AF3EE - aqua/pale blue;
  #94EA18 - green ; #E5F115 - yellow; #8825F9 - dark violet; #d95ea1 - light pink; #6b8e23 - cow green, #C0C4E0-  light grey
  clusterColor.index <- c("#C9197A", "#94EA18", "#0AF3EE", "#808080", "#E5F115", "#8825F9", "#CE9A89",  "#fd991c", "black", "#6b8e23", "#C0C4E0", "#669999", "blue", "pink")  #distinctColorPalette(25)
  clusterColor.index[20] <-c("black")

  outcomeColor.index <- c("magenta", "sky blue", "orange", "yellow", "brown")
  #outcomeColor.index[20] <- c("black")

  if (ncol(explode.net) == 5){
    vcolors <- switch(input$nodeColor,
                              #"1"=sapply(as.data.frame(explode.net)[,6], function(x) outcomeColor.index[x]),
                              "2"=sapply(as.data.frame(explode.net)[,4], function(x) clusterColor.index[x]),
                              "3"=sapply(as.data.frame(explode.net)[,5],  function(x) entityColor.index[x])
    )
  }
  else{
    vcolors <- switch(input$nodeColor,
                              "1"=sapply(as.data.frame(explode.net)[,4], function(x) { if(x != 20){ return (outcomeColor.index[x])}
                                else {return ("black")}

                              }
                              ),
                              "2"=sapply(as.data.frame(explode.net)[,5], function(x) clusterColor.index[x]),
                              "3"=sapply(as.data.frame(explode.net)[,6],  function(x) entityColor.index[x])
    )


  }

  return (vcolors)

}

getVertexShape <- function(explode.net, input, rvalues){

  entity <- as.numeric(explode.net[,6])
  shape.index = c(21, 24)
  vertexShape <- sapply(entity, function(x) shape.index[x])

  return(vertexShape)

}

getDataInBoundingBox <- function(data , ranges, rvalues){
    #data <- newSliderValues$explode.net
  if(!(is.null(ranges$x)) && (!is.null(ranges$y))){

      rvalues$nrows <- nrow(data[data[,6] == 1,] )
      if (rvalues$netType == 1) {
        all_var <- as.character(data[data[,6] == 1,1])
      } else {
        all_var <- as.character(data[data[,6] == 2,1])
      }
       # list of all variables
      names <- colnames(data)
  
      #retrieve data points inside bounding box
      a <- data[data[,names[2]] >= ranges$x[1], ]
      b <- a[a[,names[2]] <= ranges$x[2], ]
      c <- b[b[,names[3]] >= ranges$y[1], ]
      d <- c[c[,names[3]] <= ranges$y[2], ]

      if (rvalues$netType == 1) {
        varname <- d[d[,6] == 1, 1]
        } else {
          varname <- d[d[,6] == 2,1] #var name inside bounding box.

        }
  
      rvalues$colnames <- as.character(varname[match(all_var, varname)])
  }
  else{
      #no. of rows in network file
      if (rvalues$netType == 1) {
        rvalues$nrows <- nrow(data)
        rvalues$colnames <- as.character(data[,1])  #retrieve variable names.
      } else {
        rvalues$nrows <- nrow(data[data[,6] == 1,] ) #counts number of entities.
        rvalues$colnames <- as.character(data[data[,6] == 2,1])  #retrieve variable names.
      }
  }

  labelsData <- list('nrows' = rvalues$nrows, 'colnames' = rvalues$colnames)
  return(labelsData)
}
#plotUsingGraphics <- function(plotcord, edges, rvalues, input){
  # plot(x=NULL, y=NULL, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05),  xaxt = "n", yaxt = "n", axes = FALSE, xlab="", ylab="")
  #
  # #Adding the edges to the plot
  # segments(x0 = edges$X1, y0 = edges$Y1 , x1 = edges$X2, y1 = edges$Y2, col = "grey" )
  #
  # # #Plot the nodes
  # points(x = as.numeric(plotcord$X1), y = as.numeric(plotcord$X2), col = "black",
  #        cex = input$vertexSize/2, pch = rvalues$vertexShape1, bg = rvalues$vcolors, xaxt = "n", yaxt = "n", xlab="", ylab="")
  #
  # text(x = as.numeric(plotcord$X1), y = as.numeric(plotcord$X2), labels=c(rep(NA, rvalues$nrows), rvalues$colnames), cex= 1.5)
#}
#=========================

# gggraph <- function(g, xpos, ypos) {
#
#   require(ggplot2)
#
#   g_ <- get.edgelist(g)
#   g_df <- as.data.frame(g_)
#   g_df$id <- 1:length(g_df[,1])
#   g_df <- melt(g_df, id=3)
#   xy_s <- data.frame(value = unique(g_df$value), x = xpos, y = ypos)
#                      # x = vplace(length(unique(g_df$value))),
#                      # y = vplace(length(unique(g_df$value))))
#   g_df2 <- merge(g_df, xy_s, by = "value")
#
#   p <- ggplot(g_df2, aes(xpos, ypos)) +
#     geom_point()
#   # +
#   #   geom_line(size = 0.3, aes(group = id, linetype = id)) +
#   #   geom_text(size = 3, hjust = 1.5, aes(label = value)) +
#   #   theme_bw() +
#   #   opts(panel.grid.major=theme_blank(),
#   #        panel.grid.minor=theme_blank(),
#   #        axis.text.x=theme_blank(),
#   #        axis.text.y=theme_blank(),
#   #        axis.title.x=theme_blank(),
#   #        axis.title.y=theme_blank(),
#   #        axis.ticks=theme_blank(),
#   #        panel.border=theme_blank(),
#   #        legend.position="none")
#
#   p
#
# }

# ggbigraph <- function(g) {
#
#   require(ggplot2)
#
#   g_ <- get.edgelist(g)
#   g_df <- as.data.frame(g_)
#   g_df$id <- 1:length(g_df[,1])
#   g_df <- melt(g_df, id=3)
#   xy_s <- data.frame(value = unique(g_df$value),
#                      x = c(rep(2, length(unique(g_df$value))/2), rep(4, length(unique(g_df$value))/2)),
#                      y = rep(seq(1, length(unique(g_df$value))/2, 1), 2))
#   g_df2 <- merge(g_df, xy_s, by = "value")
#
#   p <- ggplot(g_df2, aes(x, y))
  # +
  #   geom_point() +
  #   geom_line(size = 0.3, aes(group = id, linetype = id)) +
  #   geom_text(size = 3, hjust = 1.5, aes(label = value)) +
  #   theme_bw()
  # +
  #   opts(panel.grid.major=theme_blank(),
  #        panel.grid.minor=theme_blank(),
  #        axis.text.x=theme_blank(),
  #        axis.text.y=theme_blank(),
  #        axis.title.x=theme_blank(),
  #        axis.title.y=theme_blank(),
  #        axis.ticks=theme_blank(),
  #        panel.border=theme_blank(),
  #        legend.position="none")

 # p

#}

# plotg <- function(net, value = NULL) {
#   m <- as.matrix.network.adjacency(net)  # get sociomatrix
#
#   # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
#   plotcord <- data.frame(gplot.layout.fruchtermanreingold(m, NULL))
#
#   colnames(plotcord) = c("X1", "X2")
#   edglist <- as.matrix.network.edgelist(net)
#   edges <- data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
#   plotcord$elements <- as.factor(get.vertex.attribute(net, "elements"))
#   colnames(edges) <- c("X1", "Y1", "X2", "Y2")
#
#   edges$midX <- (edges$X1 + edges$X2)/2
#   edges$midY <- (edges$Y1 + edges$Y2)/2
#
#   pnet <- ggplot()
#     + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2), data = edges, size = 0.5, colour = "grey")
#     + geom_point(aes(X1, X2, colour = elements), data = plotcord)
#     + scale_colour_brewer(palette = "Set1")
#     + scale_x_continuous(breaks = NULL)
#     + scale_y_continuous(breaks = NULL)
#     + theme(panel.background = element_blank()) + theme(legend.position = "none") # discard default grid + titles in ggplot2
#     + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
#     + theme(legend.background = element_rect(colour = NA))
#     + theme(panel.background = element_rect(fill = "white", colour = NA))
#     + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
#   return(print(pnet))
# }
