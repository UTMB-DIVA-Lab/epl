#'The 'search' step of the ExplodeLayout(EL) algorithm.
#'
#'This function basically takes as input the network and cluster information and generate the
#'networks which are exploded using the EL algorithm.
#'
#'Input:
#'There are mainly two data sets the user needs to have ready before starting:
#'1. A 'nodelist' which containing the lables, original coordinates,
#'outcome (like case or control), entity, and cluster membership for all
#'nodes in the network.
#'2. A 'network' which represents the network.
#'
#'The description of the data format:
#'The nodelist must be a dataframe with 8 columns, named as 'Label', 'FRX', 'FRY'
#''KKX', 'KKY', 'Outcome', 'Cluster', 'Entity'.
#'The network must be a dataframe or matrix. When it's a bipartite network, the matrix should be a
#'incidence matrix, and the column names should be the variable names and row names should be the subjects
#'names. When it's unipartite, the matrix should be an adjacent matrix, and the column names and row names
#'are the same. The first column of the network file should be the outcome of the nodes.
#'
#'Please note: the vertices names in the network MUST be consistent with the ones in the network!
#'
#'As the two data sets are not trivial to make manually, it's highly recommended that
#'user use the function in our package to generate them.
#'
#'Output:
#'The search function will return a list containing two dataframes. One called 'coordinates'
#'mainly contains the coordinates of all vertices in networks with different radius ranging from 0.1 to 2.0
#'(radius is a parameter in EL algorith). Another called 'edgelist' is the edgelist of the input network. These two dataframes
#'would serve as the input of the next visuallization step.
#'
#'By default, the serach function would save the 'coordinates' and 'edgelist' dataframes into two local files.
#'If the user prefers just save them in memory or output as .Rdata files, he/she can choose to turn off the
#'output parameter.
#'
#'The calculated CCS scores with different network radius are also saved into a local 'statList.dat' file.
#'Users can refer to the CCS scores when trying to determine which graph(s) they'd like to visualize in
#'the next step.
#'
#'@param nodelistFile  a nodelist object, which is a dataframe with first column as the name of nodes
#'@param networkFile a network object, which is a matrix or dataframe with first column as the outcome
#'@param projName a string which is the project name or data name specified by user. This name will be used in saved export files.
#'@param coordfileCheckBox a boolean variable which specifies whether the nodelistFile contains initial coordinates for all vertices.
#'If the nodelist is generated using our functions, it will contain two sets of original coordinates which are resulted from the
#'Fruchterman-Reingold and Kamada-Kawai layout.
#'@param initLayout a string which can be set as 'fr' or 'kk' when you have Fruchterman-Reingold and Kamada-Kawai layout initial coordinates.
#'@param radialOrEquiDist a string which is either 'equidist' or 'radial', which is a parameter used in the EL algorithm. For details please
#'read the reference paper below.
#'@param shiftOrRotate a string which is either 'shift' or 'rotate'. This is also a parameter used in EL algorithm.
#'@param selectAlgoForCentroids a string which can be 'mean', 'median' or 'minmax'. This is a parameter in EL algorithm to determine
#'the centroids of clusters and networks.
#'@param rings a interger variable which is the number of imaginary circles used in EL. It's normally set as 1 unless user has too many clusters.
#'@param output a boolean variable to specify whether save the output of the function into local files.
#'
#'@references \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5543384/}
#'
#'@export
search <- function(nodelistFile, networkFile, projName = 'defaultName', coordfileCheckBox = TRUE,
                   initLayout = 'fr', radialOrEquiDist = 'equidist', shiftOrRotate = 'rotate',
                   selectAlgoForCentroids = 'median', rings = 1, output = TRUE) {
  #setting default parameters for explodeLayout
  input <- list()
  rvalues <- list()

  #true means there's coordinates file to input
  input$coordfileCheckBox <- coordfileCheckBox
  #choose layout algorithm if no initial coodinates are read in
  input$layoutAlgoType <- initLayout
  #choose equil distance centroids on the explodeLayout circle
  input$radialOrEquiDist <- radialOrEquiDist
  #number of circles for explodeLayout
  input$rings <- rings
  #choose how to decide centroids for each cluster
  input$selectAlgoForCentroids <- selectAlgoForCentroids
  #choose whether to shift or rotate clusters
  input$shiftOrRotate <- shiftOrRotate

  #choose view types (2 means the standard explodeLayout, other views would be added later)
  input$viewTypes <- '2' #ExplodeLayout View
  input$displace <- 0.4  #???


  nodelist <- convertNodeListFile(nodelistFile, input)
  rvalues$netType <- length(unique(nodelist$entity$Entity))
  if (rvalues$netType > 2) {stop("Data not supported. Make sure your input has only one or two entities.")}

  #get data from file, edgelist.dat file is saved
  #rvalues$dataFromFile <- getData(basePath, input, rvalues, nodelist$modules)
  rvalues$dataFromFile <- getDataFromFile(networkFile, nodelist$modules, rvalues)

  rvalues$entity <- nodelist$entity    #as.matrix(entity)

  #rvalues$colnames <- colnames(rvalues$dataFromFile$df)
  rvalues$colnames <- nodelist$coords[,1]

  rvalues$nrows <- nrow(rvalues$dataFromFile$df)

  rvalues$outcome <- nodelist$outcomes

  rvalues$coord <- getOriginalNetworkCoordinates(input, rvalues$dataFromFile$g, nodelist$coords)

  res <- findClosestRadius(input, rvalues$dataFromFile, rvalues$coord)
  radius <- res$radius

  rvalues$originalScore <- res$originalScore
  rvalues$ClosestRadius <- radius

  GOEresult <- findRadiusWithMaxGOE(input, rvalues$dataFromFile, rvalues$coord)
  maxRadius <- GOEresult[1]
  rvalues$displace <- maxRadius

  goodnessscore <-  as.numeric(sprintf("%.4f",as.numeric(GOEresult[2])))

  #calculate coordinates for networks using ExplodeLayout at different radius, then store them
  #coordinates.dat and statList.dat files are saved
  returnExplodeNetwork <- storeCoordinatesForExplodeLayout(rvalues$dataFromFile, input, rvalues$coord, rvalues)

  explode.net <- getCoordinatesForExplodeLayout(rvalues$dataFromFile, input, rvalues$coord, rvalues)

  if(rvalues$dataFromFile$nclust != 1){
    explodeGoodness <- calculateScore(input, explode.net, rvalues$dataFromFile$nclust)
  } else{
    explodeGoodness <- "NA"
  }

  explode.net <- formatExplodeNetMatrix(explode.net, rvalues, input)
  #cat(paste(input$displace, ";", explodeGoodness,"\n", sep=""))

  edgelist <- rvalues$dataFromFile$edgelist
  colnames(edgelist) <- c('V1', 'V2', 'Weight')
  write.table(edgelist, paste(projName,'edgelist.dat', sep = '_'), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(returnExplodeNetwork, paste(projName,'coordinates.dat', sep = '_'), sep = "\t", quote = FALSE, row.names = FALSE )


  ex_output <- list()
  ex_output$coordinates <- returnExplodeNetwork
  ex_output$edgelist <- edgelist
  return (ex_output)
}
