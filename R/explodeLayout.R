#'Layout a network using the ExplodeLayout Algorithm.
#'
#'The explodelayout function applies the ExplodeLayout Algorithm to help user visuallize
#'and analyze their network data.
#'
#'
#'The explodelayout function takes as input two main files:
#'1. A nodelist file which containing the lables, original coordinates,
#'outcome (like case or control), entity, and cluster membership for all
#'nodes in the network.
#'2. A network file which represents the network.
#'
#'Data Format:
#'The nodelist file must be a dataframe with 8 columns, named as 'Label', 'FRX', 'FRY'
#''KKX', 'KKY', 'Outcome', 'Cluster', 'Entity'
#'The network file must be a dataframe or matrix. When it's a bipartite network, the column names
#'should be the variable names and row names should be the subjects names. When it's unipartite,
#'column names and row names are the same. The first column of the network file should be the outcome
#'of the nodes.
#'
#'Please note: the names in the network files MUST be consistent with the names in the network files.
#'
#'Output: the function will export the layout of network with different radius into 21 png files, which
#'includes 1 original layout graph and 20 layout graphs with radius of 0.1, 0.2, ... 1.9, 2.0. The
#'corresponding CCS scores of different radius will be saved into a 'statList.dat' file, which would
#'enable user to pick the best CCS score graph.

#'@param nodelistFile  a nodelist, which is a dataframe with first column as the name of nodes
#'@param networkFile a network, which is a matrix or dataframe with first column as the outcome
#'@param data_name name your data so that we know how to name your results files
#'@param labelOn labels of nodes on or not. should be 'yes' or 'no'
#'@examples
#'explodelayout(nl1, net1, 'no','sample_data')
#'
#'@references \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5543384/}
#'
#'
#'
#'@export
explodelayout <- function(nodelistFile, networkFile, labelOn = 'no', data_name = 'data') {
#para <- read.table('explodeLayout_input.txt', header = F, stringsAsFactors = F)

#data_name <- para[1,]
#labelOn <- para[2,]

################ DO NOT modify anything below this line ###############################
#######################################################################################

#devtools::load_all()

#import library functions written in separated scripts
#libPath1 <- paste(getwd(), 'ExplodeLayout_Library','libraries.R', sep = '/')
#libPath2 <- paste(getwd(), 'ExplodeLayout_Library','functions.R', sep = '/')
#libPath3 <- paste(getwd(), 'ExplodeLayout_Library','explodeLayoutMethodLibrary_uni&bi.R', sep = '/')
#source(libPath1)
#source(libPath2)
#source(libPath3)


#setting default parameters for explodeLayout
input <- list()
rvalues <- list()


#true means there's coordinates file to input
input$coordfileCheckBox <- TRUE
#choose layout algorithm if no initial coodinates are read in
input$layoutAlgoType <- 'fr'
#choose equil distance centroids on the explodeLayout circle
input$radialOrEquiDist <- "equidist"
#number of circles for explodeLayout
input$rings <-1
#choose how to decide centroids for each cluster
input$selectAlgoForCentroids <- 'median'
#choose whether to shift or rotate clusters
input$shiftOrRotate <- "rotate"
#choose view types (2 means the standard explodeLayout, other views would be added later)
input$viewTypes <- '2' #ExplodeLayout View
input$displace <- 0.4  #???
#vertex size
input$vertexSize <- 3
#vertex label size
input$labelSize <- 4
#show node degree or not
input$nodeDegree <- FALSE
#show edge weight or not
input$edgeWeight <- FALSE
#choose node color mode (color by clusters, outcome, etc..)
input$nodeColor <- '2' #cluster color
#set default edge thickness or not
input$edgeThickness <- 0.1
#select the picture format to export
input$selectExportFileType <- 'png'



#specify a few path for input and output files
#basePath   <- c(wd=paste(getwd(),'/Project', sep = ''))
#exportPath <- c(wd=normalizePath("~"))


#userSelectedDataFolder <- data_name
#userSelectedDataFolderFullPath <- paste(basePath, '/', userSelectedDataFolder, sep = "")

#dataFolder <- paste(userSelectedDataFolderFullPath,'/','Data',sep = "")
#rvalues$dataFolder  <- dataFolder
#rvalues$resultsFolder <- paste(userSelectedDataFolderFullPath,'/','Results',sep = "")


#import data files provided by users
#input$nodeListFile <- paste(userSelectedDataFolder, '/Data/', userSelectedDataFolder, '_nodelist.dat',sep = "")
#input$networkFile <- paste(userSelectedDataFolder, '/Data/', userSelectedDataFolder, '_subset.txt',sep = "")

#nodelistFile <- paste(basePath, '/',input$nodeListFile, sep = "")
#nodelist_name <- paste(data_name, '_nodelist.dat',sep = "")
#network_name <- paste(data_name, '_subset.txt',sep = "")
#nodelistFile <- system.file("extdata", nodelist_name, package = "epl")  #full path of nodelist
#networkFile <- system.file("extdata", network_name, package = "epl")  #full path of network

#put nodelistFile and networkFile as two parameters

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



######################################## generating graphs  ###############################################
ranges <- list()


#input$visualDirectorypath <- rvalues$resultsFolder
#if (file.exists('coordinates.data')) {
#returnExplodeNetwork <- data.frame(data.table::fread('coordinates.dat', sep = "\t"))}
columnNames <- colnames(returnExplodeNetwork)


coordX <- columnNames[2]
coordY <- columnNames[3]

explodeNet <- cbind.data.frame(returnExplodeNetwork[,1], returnExplodeNetwork[,colnames(returnExplodeNetwork) == coordX],
                               returnExplodeNetwork[,colnames(returnExplodeNetwork) == coordY], returnExplodeNetwork[,4],
                               returnExplodeNetwork[,5], returnExplodeNetwork[,6])

colnames(explodeNet) <- c(columnNames[1], coordX, coordY, columnNames[4], columnNames[5], columnNames[6] )

returnList <- list("explode.net" = explodeNet)#, "edgeList" = edgeList)

newSliderValues <- returnList

explode.net <- newSliderValues$explode.net

ranges$x <- c(min(explode.net[,2]), max(explode.net[,2]))
ranges$y <- c(min(explode.net[,3]), max(explode.net[,3]))

labelsData <- getDataInBoundingBox(newSliderValues$explode.net , ranges, rvalues)
lab <- labelsData$colnames
if(rvalues$netType == 2) {
  lab <- c(rep(NA, labelsData$nrows), labelsData$colnames)
}

p <-  suppressWarnings(plotnetwork(input, newSliderValues$explode.net, rvalues, 1, labelsData))

plotcord <- as.data.frame(newSliderValues$explode.net[,2:3])
colnames(plotcord) = c("X1", "X2")



# generate and save the original graph (before applying the ExplodeLayout algorithm)
p <- p + ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y)

if(labelOn == 'yes' && ((input$labelSize != 0) && (!is.na(input$labelSize)))){
  p <- p + suppressWarnings( #rvalues$colnames
    ggrepel::geom_label_repel(ggplot2::aes(plotcord$X1, plotcord$X2,label=lab, alpha=0.5),
                     na.rm = TRUE,
                     size = input$labelSize,
                     fontface = "bold" ,
                     box.padding = ggplot2::unit(0.5, "lines"),
                     point.padding = ggplot2::unit(0.5, "lines"),
                     segment.color = "#0000FF"
                     # remove lollipop stick : segment.size = 0
    )
  )
}

rvalues$plot <- p
exportFile <- paste(data_name, '_original','.',input$selectExportFileType, sep="")
ggplot2::ggsave(exportFile, rvalues$plot, device = input$selectExportFileType, width = 15, height = 15, units = 'cm' )


# generate and save 20 graphs using ExplodeLayout with radius from 0.1 to 2.0, step size is 0.1
for (i in 1 : 20) {
  coordX <- columnNames[2 * i + 5]
  coordY <- columnNames[2 * i + 6]

  explodeNet <- cbind.data.frame(returnExplodeNetwork[,1], returnExplodeNetwork[,colnames(returnExplodeNetwork) == coordX],
                                 returnExplodeNetwork[,colnames(returnExplodeNetwork) == coordY], returnExplodeNetwork[,4],
                                 returnExplodeNetwork[,5], returnExplodeNetwork[,6])

  colnames(explodeNet) <- c(columnNames[1], coordX, coordY, columnNames[4], columnNames[5], columnNames[6] )

  returnList <- list("explode.net" = explodeNet)#, "edgeList" = edgeList)

  newSliderValues <- returnList

  explode.net <- newSliderValues$explode.net
  #ranges$x <- c(min(explode.net[,2]), max(explode.net[,2]))
  #ranges$y <- c(min(explode.net[,3]), max(explode.net[,3]))

  p <-  suppressWarnings(plotnetwork(input, newSliderValues$explode.net, rvalues, 1, labelsData))


  #plotcord <- as.data.frame(matrix(as.numeric(sliderValues[,2:3]), nrow(sliderValues), 2, byrow=F))
  plotcord <- as.data.frame(newSliderValues$explode.net[,2:3])
  colnames(plotcord) = c("X1", "X2")


  p <- p + ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y)

  if(labelOn == 'yes' && ((input$labelSize != 0) && (!is.na(input$labelSize)))){
    p <- p + suppressWarnings( #rvalues$colnames
      ggrepel::geom_label_repel(ggplot2::aes(plotcord$X1, plotcord$X2,label=lab, alpha=0.5),
                       na.rm = TRUE,
                       size = input$labelSize,
                       fontface = "bold" ,
                       box.padding = ggplot2::unit(0.5, "lines"),
                       point.padding = ggplot2::unit(0.5, "lines"),
                       segment.color = "#0000FF"
                       # remove lollipop stick : segment.size = 0
      )
    )
  }

  rvalues$plot <- p
  exportFile <- paste(data_name,'_radius',as.double(i/10),'.',input$selectExportFileType, sep="")
  ggplot2::ggsave(exportFile, rvalues$plot, device = input$selectExportFileType, width = 15, height = 15, units = 'cm' )
}
}
