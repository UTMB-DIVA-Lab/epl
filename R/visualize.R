#'The 'visualize' step of the ExplodeLayout(EL) algorithm.
#'
#''This function takes as input the (1) coordinate file, and (2) edgelist file for each of the exploded
#'networks generated from the 'search' function. The visualize function saves the exploded networks in
#'local folders named by the value of the explode radius (e.g., 2).
#'
#'Input:
#'
#'(1) A 'coordinates' data object containing the coordinates of all nodes in network with different
#'radiuses. The default radius range is 0.1 to 2.0, with increments of 0.1.
#'
#'(2) An 'edgelist' data object containing all edges connecting pairs of nodes in the network,
#'and their weight (if the original network is unweighted, the weight of all edges should be 1).
#'
#'Please note: This 'visualize' function is designed to visualize networks layouts generated using the
#'ExplodeLayout algorithm. Please use the output of the 'search' step as input for the 'visualize' function.
#'
#'Output:
#'
#'The 'visualize' function will generate a network of a particular radius (default between 0.1 - 2.0, at 0.1 increments),
#'or all of them based on parameters specified by the user. Users may refer to the CCS scores saved in the
#''statList' file to select a network based on the CCS scores generated from the search algorithm in the
#''search' function.
#'
#'@param coordinatesFile A dataframe containing coordinates of all nodes for networks with radius
#'ranging from 0.1 to 2.0
#'@param edgelistFile A dataframe containing three columns: node-1, node-2, and edgeweight.
#'@param projName A string specifying the project name used to export files. Default projName = 'defaultName'
#'@param radius A numeric specifying the radius (ranging from 0.1-2.0) of the ExplodeLayout circle. For all
#'radius ranging from 0.1-2.0, set radius=0. Default = 3.
#'@param nodeSize An integer specifying the size of nodes. Default = 3.
#'@param entityDegree A boolean variable specifying whether the nodes should be sized by its degree (number
#'of edges connected to that node). Default = FALSE.
#'@param labelSize An integer specifying the size of the node labels. For no labels, set labelSize=0. Default = 0.
#'@param edgeThickness A numeric specifying the thickness of edges in the visualization. For weighted networks,
#'this value will be multiplied by the edge weights provided, to generate the final edge weight. Default = 0.1.
#'@param edgeWeight A boolean variable specifying whether to display the edge weight (using line thickness) in
#'weighted networks. If edgeWeight = FALSE, then all edge thicknesses will be equal. Default = FALSE.
#'@param clusterColor A boolean variable specifying whether to show clusters using different colors or a uniform
#'color. Default clusterColor = TRUE.
#'@param picFormat A string specifying the network layout output with values 'png', 'svg', and 'jpeg', Default = 'png'
#'
#'
#'@references \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5543384/}
#'
#'@export
visualize <- function(coordinatesFile, edgelistFile, proj_name = 'defaultName', radius = 0, nodeSize = 3,
                      entityDegree = FALSE, labelSize = 0, edgeThickness = 0.1,
                      edgeWeight = FALSE, clusterColor = TRUE, picFormat = 'png') {
  #setting default parameters for explodeLayout
  input <- list()
  rvalues <- list()

  returnExplodeNetwork <- coordinatesFile
  #project/data name
  input$data_name <- proj_name
  #vertex size
  input$vertexSize <- nodeSize
  #show node degree or not
  input$nodeDegree <- entityDegree
  #vertex label size
  input$labelSize <- labelSize
  #set default edge thickness or not
  input$edgeThickness <- edgeThickness
  #show edge weight or not
  input$edgeWeight <- edgeWeight
  #choose node color mode (color by clusters, outcome, etc..)
  input$nodeColor <- '2' #cluster color
  if(!clusterColor) {input$nodeColor <- '3'}
  #select the picture format to export
  input$selectExportFileType <- picFormat

  dataFromFile <- list()
  dataFromFile$edgelist <- edgelistFile
  rvalues$dataFromFile <- dataFromFile
  rvalues$netType <- length(unique(returnExplodeNetwork[,6]))
  if (rvalues$netType > 2) {stop("Data not supported. Make sure your input has only one or two entities.")}

  #plot the original graph
  plotOne(returnExplodeNetwork, input, rvalues, 2, paste(proj_name,'_original', sep = ""))

  #plot single or all graphs
  if (radius == 0) {
    for (i in 1 : 20) {
      name <- paste(proj_name, '_radius', as.double(i/10), sep = "")
      plotOne(returnExplodeNetwork, input, rvalues, 2*i+5, name)
    }
  } else {
    rad <- radius * 10
    name <- paste(proj_name, '_radius', radius, sep = "")
    plotOne(returnExplodeNetwork, input, rvalues, 2*rad+5, name)
  }
}


plotOne <- function(returnExplodeNetwork, input, rvalues, radiusPos, tag) {
  columnNames <- colnames(returnExplodeNetwork)
  radiusPos <- as.integer(radiusPos)
  coordX <- columnNames[radiusPos]
  coordY <- columnNames[radiusPos+1]

  explode.net <- cbind.data.frame(returnExplodeNetwork[,1], returnExplodeNetwork[,columnNames == coordX],
                                 returnExplodeNetwork[,columnNames == coordY], returnExplodeNetwork[,4],
                                 returnExplodeNetwork[,5], returnExplodeNetwork[,6])

  colnames(explode.net) <- c(columnNames[1], coordX, coordY, columnNames[4], columnNames[5], columnNames[6] )

  ranges <- list()
  ranges$x <- c(min(explode.net[,2]), max(explode.net[,2]))
  ranges$y <- c(min(explode.net[,3]), max(explode.net[,3]))

  labelsData <- getDataInBoundingBox(explode.net , ranges, rvalues)
  lab <- labelsData$colnames
  if(rvalues$netType == 2) {
    lab <- c(rep(NA, labelsData$nrows), labelsData$colnames)
  }

  p <-  suppressWarnings(plotnetwork(input, explode.net, rvalues, 1, labelsData))


  #plotcord <- as.data.frame(matrix(as.numeric(sliderValues[,2:3]), nrow(sliderValues), 2, byrow=F))
  plotcord <- as.data.frame(explode.net[,2:3])
  colnames(plotcord) = c("X1", "X2")


  p <- p + ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y)

  if((input$labelSize != 0) && (!is.na(input$labelSize))){
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

  exportFile <- paste(tag, '.',input$selectExportFileType, sep="")
  ggplot2::ggsave(exportFile, p, device = input$selectExportFileType, width = 15, height = 15, units = 'cm')

}
