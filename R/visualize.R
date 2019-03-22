#'The 'visualize' step of the ExplodeLayout(EL) algorithm.
#'
#'This function takes as input the coordinates file and edgelist file from previous step and visuallize
#' the networks generated from EL algorithm in saved local graph files.
#'
#'Input:
#'The input of the visualize function includes mainly two data sets, which should be the output of the
#'previous 'search' step: 1. A 'coordinates' data containing the coordinates of all vertices in networks
#'with different radius(a parameter in EL algorithm). The radius range is from 0.1 to 2.0 with step of 0.1.
#'2. An 'edgelist' file which containing all the edges in the network and their weight.
#'
#'Please note: This 'visualize' function is designed to visualize networks layouts generated using the
#'ExplodeLayout algorithm. The input data sets should always be the output of the previous 'search' step.
#'
#'Output:
#'This function will generate network of particular radius(between 0.1 and 2.0) or all of them
#'at once depending on user's need. Plus, the initial layout of the network would always be saved
#'as a comparison to the 'exploded' networks. User may want to refer to the CCS scores saved in
#'the 'statList' file generated previously to determine the radius they'd like to set.
#'
#'@param coordinatesFile a dataframe containing coordinates of all vertices for networks with radius
#'ranging from 0.1 to 2.0
#'@param edgelistFile a datafram containing three columns representing all edges in the network and
#'their weight
#'@param projName a string which is the project name or data name specified by user. This name will be used in saved export files.
#'@param radius a numeric that represents the radius of the network user would like to visualize. It should be
#'a value in 0.1, 0.2 ... 2.0. If the user wants to output all the possible graphs, set radius equal to 0.
#'@param nodeSize a numeric represents the size of vertices in the network visualization
#'@param entityDegree a boolean variable specifying whether the vertice size should reflect its degree(vertice with more edges would be larger
#'and vice versa)
#'@param labelSize a numeric represents the size of vertices label. If user prefers not showing the label, set it as 0.
#'@param edgeThickness a numeric specifying the thickness of edges in the output graphs
#'@param edgeWeight a boolean variable specifying whether to show the weight of edges. If disabled, all edges in the output graph have the same thickness.
#'@param clusterColor a boolean variable specifying whether to show clusters using different colors or just give the same color to the whole network
#'@param picFormat a string specifying the graph file format, currently we support 'png', 'svg', and 'jpeg'
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
