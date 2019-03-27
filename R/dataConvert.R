#'Function to convert data into the format required by the ExplodeLayout (EL) algorithm
#'
#'This function takes (1) a matrix representing networks, (2) cluster membership, (3) outcome  and (4)
#'entity of each node (i.e. type of each node. In a bipartite network, 1 for the first type and 2 for the
#'other. In a unipartite network, all nodes has entity value of 1.)
#'
#'Input:
#'
#'(1) The input network matrix is an adjacency matrix for a unipartite network or biadjacency for a
#'bipartite network. (see definition of adjacency matrix from wikipedia:\url{https://en.wikipedia.org/wiki/Adjacency_matrix}.)
#'
#'(2) Cluster membership is an integer vector which has the length equal to the number of nodes in the
#'network. Each element in this vector represents the cluster that the corresponding node belongs to.
#'
#'(3) The outcome is an integer vector with the same length as cluster membership vector. For example, in
#'a patient-variable data sets, some patients are selected as cases(unhealthy) and some are
#'controls(healthy), then the outcome could be set as 1 for case and 2 for control. If outcome information
#'is not provided, it would be set to all 1s by default.
#'
#'(4) The entity is an integer vector with the same length as cluster membership and outcome. In a
#'unipartite network, entity for each node will set as 1. In a bipartite network, entity is 1 for one type of
#'nodes and 2 for the other.
#'
#'Output:
#'
#'The output of the function would be a list containing two dataframes, 'nodelist' and 'net', which will
#'serve as the input for the next 'search' step for EL algorithm. The 'nodelist' contains the following
#'information:label names for all nodes (will be set to default as number if not specified in the input
#'data); initial coordinates of all nodes using the Fruchterman-Reingold and Kamada-Kawai layout
#'algorithms; outcome, cluster membership and entity for each node. The 'net' dataframe
#'stores the network.
#'
#'
#'@param net_matrix A numeric matrix or dataframe storing the network as (bi)adjacency matrix, row and column names
#'within the data is recommended
#'@param cluster A numeric vector containing the cluster membership for each nodes. Make sure the length
#' matches the data matrix file.
#'@param outcome A numeric vector containing the outcome for each vertice, length should be the same
#'as 'cluster'
#'@param entity A numeric vector representing the entity for each vertice , length should be the same as 'cluster'
#'and 'outcome'.
#'
#'
#'@export
dataConvert <- function(net_matrix, cluster, outcome, entity) {
  netType <- length(unique(entity))
  if (netType > 2) {stop("Data not supported. Make sure your network is either unipartite or bipartite.")}

  if (netType == 1) {
    x <- dim(net_matrix)[1]
    y <- dim(net_matrix)[2]
    if (x != y) {stop("Data not supported. Matrix of unipartite network must have same number of rows and columns!")}
    if (x != length(cluster) || x != length(outcome) || x != length(entity)) {
      stop("Data not supported. Length of cluster, outcome or entity does not match number of nodes in the network!" )
    }
    mode <- isDirected(net_matrix)
    graph <- igraph::graph_from_adjacency_matrix(as.matrix(net_matrix), mode = mode, weighted = T)
  } else {
    graph <- igraph::graph_from_incidence_matrix(as.matrix(net_matrix),directed = F, weighted = T)
  }
  KKcoordinates <- igraph::layout_with_kk(graph, dim=2, set.seed(42))
  KKcoordinates <- apply(data.matrix(KKcoordinates), 2, function(x) sapply(x, function(y) (y-min(x))/(max(x)-min(x))))
  FRcoordinates <- igraph::layout_with_fr(graph, dim=2, set.seed(42),grid="nogrid")
  FRcoordinates <- apply(data.matrix(FRcoordinates), 2, function(x) sapply(x, function(y) (y-min(x))/(max(x)-min(x))))
  Label <- row.names(net_matrix)
  if (netType == 1) colnames(net_matrix) <- Label
  if (netType == 2) {
    Label <- c(Lable, colnames(net_matrix))
  }
  epl_data <- list()
  epl_data$net <- cbind(outcome, net_matrix)
  nodes <- data.frame(Label, FRx=FRcoordinates[,1], FRy=FRcoordinates[,2],
                         KKx=KKcoordinates[,1], KKy=KKcoordinates[,2],
                         Outcome = outcome, Cluster = cluster, Entity = entity)
  epl_data$nodelist <- nodes
  return (epl_data)
}

isDirected <- function(net_matrix) {
  net <- as.matrix(net_matrix)
  n <- dim(net)[1]
  for (i in 1 : n) {
    for (j in 1 :n) {
      if (net[i,j] != net[j, i]) return("undirected")
    }
  }
  return ("directed")
}
