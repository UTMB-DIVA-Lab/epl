#' Converting data into format used for the ExplodeLayout(EL) Algorithm
#'
#' This function takes user's data, including the adjacency or incidence matrix representing networks,
#' and vectors containing information of the vertices' cluster membership, entity, etc.
#'
#' The output of the data would be a list containing two dataframes, 'nodelist' and 'net', which would
#' serve as the input of the next 'search' step for EL algorithm.
#'
#'@param net_matrix a numeric matrix or dataframe storing the network as adjacency matrix(for unipartite networks)
#' or incidence matrix(for bipartite networks). It would be great if you have vertices names contained.
#'@param cluster a numeric vector containing the cluster membership for each vertices. Make sure the length
#' matches the data matrix file.
#'@param outcome a numeric vector containing the outcome for each vertice(such as 1 for cases and 2 for controls)
#'@param entity a numeric vector representing the entity for each vertice (if it's unipartite, put all 1, if its'
#'bipartite network, put one type vertice as 1, the other as 2)
#'
#'
#'@export
dataConvert <- function(net_matrix, cluster, outcome, entity) {
  netType <- length(unique(entity))
  graph <- igraph::graph_from_incidence_matrix(net_matrix, dircted = FALSE)
  if (netType == 1) {
    graph <- igraph::graph_from_adjacency_matrix(net_matrix, directed = FALSE)
  }
  KKcoordinates <- igraph::layout_with_kk(graph, dim=2, set.seed(42))
  KKcoordinates <- apply(data.matrix(KKcoordinates), 2, function(x) sapply(x, function(y) (y-min(x))/(max(x)-min(x))))
  FRcoordinates <- igraph::layout_with_fr(graph, dim=2, set.seed(42),grid="nogrid")
  FRcoordinates <- apply(data.matrix(FRcoordinates), 2, function(x) sapply(x, function(y) (y-min(x))/(max(x)-min(x))))
  Label <- row.names(net_matrix)
  if (netType == 2) {
    Label <- c(Lable, colnames(net_matrix))
  }
  epl_data <- list()
  epl_data$net <- cbind(outcome, net_matrix)
  vertices <- data.frame(Label, FRx=FRcoordinates[,1], FRy=FRcoordinates[,2],
                         KKx=KKcoordinates[,1], KKy=KKcoordinates[,2],
                         Outcome = outcome, Cluster = cluster, Entity = entity)
  epl_data$nodelist <- vertices
  return (epl_data)
}
