data_from_matrix <- function(net_matrix, cluster, outcome, entity) {
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

data_from_edgelist <- function(net_edgelist, cluster, outcome, entity) {
  net <- 
  return (data_convert_matrix())
}