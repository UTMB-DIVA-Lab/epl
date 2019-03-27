#'\pkg{epl} ExplodeLayout Algorithm.
#'
#'This package provides functions for the ExplodeLayout (EL) algorithm, described in the following
#'publication: \emph{Bhavnani S.K., Chen, T., Ayyaswamy, A., Visweswaran, S., Bellala, G.} \strong{Enabling Comprehension
#' of Patient Subgroups and Characteristics in Large Bipartite Networks: Implications for Precision Medicine.}
#' \emph{Proceedings of AMIA Summit on Translational Bioinformatics (2017).}
#'
#'Motivation: Despite strong and significant clustering, many networks look like "hairballs" when visualized
#'using standard layout algorithms such as Kamada Kawai and Fruchterman Rheingold. The ExplodeLayout algorithm
#'addresses this problem by separating given clusters using the analogy of an explode layout drawing utilized
#'to understand complex assemblies in mechanical engineering.
#'
#'Overview of the EL Algorithm: The EL algorithm takes as input (1) the network data (unipartite or bipartite),
#'(2) node cluster membership (generated from a cluster algorithm such as modularity, or hierarchical clustering),
#'and (3) node layout coordinates (generated from a layout algorithm such as Kamada-Kawai and Fruchterman-Reingold).
#'The algorithm uses the above input to calculate the centroid of each cluster, and moves all nodes in a cluster
#'radially outward such that the centroid is incident to a given circle, and the distances of the nodes to their
#'respective centroid are preserved. Because this movement alters the rotation of the clusters relative to the original
#'layout, the clusters are rotated to match their original orientation. The algorithm generates the radius of the
#'circle through a search, which optimizes the ratio of the overlap of the bounding boxes defining each cluster, to
#'the bounding box defining the entire network. The algorithm uses this optimal explode layout radius to generate a
#'range of layouts in its close vicinity, and displays the network by connecting the nodes based on the inputted network.
#'The EL algorithm has a computational complexity of n.
#'
#'The algorithm provides parameters that enable the user to (1) input data without layout coordinates and to instruct
#'EL to generate those coordinates using a standard layout algorithm, (2) modify default methods for generating the
#'centroid, exploded coordinates, and range of circle radiuses, and (3) modify how the network is visualized.
#'
#'#The EL package provides three functions: (1) 'dataConvert' converts user's network data into a format required by EL,
#'(2) 'search' generates a range of exploded network layouts, and (3) 'visualize' generates visualizations of the
#'exploded network layouts.
#'
#'A sample dataset in the EL input format can be obtained from {https://github.com/UTMB-DIVA-Lab/epl}) Full documentation
#' can be obtained from: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5543384/}
#'
#'@docType package
#'@name epl
NULL
