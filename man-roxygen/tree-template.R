#' @param tree An object encoding the tree structure.  Can be one of three
#'   formats: (1) an `hclust` object (if tree is binary), (2) a `dendrogram`, or
#'   (3) a generalization of an hclust object to the case of non-binary trees,
#'   which we call an hc_list object. An hc_list object is a list of
#'   length-\code{num_interior_nodes} where the ith item in the list contains
#'   the child nodes of the ith node in the tree. The negative values in the
#'   list indicate leaf nodes.
