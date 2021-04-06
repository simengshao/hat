#' @param tree An hclust (if binary tree), or dendrogram, or hc_list object that
#'   stores the tree structure. An hc_list object is a list of
#'   length-\code{num_interior_nodes} where the ith item in the list contains
#'   the child nodes of the ith node in the tree. The negative values in the
#'   list indicate leaf nodes.
