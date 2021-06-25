#' Find the leaves under a node
#'
#' This function finds all leaves that are descendants of a node in a tree (represented by an \code{hc_list} object).
#' @param ind Index of a tree node.
#' @param hc_list A list of length-\code{num_interior_nodes} where the ith item in the list contains the child nodes of the ith node in the tree. The negative values in the list indicate leaf nodes.
#' @return Returns the leaf indices.
#' @examples
#' ## See data example.
#' @export
find_leaves_in_list = function(ind, hc_list){
  # This function recursively finds the leaves under a node
  if(ind<0){
    # leaf
    v = c()
    v[1] = -ind
    return(v)
  }else{
    v = c()
    for(i in hc_list[[ind]]){
      parts = find_leaves_in_list(i, hc_list)
      v = unique(c(v, parts))
    }
  }
  return(v)
}
