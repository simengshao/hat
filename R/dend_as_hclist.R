#' Transform a dendrogram into a list object
#'
#' This function transforms a dendrogram into a hc_list object
#' @param dend An object of class \code{dendrogram}.
#' @return A list object of length-\code{num_interior_nodes}. The \code{i}-th item in the list contains the child nodes of the \code{i}-th node in the tree. The negative values in the list indicate leaf nodes, and positive values indicate interior node. The interior nodes are numbered reversely from the root along the branch, i.e., the number corresponding to a node is always smaller than the number corresponding to its parent.
#' \item{hc_list}{The list object that represents the tree structure,}
#' @examples
#' hc <- hclust(dist(USArrests), "ave")
#' dend = as.dendrogram(hc)
#' dah = dend_as_hclist(dend)
#' ## check length of the list
#' length(dah)
#' ## leaf labels
#' attr(dah, "leaf_labels")
#' @importFrom stats dendrapply as.dendrogram na.omit
#' @export
dend_as_hclist = function(dend){
  p = dendextend::nleaves(dend)
  m = dendextend::nnodes(dend)
  leaf_labels = dendextend::get_leaves_attr(dend, "label")
  not_labeled = is.null(attr(dend, "label"))
  if(not_labeled)
  {
    dend = dendextend::assign_values_to_nodes_nodePar(dend, seq(1:m))
  }

  i_node_counter <- 0
  hc_list = vector("list", m)
  leaf_nodes = c()
  go_through_nodes <- function(dend_node, hc_list) {
    i_node_counter <<- i_node_counter + 1
    if(!is.leaf(dend_node)){
      n_children = length(dend_node)
      tmp = c()
      for(i in 1:n_children){
        if(not_labeled){
          tmp = c(tmp, attr(dend_node[[i]], "nodePar")["pch"])
        }else{
          tmp = c(tmp, attr(dend_node[[i]], "label"))
        }

      }
      hc_list[[i_node_counter]] <<- tmp
    }else{
      leaf_nodes <<- c(leaf_nodes, i_node_counter)
    }
  }

  dendrapply(dend, go_through_nodes, hc_list)


  # Remove leaf nodes
  hc_list = Filter(function(a) any(!is.na(a)), hc_list)

  if(not_labeled){
    all_nodes = 1:m
    # Recode nodes with positive 1:(m-p) corresponding to interior nodes
    # and negative (1:p) corresponding to leaf nodes
    non_leaf_nodes = all_nodes[!all_nodes %in% leaf_nodes]
    non_leaf_recoding = cbind(org = non_leaf_nodes, new = length(non_leaf_nodes):1)
    leaf_recoding = cbind(org = leaf_nodes, new = - (1:length(leaf_nodes)))
    recoding = data.frame(rbind(non_leaf_recoding, leaf_recoding))
    hc_list = rev(lapply(hc_list, function(x) {recoding$new[match(x, recoding$org)]}))
  }else{
    node_in_order = as.numeric(na.omit(dendextend::get_nodes_attr(dend, "label", include_leaves = FALSE)))
    hc_list_new = hc_list
    for(i in length(hc_list):1){
      hc_list_new[[node_in_order[i]]] = hc_list[[i]]
    }
    hc_list = hc_list_new
  }

  attr(hc_list, "leaf_labels") = leaf_labels
  return(hc_list)
}
