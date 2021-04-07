#' Visualize the Tree Structure
#'
#' This function plots a tree structure. The tree structure can be represented by a variety of object formats: \code{dendrogram}, \code{hclust} or a list (\code{hc_list}).
#' @template tree-template
#' @param iflabel A boolean variable indicating whether leaves should be labeled.
#' @param labels A vector of length-\code{n-leaf} of labels corresponding to each leaf.
#' @examples
#' ## Example 1: Plot an hc object
#' hc = hclust(dist(USArrests), "ave")
#' labels = rownames(USArrests)[hc$order]
#' plot_tree(tree = hc, iflabel = TRUE)
#' ## Example 2: Plot a dendrogram object
#' dend = as.dendrogram(hc)
#' plot_tree(tree = dend, iflabel = TRUE)
#' ## Example 3: Plot a hc_list object
#' dah = dend_as_hclist(dend)
#' labels = attr(dah, "leaf_labels")
#' plot_tree(tree = dah, iflabel = TRUE, labels = labels)
#' @importFrom graphics par plot segments text
#' @importFrom stats is.leaf runif as.dendrogram
#' @export
plot_tree = function(tree = NULL, iflabel = FALSE, labels = NULL){

  # transform dendrogram to a list of length |T\L| (number of interior nodes). Item i in the list stores the children node of node i.
  if(class(tree) == "dendrogram"){
        dah = dend_as_hclist(tree)
        hc_list = dah
        labels = attr(dah, "leaf_labels")
    }else if(class(tree) == "hclust"){
        dend = as.dendrogram(tree)
        dah = dend_as_hclist(dend)
        hc_list = dah
        labels = attr(dah, "leaf_labels")
    }else if(class(tree) == "list"){
      hc_list = tree
    }else{
      stop("No proper tree structure is provided.")
  }
  nolabel = all(is.na(labels))
  if(iflabel == TRUE){
    if(nolabel == TRUE)
      stop("No labels are provided.")
  }else{
    nolabel = TRUE
  }

  p = -min(unlist(hc_list))
  n_interiornodes = length(hc_list)
  depths = find_depths_in_list(hc_list)
  num_each_depth = as.numeric(table(depths))
  leaves_in_order = find_leaves_in_list(length(hc_list),hc_list)


  orders = match(1:p,leaves_in_order)
  nodes_x = rep(0, n_interiornodes)

  for(i in 1:n_interiornodes){
    children = hc_list[[i]]
    nodes_x[i] = mean(sapply(children, function(x){ifelse(x>0, nodes_x[x], orders[-x])}))
  }


  nodes_y = rev(seq(5, 150, length.out = max(depths)))[depths[-(1:p)]]
  nodes_y = nodes_y + runif(length(nodes_y), -5, 5)
  leaves_y = rev(seq(5, 150, length.out = max(depths)))[depths[1:p]]
  leaves_y = leaves_y + runif(length(leaves_y), -5, 5)

  par(xpd = FALSE, mar = c(0,0,0,0))
  plot(NA,  xaxt = 'n', yaxt = 'n',bty = 'n', pch = '', ylab = '', xlab = '', xlim=c(0, p+1), ylim=c(0,155))

  for(i in 1:(n_interiornodes-1)){
    segments(nodes_x[i], nodes_y[i], nodes_x[i], nodes_y[find_parent_in_list(i+p, hc_list)-p])
  }

  for(i in 1:n_interiornodes){
    for(h in hc_list[[i]]){
      if(h>0){
        segments(nodes_x[h], nodes_y[i], nodes_x[i], nodes_y[i])
      }else{
        segments(orders[-h], nodes_y[i], nodes_x[i], nodes_y[i])
      }
    }
  }

  for(i in 1:p){
    parent_i = find_parent_in_list(i, hc_list)
    segments(orders[i], leaves_y[i], orders[i], nodes_y[parent_i-p])
    text(orders[i], leaves_y[i], labels[i], pos = 1, srt = 90, offset = 0.11*nchar(labels[i]), cex = 0.5)

  }

}








#' Visualize Aggregation on the Tree
#'
#' This function plots a tree structure with branches and leaves colored by group memberships. The tree structure can be represented by a variety of object formats: \code{dendrogram}, \code{hclust} or a a list (\code{hc_list}). The aggregation is determined by a vector corresponding to each interior nodes
#' @param rejections A vector of length-\code{num_interior_nodes} indicating if each node is rejected. When coloring with true aggregation, this vector indicate whether each interior node is a non-null node.
#' @param iflabel A boolean variable indicating whether leaves should be labeled.
#' @param labels A vector of length-\code{n-leaf} of labels corresponding to each leaf.
#' @template tree-template
#' @examples
#' ## Example 1: Plot an hc object
#' hc = hclust(dist(USArrests), "ave")
#' labels = rownames(USArrests)[hc$order]
#' plot_aggregation(rejections = c(rep(FALSE, 45), rep(TRUE, 4)), iflabel = TRUE, tree = hc)
#'
#'
#' ## Example 2: Plot a dendrogram object
#' hc = hclust(dist(USArrests), "ave")
#' dend = as.dendrogram(hc)
#' plot_aggregation(rejections = c(rep(FALSE, 45), rep(TRUE, 4)), iflabel = TRUE, tree = dend)
#'
#' ## Example 3: Plot a hc_list object
#' hc = hclust(dist(USArrests), "ave")
#' dend = as.dendrogram(hc)
#' dah = dend_as_hclist(dend)
#' labels = attr(dah, "leaf_labels")
#' plot_aggregation(rejections = c(rep(FALSE, 45), rep(TRUE, 4)), iflabel = TRUE, labels = labels, tree = dah)
#' @importFrom graphics par plot segments text
#' @importFrom stats is.leaf runif
#' @export
plot_aggregation = function(rejections, iflabel = FALSE, labels = NULL, tree=NULL){


  # transform dendrogram to a list of length |T\L| (number of interior nodes). Item i in the list stores the children node of node i.
  if(class(tree) == "dendrogram"){
    dah = dend_as_hclist(tree)
    hc_list = dah
    labels = attr(dah, "leaf_labels")
  }else if(class(tree) == "hclust"){
    dend = as.dendrogram(tree)
    dah = dend_as_hclist(dend)
    hc_list = dah
    labels = attr(dah, "leaf_labels")
  }else if(class(tree) == "list"){
    hc_list = tree
  }else{
    stop("No proper tree structure is provided.")
  }

  nolabel = all(is.na(labels))
  if(iflabel == TRUE){
    if(nolabel == TRUE)
      stop("No labels are provided.")
  }else{
    nolabel = TRUE
  }



  p = -min(unlist(hc_list))
  n_interiornodes = length(hc_list)
  depths = find_depths_in_list(hc_list)
  num_each_depth = as.numeric(table(depths))

  leaves_in_order = find_leaves_in_list(length(hc_list),hc_list)
  orders = match(1:p,leaves_in_order)
  nodes_x = rep(0, n_interiornodes)

  for(i in 1:n_interiornodes){
    children = hc_list[[i]]
    nodes_x[i] = mean(sapply(children, function(x){ifelse(x>0, nodes_x[x], orders[-x])}))
  }

  if(!nolabel){
    nodes_y = rev(seq(15, 150, length.out = max(depths)))[depths[-(1:p)]]
    leaves_y = rev(seq(15, 150, length.out = max(depths)))[depths[1:p]]
  }else{
    nodes_y = rev(seq(5, 150, length.out = max(depths)))[depths[-(1:p)]]
    leaves_y = rev(seq(5, 150, length.out = max(depths)))[depths[1:p]]
  }
  nodes_y = nodes_y + runif(length(nodes_y), -5, 5)
  leaves_y = leaves_y + runif(length(leaves_y), -5, 5)

  par(xpd = FALSE, mar = c(0,0,0,0))
  plot(NA,  xaxt = 'n', yaxt = 'n',bty = 'n', pch = '', ylab = '', xlab = '', xlim=c(0, p+1), ylim=c(0,155))


  for(i in 1:(n_interiornodes-1)){
    segments(nodes_x[i], nodes_y[i], nodes_x[i], nodes_y[find_parent_in_list(i+p, hc_list)-p], lwd = 0.8)
  }

  for(i in 1:n_interiornodes){
    for(h in hc_list[[i]]){
      if(h>0){
        segments(nodes_x[h], nodes_y[i], nodes_x[i], nodes_y[i], lwd = 0.8)
      }else{
        segments(orders[-h], nodes_y[i], nodes_x[i], nodes_y[i], lwd = 0.8)
      }
    }
  }


  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package \"RColorBrewer\" needed for this function to work. Please install it.",
         call. = FALSE)
  }



  # transform rejection to clusters
  rejected_nodes = which(rejections==TRUE)
  all_end_nodes = sort(unique(unlist(sapply(rejected_nodes, function(x){
    hc_list[[x]]
  }))))
  end_nodes = all_end_nodes[!all_end_nodes %in% rejected_nodes]

  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  mycolors = sample(sort(col_vector[-c(53, 22, 40, 28)], decreasing = TRUE), length(end_nodes))

  count = 1
  for(i in end_nodes){
    if(i>0){
      positive_desc = find_descendants_interior_in_list(hc_list, i)
      for(j in positive_desc){
        segments(nodes_x[j], nodes_y[j], nodes_x[j], nodes_y[find_parent_in_list(j+p, hc_list)-p], col = mycolors[count])
        for(h in hc_list[[j]]){
          if(h>0){
            segments(nodes_x[h], nodes_y[j], nodes_x[j], nodes_y[j], col =  mycolors[count])
          }else{
            segments(orders[-h], nodes_y[j], nodes_x[j], nodes_y[j], col =  mycolors[count])
            segments(orders[-h], leaves_y[-h], orders[-h], nodes_y[j], col = mycolors[count])
            if(!nolabel){
              text(orders[-h], leaves_y[-h], labels[-h], pos = 1, srt = 90, offset = 0.11*nchar(labels[-h]),col = mycolors[count], cex = 0.5)
            }
          }
        }
      }
    }else{
      parent_i = find_parent_in_list(-i, hc_list)
      segments(orders[-i], leaves_y[-i], orders[-i], nodes_y[parent_i-p], col = mycolors[count])
      if(!nolabel){
        text(orders[-i], leaves_y[-i], labels[-i], pos = 1, srt = 90, offset = 0.11*nchar(labels[-i]), col = mycolors[count], cex = 0.5)
      }
    }
    count = count+1
  }

}


