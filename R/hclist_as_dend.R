#' Transform a hc_list object into a dendrogram
#'
#' This function transforms a dendrogram into a hc_list object
#' @param hc_list An object of class \code{hc_list}.
#' @param leaf_labels A length-\code{n-leaf} of labels of leaves. If null, the leaves will be coded from 1 to \code{n-leaf}.
#' @return An object of class \code{dendrogram} that stores a tree structure.
#' \item{dend}{An object of class \code{dendrogram}.}
#' @examples
#' hc <- hclust(dist(USArrests), "ave")
#' dend = as.dendrogram(hc)
#' hclist_obj = dend_as_hclist(dend)
#' dend2 = hclist_as_dend(hc_list = hclist_obj$hc_list, hclist_obj$leaf_labels)
#' plot(dend2, center = TRUE)
#' # can compare with plot(dend)
#' @export
hclist_as_dend = function(hc_list, leaf_labels = NULL){
  p = abs(min(unlist(hc_list)))
  if(is.null(leaf_labels)){
    leaf_labels = -(1:p)
  }
  depths = find_depths_in_list(hc_list)
  max_depth = max(depths)
  ys = rev(seq(5, 150, length.out = max(depths)))[depths] + runif(length(depths), -5, 5)
  go_through_list = function(ith, hc_list, ys, depths){
    branch_ith = list()
    if(ith<0){
      attributes(branch_ith)<-list(members = 1, label = leaf_labels[-ith],
                                    leaf = TRUE, height = ys[-ith])


      class(branch_ith) <-"dendrogram"
      return(branch_ith)
    }else{
      attributes(branch_ith)<- list(members = length(find_leaves_in_list(ith, hc_list)),
                                    label = ith,
                                    height = ys[ith+p],
                                    leaf = FALSE)
      for(j in 1:length(hc_list[[ith]])){
        branch_ith[[j]] = go_through_list(hc_list[[ith]][j], hc_list, ys, depths)
      }
      class(branch_ith) <-"dendrogram"
      return(branch_ith)
    }
  }
  return(go_through_list(length(hc_list), hc_list, ys = ys, depths))
}

