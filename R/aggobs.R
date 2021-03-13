pvalue_anova_all = function(y, hc_list, sigma_constant = NULL, simes = TRUE){
  # y: vector of observed value
  # hc_list: list length of $T\L$ for storing tree structure
  # sigma_constant: standard deviation of noise

  n_interiornodes = length(hc_list)
  ps = rep(0, n_interiornodes)

  for(i in length(ps):1){
    n_deg = length(hc_list[[i]])
    n_leaves = length(find_leaves_in_list(i, hc_list))
    if(n_deg == 1){
      ps[i] = 1
      next
    }
    y_i_bar = rep(0, n_deg)
    n_i = rep(0, n_deg)
    sum_squares_i = rep(0, n_deg)
    for(j in 1:n_deg){
      leaves_child_i = y[find_leaves_in_list(hc_list[[i]][j], hc_list)]
      y_i_bar[j] = mean(leaves_child_i)
      n_i[j] = length(leaves_child_i)
      sum_squares_i[j] = sum((leaves_child_i - y_i_bar[j])^2)
    }
    y_bar = sum(y_i_bar*n_i)/sum(n_i)

    sigma_est = sd(y)
    if(is.null(sigma_constant)){

      # calculate the F-statistic
      sb_2 = sum(n_i * (y_i_bar - y_bar)^2)/(n_deg-1)
      sw_2 = sum(sum_squares_i)/(n_leaves - n_deg)

      if(n_leaves == n_deg){
        p_value = pchisq(sb_2/sigma_est^2, df= n_deg-1 , lower.tail=FALSE)
      }else{
        p_value = stats::pf(sb_2/sw_2, n_deg-1, n_leaves - n_deg); if(p_value>.5) p_value = 1-p_value
      }
    }else{
      # calculate the chi-statistic and p-value
      sb_2 = sum(n_i * (y_i_bar - y_bar)^2)/(n_deg-1)
      p_value = stats::pchisq(sb_2/sigma_constant^2, df= n_deg-1 , lower.tail=FALSE)
    }
    ps[i] = p_value
  }

  # Simes' procedure
  if(simes){
    for(i in n_interiornodes:1){
      descendants_i = find_descendants_interior_in_list(hc_list, i)
      ps_descendants = sort(ps[descendants_i])
      ps[i] = min(ps_descendants * length(ps_descendants) / seq(1, length(ps_descendants), by = 1))
    }
  }
  return(ps = ps)
}



#' Aggregate observations hierarchically with False Split Rate control
#'
#' This function aggregate observations with the same means while simultaneously controlling False Split Rate under a target level.
#' The aggregation is achieved by two steps: (1) Generate p-values for each interior node through ANOVA test (2) Sequentially test on the tree.
#' @param y A length-\code{n-observation} vector of observations
#' @param sigma_constant Standard deviation of noise.
#' @param hc An object of class \code{hclust}.
#' @param dend An object of class \code{dendrogram}.
#' @param hc_list A list of length-\code{num_interior_nodes} where the ith item in the list contains the child nodes of the ith node in the tree. The negative values in the list indicate leaf nodes. When \code{hc_list} is \code{NULL}, function will learn it from \code{hc} or \code{dend}.
#' @param alpha_level A use-specified target FSR level
#' @return Returns the aggregation result.
#' \item{alpha}{The target FSR level.}
#' \item{groups}{A length-\code{n-observation} vector of integers indicating the cluster to which each observation is allocated.}
#' \item{rejections}{A length-(\code{num_interior_nodes}) vector indicating whether each node is rejected.}
#' @examples
#' ## See data example.
#' @importFrom stats pchisq pf
#' @export
aggregate_observations = function(y, sigma_constant = NULL, hc= NULL, dend = NULL, hc_list = NULL, alpha_level){
  if(is.null(hc_list)){
  # transform dendrogram to a list of length |T\L| (number of interior nodes). Item i in the list stores the children node of node i.
    if(is.null(hc)){
      if(class(dend) != "dendrogram"){
        stop("The tree needs to be an dendrogram object.")
      }else{
        hc_list = dend_as_hclist(dend)$hc_list
      }
    }else{
     if(class(hc) != "hclust"){
        stop("The tree needs to be an hclust object.")
      }else{
        dend = stats::as.dendrogram(hc)
        hc_list = dend_as_hclist(dend)$hc_list
      }
    }
  }

  p_vals = pvalue_anova_all(y, hc_list, sigma_constant, simes = TRUE)

  result = hierarchical_test(hc_list = hc_list, p_vals = p_vals, alpha_level = alpha_level, independent = FALSE)
  groups = result$groups
  rejections = result$rejections
  return(list(alpha = alpha_level, groups = groups, rejections = rejections))

}
