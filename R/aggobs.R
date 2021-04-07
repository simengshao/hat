pvalue_anova_all = function(y, hc_list, sigma = NULL, simes = TRUE){
  # y: vector of observed value
  # hc_list: list length of $T\L$ for storing tree structure
  # sigma: standard deviation of noise

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
    if(is.null(sigma)){

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
      p_value = stats::pchisq(sb_2/sigma^2, df= n_deg-1 , lower.tail=FALSE)
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
#' @param sigma Standard deviation of noise. If given, the algorithm will compute nodewise p-values with chi-squared statistics. If sigma is unkown, the algorithm will compute p-values with F-test statistics.
#' @template tree-template
#' @param alpha A use-specified target FSR level
#' @return Returns the aggregation result.
#' \item{alpha}{The target FSR level.}
#' \item{groups}{A length-\code{n-observation} vector of integers indicating the cluster to which each observation is allocated.}
#' \item{rejections}{A length-(\code{num_interior_nodes}) vector indicating whether each node is rejected.}
#' @examples
#' set.seed(123)
#' hc = hclust(dist((1:20) + runif(20)/20), method = "complete")
#' k = 4 # 4 true groups
#' groups = cutree(hc, k = 4)
#' theta = runif(k, 0, 10)[groups]
#' y = theta + runif(20, 0, 1)
#' aggregate_observations(y, sigma = 1, tree= hc, alpha = 0.1)
#' @importFrom stats pchisq pf
#' @export
aggregate_observations = function(y, sigma = NULL, tree= NULL, alpha){

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

  p_vals = pvalue_anova_all(y, hc_list, sigma, simes = TRUE)

  result = hierarchical_test(tree = hc_list, p_vals = p_vals, alpha = alpha, independent = FALSE)
  groups = result$groups
  rejections = result$rejections
  return(list(alpha = alpha, groups = groups, rejections = rejections))

}
