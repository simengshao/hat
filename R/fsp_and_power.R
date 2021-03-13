#' Calculate False Split Proportion (FSP)
#'
#' This function calculate False Split Proportion (FSP) of an achieved aggregation given the true aggregation.
#' @param theta_est A vector of integers indicating the cluster to which each item is allocated by an achieved aggregation.
#' @param theta_true A vector of integers indicating the cluster to which each item is allocated by the true aggregation.
#' @return False Split Proportion.
#' @examples
#' v1 = c(1,1,2,2,2,3,3,3,3,4,4,4)
#' v2 = c(1,1,1,1,1,2,2,2,3,4,4,4)
#' calculate_fsp(v1, v2)
#' ## 0.3333333
#' @export
calculate_fsp = function(theta_est, theta_true){
  if(length(theta_est)!=length(theta_true)){
    stop("The two vectors need to be of the same length.")
  }
  true_num_groups = length(unique(theta_true))
  R = 0  # number of positive
  V = 0  # number of false positive
  for (i in 1:true_num_groups){
    V = V + length(unique(theta_est[which(theta_true == unique(theta_true)[i])])) - 1
  }
  R = length(unique(theta_est)) - 1
  return(ifelse(R>0, V/R, 0))
}

#' Calculate Power
#'
#' This function calculate power of an achieved aggregation given the true aggregation.
#' @param theta_est A vector of integers indicating the cluster to which each item is allocated by an achieved aggregation.
#' @param theta_true A vector of integers indicating the cluster to which each item is allocated by the true aggregation.
#' @return Power.
#' @examples
#' v1 = c(1,1,2,2,2,3,3,3,3,4,4,4)
#' v2 = c(1,1,1,1,1,2,2,2,3,4,4,4)
#' calculate_power(v1, v2)
#' ## 0.6666667
#' @export
calculate_power = function(theta_est, theta_true){
  if(length(theta_est)!=length(theta_true)){
    stop("The two vectors need to be of the same length.")
  }
  true_num_groups = length(unique(theta_true))
  ALT = true_num_groups - 1 # number of alternatives

  S = 0 # number of rejected alternatives
  est_num_groups = length(unique(theta_est))
  for (i in 1:est_num_groups){
    S = S + length(unique(theta_true[which(theta_est == unique(theta_est)[i])])) - 1
  }
  return(ifelse(ALT>0, 1- S/ALT, 1))
}
