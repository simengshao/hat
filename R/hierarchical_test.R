
# This function performs the rejection step at a specific depth d.
make_rejection = function(d, nodes, alpha_level, reshaping_func,
                          p, D, Deltas, p_vals,
                          leaves, depths, degrees, degrees_rej,
                          rejections){

  # Inputs:
  # d: current depth
  # nodes: nodes at depth d that are eligible for testing
  # alpha_level: target FSR level
  # reshaping_func: if the thresholds need to be reshaped
  # p: number of leaves in total
  # D: max depth of non-leaf nodes
  # Delta: maximum degrees of nodes on the tree
  # p_vals: vector of p-values
  # leaves: vector of number of leaves in the subtree rooted at each node
  # degrees: degree of each node on T
  # degree_rej: degree of each node on T_rej
  # rejections: vector that encodes current rejection status


  # Outputs:
  # rejected_d: rejections on depth d
  # crit_func_d: computed threshold values on depth d

  if(length(nodes)==0)
    return(0)


  leaves_d = leaves[nodes]
  p_vals_d = p_vals[nodes-p]
  Deltas = Deltas[nodes]

  delta_min = min(degrees[degrees>1])

  # Initialize rejection result vector
  rejected_d = rep(FALSE, length(nodes))

  possible_r = rev(1:((sum(degrees[nodes]) - length(nodes))))

  for(r in possible_r)
  {
    # a) compute cutoffs across multiple r values
    harmonics_u = 1 + harmonic_diff(p -1 - (sum(degrees[nodes]) - length(nodes) -r), sum(rejections* (degrees - degrees_rej)) + r)
    if(reshaping_func == "ID"){
      numerator = alpha_level * leaves_d * (sum(rejections* (degrees - degrees_rej)) - 1+r)
      term = p*(1-1/(Deltas)^2) * harmonics_u
      crit_func_d = 1/Deltas * numerator / (term + numerator)
      crit_func_d[which(is.na(crit_func_d))] = 1
    }else if(reshaping_func == "BY")
    {
      denoms = harmonic_diff(sum(degrees)-1, d*(delta_min-1))
      crit_func_d =  alpha_level *leaves_d *(sum(rejections* (degrees - degrees_rej)) - 1+r)/
        (p*(Deltas-1/Deltas) * D)/denoms
      crit_func_d[which(is.na(crit_func_d))] = 1
    }

    # b) decide r and who gets rejected
    if(r<=sum(as.numeric(p_vals_d <= crit_func_d) ))
    {
      rejected_d = (p_vals_d <= crit_func_d)
      print(paste0("The rejected nodes in level ", d,
                   " are ", nodes[rejected_d]))
      print(paste0("The critical function at nodes in level ", d,
                   " are ", crit_func_d))
      break
    }

  }
  return(list(rejected_d = rejected_d, crit_func_d = crit_func_d))

}


#' Test on a given tree to achieve aggregation of leaves
#'
#' This function sequentially tests on a tree in a top-down manner. The testing result determines aggregation of the leaves and make sure the False Split Rate is controlled under a target level.
#' @param hc An object of class \code{hclust}.
#' @param dend An object of class \code{dendrogram}.
#' @param hc_list A list of length-\code{num_interior_nodes} where the ith item in the list contains the child nodes of the ith node in the tree. The negative values in the list indicate leaf nodes. When \code{hc_list} is \code{NULL}, function will learn it from \code{hc} or \code{dend}.
#' @param p_vals A length-\code{num_interior_nodes} vector of p-values corresponding to interior nodes.
#' @param alpha_level A use-specified target FSR level
#' @param independent Whether the p-values are independent (default = TRUE).
#' @return Returns the testing result and calculated threshold values.
#' \item{alpha}{The target FSR level.}
#' \item{rejections}{A length-\code{num_interior_nodes} vector indicating whether each node is rejected.}
#' \item{threshold_functions}{A length-(\code{nnodes}-\code{nleaves}) vector of computed threshold value at each node. NA if not tested.}
#' \item{groups}{A length-\code{n-feature} vector of integers indicating the cluster to which each leaf is allocated.}
#' @examples
#' set.seed(123)
#' hc = hclust(dist((1:10) + runif(10)/10), method = "complete")
#' p_vals = c(runif(5), rbeta(4, 1, 60))
#' alpha_level = 0.3
#' independent = TRUE
#' hierarchical_test(hc = hc, dend = NULL, p_vals = p_vals, alpha_level = alpha_level, independent = independent)
#' @export
hierarchical_test = function(hc = NULL, dend = NULL, hc_list = NULL, p_vals, alpha_level, independent = TRUE){

  # transform dendrogram to a list of length |T\L| (number of interior nodes). Item i in the list stores the children node of node i.
  if(is.null(hc_list)){
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
        dend = as.dendrogram(hc)
        hc_list = dend_as_hclist(dend)$hc_list
      }
    }
  }


  n_interiornodes = length(hc_list) # number of interior nodes
  p = -min(unlist(hc_list)) # number of leaves
  m = n_interiornodes + p
  all_nodes = 1: m #including root node and interior nodes

  # Initialize rejections (the output of the algorithm.)
  rejections = rep(FALSE, m)
  threshold_functions = c(rep(NA, m), 1)

  # Find depth of each node (root has depth 1).
  depths = find_depths_in_list(hc_list)
  max_depth = max(depths)

  # Find the number of leaves of subtree rooted at each node
  leaves = number_leaves_in_list(hc_list)

  # Find degree_T and initilize degree_rej for each node
  degrees = c(rep(0,p),sapply(hc_list, function(x) length(x)))
  degrees_rej = rep(0,m)

  # upper bound of degree of the tree
  Deltas = find_max_degree_in_list(hc_list)

  # reshaping function
  # if p-values are independent, then use threshold function as it is
  # if p-values are dependent, then use BY-reshaped threshold function
  reshaping_func = ifelse(independent, "ID", "BY")

  for(d in 1:max_depth)
  {
    # Get the index of nodes all depth d
    nodes_depth_d = all_nodes[depths == d]

    if(d == 1)
    {
      rejections[nodes_depth_d] = TRUE
      # We always reject the root node
      print(paste0("Nodes ", nodes_depth_d ," at level ", d, " are rejected."))
      print(paste0("Now moving to depth ", d+1))
    }
    # Discard the nodes one of whose parents has not been rejected.
    if(d>1)
    {

      nodes_depth_d = nodes_depth_d[nodes_depth_d>p]
      if(length(nodes_depth_d)>0){

        # Only test interior nodes whose parent is rejected
        nodes_to_test_depth_d = nodes_depth_d[rejections[sapply(nodes_depth_d, find_parent_in_list, hc_list)]]
        nodes_to_test_depth_d = nodes_to_test_depth_d[p_vals[nodes_to_test_depth_d-p]!=1]

        # Performs the rejection step at depth d.
        if(length(nodes_to_test_depth_d)>0){
          ifrejected = make_rejection(d, nodes_to_test_depth_d, alpha_level, reshaping_func,
                                      p, D=max_depth-1, Deltas, p_vals,
                                      leaves, depths, degrees, degrees_rej,
                                      rejections)

          # determine who gets rejected
          rejected_nodes_depth_d = nodes_to_test_depth_d[ifrejected$rejected_d]
          if(length(rejected_nodes_depth_d)!=0){
            # update degrees_rej
            rejected_nodes_parent_table = table(sapply(rejected_nodes_depth_d, find_parent_in_list, hc_list))
            for(i in 1:length(rejected_nodes_parent_table)){
              degrees_rej[as.numeric(names(rejected_nodes_parent_table)[i])] = rejected_nodes_parent_table[i]
            }

            # update rejections and thresholds
            rejections[rejected_nodes_depth_d] = TRUE
            threshold_functions[nodes_to_test_depth_d] = ifrejected$crit_func_d
          }
        }
        else{
          print(paste0("Stop moving to next level because no nodes at depth ",
                       d, " are rejected."))
          break
        }
        print(paste0("Now moving to depth ", d+1))

      }else{
        print(paste0("Finished testing at depth", d))
        break
      }
    }
  }
  rejections = rejections[-(1:p)]
  groups = determine_aggregation(hc_list, rejections)
  return(list(alpha = alpha_level,
              rejections = rejections,
              threshold_functions = threshold_functions[-(1:p)],
              groups = groups))
}
