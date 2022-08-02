
# This function performs the rejection step at a specific depth d.
make_rejection = function(d, nodes_depth_d, nodes, alpha, reshaping_func,
                          p, D, Delta, p_vals,
                          leaves, depths, degrees,
                          rejections, R_to_d_minus_1){

  # Inputs:
  # d: current depth
  # nodes_depth_d: nodes at depth d that are non-leaf
  # nodes: nodes at depth d that are eligible for testing
  # alpha: target FSR level
  # reshaping_func: if the thresholds need to be reshaped
  # p: number of leaves in total
  # D: max depth of non-leaf nodes
  # Delta: maximum degrees of nodes on the tree
  # p_vals: vector of p-values
  # leaves: vector of number of leaves in the subtree rooted at each node
  # degrees: degree of each node on T
  # rejections: vector that encodes current rejection status
  # R_to_d_minus_1: R^{1:(d-1)} defined in equation 14


  # Outputs:
  # rejected_d: rejections on depth d
  # crit_func_d: computed threshold values on depth d

  if(length(nodes)==0)
    return(0)


  leaves_d = leaves[nodes]
  p_vals_d = p_vals[nodes-p]


  delta_min = min(degrees[degrees>1])

  # Initialize rejection result vector
  rejected_d = rep(FALSE, length(nodes))

  possible_r = rev(1:((sum(degrees[nodes]) - length(nodes))))


  for(r in possible_r)
  {
    # a) compute cutoffs across multiple r values
    harmonics_u = 1 + harmonic_diff(p -1 - (sum(degrees[nodes_depth_d]) - length(nodes_depth_d) -r), R_to_d_minus_1 + r)
    if(reshaping_func == "ID"){
      numerator = alpha * leaves_d * (R_to_d_minus_1 + r)
      term = p*(1-1/(Delta)^2) * harmonics_u
      crit_func_d = 1/Delta * numerator / (term + numerator)
      crit_func_d[which(is.na(crit_func_d))] = 1
    }else if(reshaping_func == "BY")
    {
      denoms = harmonic_diff(sum(degrees[nodes_depth_d])-1, d*(delta_min-1))
      crit_func_d =  alpha *leaves_d *(R_to_d_minus_1 - 1+r)/
        (p*(Delta-1/Delta) * D)/denoms
      crit_func_d[which(is.na(crit_func_d))] = 1
    }

    # b) decide r and who gets rejected
    if(r<=sum((p_vals_d <= crit_func_d)*(degrees[nodes]-1)))
    {
      rejected_d = (p_vals_d <= crit_func_d)
      R_to_d_minus_1 = R_to_d_minus_1 + r
      print(paste0("The rejected nodes in level ", d,
                   " are ", nodes[rejected_d]))
      print(paste0("The critical function at nodes in level ", d,
                   " are ", crit_func_d))
      break
    }

  }
  return(list(rejected_d = rejected_d, crit_func_d = crit_func_d, R_to_d_minus_1 = R_to_d_minus_1))

}


#' Test on a given tree to achieve aggregation of leaves
#'
#' This function sequentially tests on a tree in a top-down manner. The testing result determines aggregation of the leaves and makes sure the False Split Rate is controlled under a target level.
#' @template tree-template
#' @param p_vals A length-\code{num_interior_nodes} vector of p-values for the interior nodes. The ordering of p-values should correspond to the ordering of nodes in the hclust/hc_list object, i.e., the ith item of the vector is the p-value of the ith node. If a dendrogram is provided that stores the tree structure, the p-values should be ordered reversely as how each node of a dendrogram is reached by a recursive algorithm (See Examples for an example using dendrogram).
#' @param alpha A use-specified target FSR level
#' @param independent Whether the p-values are independent (default = TRUE).
#' @return Returns the testing result and calculated threshold values.
#' \item{alpha}{The target FSR level.}
#' \item{rejections}{A length-\code{num_interior_nodes} vector indicating whether each node is rejected.}
#' \item{threshold_functions}{A length-(\code{nnodes}-\code{nleaves}) vector of computed threshold value at each node. NA if not tested.}
#' \item{groups}{A length-\code{n-feature} vector of integers indicating the cluster to which each leaf is allocated.}
#' @examples
#' set.seed(123)
#' ## Example 1: Test with an hclust object
#' hc = hclust(dist((1:10) + runif(10)/10), method = "complete")
#' p_vals = c(runif(5), rbeta(4, 1, 60))
#' hierarchical_test(tree = hc, p_vals = p_vals, alpha = 0.3, independent = TRUE)
#' ## Example 2: Test with a dendrogram object
#' dend = as.dendrogram(hc)
#' p_vals = c(runif(4), rbeta(1, 1, 60), runif(1), rbeta(3, 1, 60))
#' hierarchical_test(tree = dend, p_vals = p_vals, alpha = 0.3, independent = TRUE)
#' ## Example 3: Test with a hc_list object
#' hc_list = dend_as_hclist(dend)
#' p_vals = c(runif(4), rbeta(1, 1, 60), runif(1), rbeta(3, 1, 60))
#' hierarchical_test(tree = hc_list, p_vals = p_vals, alpha = 0.3, independent = TRUE)
#' @export
hierarchical_test = function(tree = NULL, p_vals, alpha, independent = TRUE){

  # transform dendrogram to a list of length |T\L| (number of interior nodes). Item i in the list stores the children node of node i.
  if(class(tree) == "dendrogram"){
    dah = dend_as_hclist(tree)
    hc_list = dah
  }else if(class(tree) == "hclust"){
    dend = as.dendrogram(tree)
    dah = dend_as_hclist(dend)
    hc_list = dah
  }else if(class(tree) == "list"){
    hc_list = tree
  }else{
    stop("No proper tree structure is provided.")
  }


  n_interiornodes = length(hc_list) # number of interior nodes
  p = -min(unlist(hc_list)) # number of leaves
  m = n_interiornodes + p
  all_nodes = 1: m #including root node and interior nodes

  # Initialize rejections (the output of the algorithm.)
  rejections = rep(FALSE, m)
  threshold_functions = c(rep(NA, m-1), 1)

  # Find depth of each node (root has depth 1).
  depths = find_depths_in_list(hc_list)
  max_depth = max(depths)

  # Find the number of leaves of subtree rooted at each node
  leaves = number_leaves_in_list(hc_list)

  # Find degree_T and initilize degree_rej for each node
  degrees = c(rep(0,p),sapply(hc_list, function(x) length(x)))

  # upper bound of degree of the tree
  Deltas = find_max_degree_in_list(hc_list)
  Delta = max(Deltas)

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
      # We always reject the root node
      rejections[nodes_depth_d] = TRUE
      # Initialize R^{1:1}
      R_to_d_minus_1 = degrees[nodes_depth_d] - 1
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
          ifrejected = make_rejection(d, nodes_depth_d, nodes_to_test_depth_d, alpha, reshaping_func,
                                      p, D=max_depth-1, Delta, p_vals,
                                      leaves, depths, degrees,
                                      rejections, R_to_d_minus_1)
          # Update R^{1:(d-1)}
          R_to_d_minus_1 = ifrejected$R_to_d_minus_1

          # determine who gets rejected
          rejected_nodes_depth_d = nodes_to_test_depth_d[ifrejected$rejected_d]
          if(length(rejected_nodes_depth_d)!=0){
            # update rejections
            rejections[rejected_nodes_depth_d] = TRUE
          }
          # update thresholds
          threshold_functions[nodes_to_test_depth_d] = ifrejected$crit_func_d
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
  return(list(alpha = alpha,
              rejections = rejections,
              threshold_functions = threshold_functions[-(1:p)],
              groups = groups))
}
