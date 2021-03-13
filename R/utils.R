
# This function output a vector indicating the depth of each node(including interior node and leaf node)
# The depths are numbered from root to leaves
# It returns a vector of 2*p-1 where the first p are leaf nodes, and the next p-1 are interior nodes
find_depths_in_list = function(hc_list){
  p = -min(unlist(hc_list))
  m = length(hc_list) + p
  all_nodes = 1:m
  depths = rep(1, m)
  for(i in (m-1):1)
  {
    depths[i] = depths[find_parent_in_list(i, hc_list)]+1
  }

  return(depths)
}


# This function outputs a sequence of number of leaves in the subtree rooted at each node
# It returns a vector of |T| where the first p entries correspond to
# leaf nodes, and the rest correspond to interior nodes
number_leaves_in_list = function(hc_list){
  p = -min(unlist(hc_list))
  m = p + length(hc_list)
  number_of_leaves = rep(1, m)
  number_of_leaves[(p+1):m] = sapply(1:length(hc_list),function(x){
    length(find_leaves_in_list(x,hc_list)) })
  return(number_of_leaves)
}


# This function finds the children of each node
add_child_interior_in_list = function(hc_list, current){
  # current: the index of a node
  # Returns the indices of all interior nodes descended from current, including
  # current.

  if(current<=0)
  {
    return(c())
  }
  if(current>0){
    v = c(current)
    for(i in hc_list[[current]]){
      parts = add_child_interior_in_list(hc_list, i)
      v = unique(c(v, parts))
    }

    return(sort(unique(c(v, current))))
  }
}


# This function outputs a vector indicating the number of descendents of each node
# including interior node and leaf node, and the node itself
# It returns a vector of |T| where the first p entries correspond to
# leaf nodes, and the rest correspond to interior nodes
number_descendents_in_list = function(hc_list){
  p = -min(unlist(hc_list))
  m = p + length(hc_list)
  number_of_descendents = rep(0, m)

  number_of_descendents[(p+1):m] = sapply(1:(m-p), function(x) {
    length(add_child_interior_in_list(hc_list, x))
  })
  return(number_of_descendents)

}


# This function gives the parent node of a given node i
# The nodes are numbered from bottom to top
# where the first p are leaf nodes, and then are interior nodes
find_parent_in_list = function(i, hc_list){
  p = -min(unlist(hc_list))
  m = p + length(hc_list)
  stopifnot(i>=1 & i<=m-1)

  if(i > p){
    return(which(sapply(hc_list, function(x){any(x %in% (i-p))})) + p)
  }else{
    return(which(sapply(hc_list, function(x){any(x %in% (-i))})) + p)
  }
}

# This function calculates 1/(a+1) + ... + 1/b
# b, a can be vectors of the same lengths
harmonic_diff = function(b, a){
  return(digamma(b+1) - digamma(a+1))
}



# This function constructs a binary matrix from a hc_list object
# The binary matrix encodes the ancestor-descendnat relationships between leaves and nodes
# as defined in Yan and Bien (2018)
tree.list.matrix = function(hc_list)
{
  p <- -min(unlist(hc_list))
  n_interior <- length(hc_list)
  A_i <- c(as.list(seq(p)), sapply(seq(n_interior), function(x) find_leaves_in_list(x,
                                                                                    hc_list)))
  A_j <- sapply(seq(length(A_i)), function(x) rep(x, len = length(A_i[[x]])))
  A <- Matrix::sparseMatrix(i = unlist(A_i), j = unlist(A_j), x = rep(1, len = length(unlist(A_i))))
  A
}




# This function determine aggregation based on rejection/non-null
determine_aggregation = function(hc_list, rejections){
  p = -min(unlist(hc_list))
  our_rejected_nodes = which(rejections == TRUE)

  theta_agg = rep(0, p)
  for(u in our_rejected_nodes){
    for(k in hc_list[[u]]){
      if(k<0){
        theta_agg[abs(k)] = theta_agg[abs(k)] + stats::runif(1)
      }else{
        theta_agg[find_leaves_in_list(k, hc_list)] = theta_agg[find_leaves_in_list(k, hc_list)] + stats::runif(1)
      }
    }
  }
  return(as.numeric(factor(theta_agg)))
}


# This function find the desendant interior nodes of a node
find_descendants_interior_in_list = function(hc_list, current){
  # hc_list: A list object that encodes the hierarchical tree structure
  # current: the index of a node

  # Returns the indices of all interior nodes descended from current, including
  # current.

  if(current<=0)
  {
    return(c())
  }
  if(current>0){
    v = c(current)
    for(i in hc_list[[current]]){
      parts = add_child_interior_in_list(hc_list, i)
      v = unique(c(v, parts))
    }

    return(sort(unique(c(v, current))))
  }
}

# This function outputs a vector of maximum degree of the subtree rooted at each node
find_max_degree_in_list = function(hc_list){
  # It returns a vector of 2*p-1 where the first p are leaf nodes, and the next p-1 are interior nodes
  p = -min(unlist(hc_list))
  m = p + length(hc_list)
  Deltas = rep(1, m)
  Deltas[1:p] = 0

  Deltas[(p+1):m] = sapply(1:(m-p), function(x) {
    children = add_child_interior_in_list(hc_list, x)
    max(sapply(children, function(y){
      length(hc_list[[y]])
    }))
  })

  return(Deltas)
}


