# Solve LP for constant lambda_n in the paper
# Use lpSolve package
optimal_mu_lp = function(p, hbetaGA, Sigmahat_hbetaGA, constr_matrix, constr_dir){
  # Sigma: Sample Covariance matrix
  # hbetaGA: hbeta_G %*% A

  constr_matrix[p+1,] = c(Sigmahat_hbetaGA, -1)
  constr_matrix[2*p+2,] = c(-Sigmahat_hbetaGA, -1)

  temp <- c(as.numeric(hbetaGA), rep(0, p-length(hbetaGA)), norm(hbetaGA, "2"))
  contr_rhs = c(temp, -temp)

  result = lpSolve::lp(direction = "min", objective.in = c(c(rep(0, p),1),c(rep(0, p),-1)),
              const.mat = cbind(constr_matrix, -constr_matrix),
              const.dir = constr_dir,
              const.rhs = contr_rhs)
  return(result$objval)
}

optimal_mu = function(Sigma, hbetaGA){
  # solve LP for constant lambda_n in the paper (Use CVXR package)

  # Sigma: Sample Covariance matrix
  # hbetaGA: hbeta_G %*% A

  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package \"CVXR\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  dim_Sigma = dim(Sigma)
  p = dim_Sigma[1]
  # find the group of vectors {e_1, .., e_p, (hbeta_G * A, 0)^T }
  if(norm(hbetaGA, "2")!=0){
    C_vectors = rbind(diag(1, p), c(as.numeric(hbetaGA),
                                    rep(0, p-length(hbetaGA)))/norm(hbetaGA, "2"))
  }else{
    C_vectors = rbind(diag(1, p), rep(0, p))
  }

  # solve linear optimization
  m = CVXR::Variable(p)
  objective <- CVXR::Minimize(max(abs(C_vectors %*% (Sigma%*%m - c(as.numeric(hbetaGA),
                                                             rep(0, p-length(hbetaGA)))))))
  problem = CVXR::Problem(objective)
  result <- CVXR::solve(problem)
  return(result$value)
}


pvalue_node_step = function(u, leaves_index, hbeta, Sigma_hat, Dmat,
                            constr_matrix, constr_dir,
                            constraint_matrix, X, y, hsigma, n, p, t){
  # This function computes p-value for a single node
  n_leaves = length(leaves_index)
  # Use the projection matrix
  D_matrix = diag(1, n_leaves) - 1/n_leaves * matrix(1, n_leaves, n_leaves)

  if(n_leaves == 2){
    A = D_matrix %*% t(D_matrix)
  }else{
    A = crossprod(D_matrix)
  }

  hbeta_G = hbeta[leaves_index]
  hbetaGA = hbeta_G %*% A
  norm_hbetaGA = norm(hbetaGA, "2")
  if(norm_hbetaGA>0){

    Sigmahat_hbetaGA =  Sigma_hat %*%
      c(as.numeric(hbetaGA), rep(0, p-length(hbetaGA)))/norm_hbetaGA

    # solve lp for mu
    mu = optimal_mu_lp(p, hbetaGA, Sigmahat_hbetaGA, constr_matrix, constr_dir)
    if(mu == 0){
      mu = optimal_mu(Sigma_hat, hbetaGA)
    }

    if(mu<10^(-15)){
      # avoid problem caused by floating precision
      mu = 0.0001
    }

    # relax the bound
    res <- 1.1
    mu2 = mu *res

    # Find the optimal projection direction u_A
    constraint_matrix[,p+1] = Sigmahat_hbetaGA
    constraint_matrix[,2*p+2] = -Sigmahat_hbetaGA

    temp <- c(as.numeric(hbetaGA), rep(0, p-length(hbetaGA)), norm(hbetaGA, "2"))
    constraint_vector = c(temp - mu2, -temp - mu2)

    result = quadprog::solve.QP(Dmat = Dmat, dvec = rep(0, p),
                                Amat = constraint_matrix, bvec = constraint_vector)

    proj_uA = result$solution
    optvalue = result$value


  }else{
    mu = 0
    proj_uA = rep(0,p)
    optvalue = 0

  }

  Q_A = sum(hbetaGA * hbeta_G) + 2/n* sum(proj_uA * crossprod(X, y - X %*% hbeta))
  V_A = 4* hsigma^2/n * optvalue + t/n
  p_val = 2*(1- pnorm(as.numeric(abs(Q_A)), mean = 0, sd = as.numeric(sqrt(V_A))))
  return(p_val)
}

pvalue_group_all = function(hc_list, X, Sigma_hat, Dmat=NULL, hbeta, hsigma, y, t=1){
  # hc_list: list length of $T\L$ for storing tree structure
  # X: n by p design matix
  # Sigma_hat: sample covariance matrix
  # D_mat: nearPD(Sigma_hat)
  # hbeta: the initial estimator of true beta
  # y: response vector
  # t: regularizing constant taking values like 0, 1 in the paper

  n = dim(X)[1]
  p = dim(X)[2]
  if(is.null(Dmat)){
    Dmat = Matrix::nearPD(Sigma_hat)
  }

  constraint_matrix = cbind(Sigma_hat, rep(0,p), -Sigma_hat, -rep(0,p))
  constr_matrix = rbind(cbind(rbind(Sigma_hat, rep(0,p)), rep(-1, p+1)),
                        cbind(-rbind(Sigma_hat, rep(0,p)), rep(-1, p+1)))
  constr_dir = rep("<=", 2*p+2)

  n_interiornodes = length(hc_list)
  ps = rep(0, n_interiornodes)
  for(u in 1:n_interiornodes){
    leaves_index = find_leaves_in_list(u, hc_list)

    ps[u] = pvalue_node_step(u, leaves_index, hbeta, Sigma_hat, Dmat,
                             constr_matrix, constr_dir,
                             constraint_matrix, X, y, hsigma, n, p, t)
  }

  return(ps = ps)
}



#' Aggregate rare features hierarchically with False Split Rate control
#'
#' This function achieves rare features aggregation while simultaneously controlling False Split Rate under a target level.
#' The aggregation is achieved by two steps: (1) Generate p-values for each interior node (2) Sequentially test on the tree.
#' @param y A length-\code{n-observation} response vector that guides the aggregation
#' @param X An \code{n-observation}-by-\code{n-feature} design matrix. Each row corresponds to a subject and each column stores the observation of a feature made by subjects.
#' @param sigma Standard deviation of noise. If not given, the algorithm will estimate sigma.
#' @template tree-template
#' @param alpha A use-specified target FSR level
#' @return Returns the aggregation result.
#' \item{alpha}{The target FSR level.}
#' \item{groups}{A length-\code{n-feature} vector of integers indicating the cluster to which each feature is allocated.}
#' \item{rejections}{A length-(\code{num_interior_nodes}) vector indicating whether each node is rejected.}
#' @examples
#' ## See vignette for a small data example.
#' @importFrom stats pnorm runif sd
#' @export
aggregate_features = function(y, X, sigma = NULL, tree= NULL, alpha){

  #center response and design matrix
  y = as.vector(y)
  y = y - mean(y)
  X = as.matrix(apply(X, 2, function(x) {
    if(sd(x) != 0){
      (x - mean(x, na.rm = T))/sd(x)
    }else{
      x
    }
  }))


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


  # Make A_matrix
  A_matrix = tree.list.matrix(hc_list)

  # Initial estimate
  fit_rare = rare::rarefit(y = y, X=X, A = A_matrix, nalpha = 5, intercept = FALSE)
  fit_rare.cv = rare::rarefit.cv(fit_rare, y=y, X, nfolds = 5)
  hbeta = fit_rare$beta[[fit_rare.cv$ibest[2]]][,fit_rare.cv$ibest[1]]

  Sigma_hat = 1/nrow(X) * crossprod(X)
  Dmat <- Matrix::nearPD(Sigma_hat)$mat

  # If sigma unkonwn
  if(is.null(sigma)){
    if (!requireNamespace("scalreg", quietly = TRUE)) {
      stop("Package \"scalreg\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    object = scalreg::scalreg(X, y)
    sigma = object$hsigma
  }

  p_vals = pvalue_group_all(hc_list, X, Sigma_hat, Dmat, hbeta, sigma, y, t=1)
  result = hierarchical_test(tree = hc_list, p_vals= p_vals, alpha = alpha, independent = FALSE)
  groups = result$groups
  rejections = result$rejections

  return(list(alpha = alpha, groups = groups, rejections = rejections))

}
