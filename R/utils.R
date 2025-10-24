

rmvnormAB <- function(B, cholA = NULL, A = NULL, stz="off") {
  # Either chol(A) or A must be provided
  if (is.null(cholA) && is.null(A)) stop("Either cholA or A must be provided")

  if (!(stz %in% c("off", "project", "basis"))){
    stop("Unknown sum-to-zero constraint type `", stz, "`")
  }
  if ("basis"==stz){
    return (rmvnormAB_stz_basis(B, cholA=cholA, A=A))
  }

  # Compute Cholesky if needed
  if (is.null(cholA)) {
    if (!is.matrix(A)) A <- as.matrix(A)
    cholA <- chol(A)  # Upper-triangular
  }
  
  n <- length(B)
  # Solve for mean: mu = A^{-1} B using triangular solves
  y <- forwardsolve(t(cholA), B)   # Solve t(cholA) y = B
  mu <- backsolve(cholA, y)        # Solve cholA mu = y
  
  # Sample from standard normal
  z <- rnorm(n)
  
  # Transform to get covariance A^{-1}: x = mu + A^{-1/2} z
  w <- forwardsolve(t(cholA), z)   # Solve t(cholA) w = z
  x <- mu + w

  if ("project"==stz){
    x <- x - mean(x)
  }
  return(x)
}

#' Construct C matrix (orthonormal basis of sum-to-zero subspace)
make.orthobasis <- function(k){
  J <- matrix(1, k, k)
  M <- diag(k) - J / k # make centering matrix
  Q <- qr.Q(qr(M))  # orthonormal basis
  C <- Q[, 1:(k - 1)]  # drop final column
  C
}

rmvnormAB_stz_basis <- function(B, cholA = NULL, A = NULL) {
  k <- length(B)
  C <- make.orthobasis(k)
  one <- rep(1, k)

  if (is.null(A)) {
    A <- t(cholA) %*% cholA
  }

  # Compute unconstrained mean and precision
  mu <- solve(A, B)

  # Project mean to constrained space
  mu_stz <- mu - mean(mu)

  # Project precision to subspace
  A_tilde <- t(C) %*% A %*% C
  cholA_tilde <- chol(A_tilde)

  # Compute projected mean in subspace
  B_tilde <- t(C) %*% A %*% mu_stz  # key fix here
  y <- forwardsolve(t(cholA_tilde), B_tilde)
  mu_v <- backsolve(cholA_tilde, y)

  # Sample
  z <- rnorm(k - 1)
  w <- forwardsolve(t(cholA_tilde), z)
  v <- mu_v + w

  # Return to original space
  u <- C %*% v
  return(as.numeric(u))
}

classify.corrmat <- function(R){
  if (!isSymmetric(R)){ return("nonsymmetric") }
  if (all(R==diag(ncol(R)))){ return("I") }
  if (all(R==diag(diag(R)))){ return("diag")}
  return("full")
}

verify.corrmat <- function(R, size=NULL){
  stopifnot(inherits(R, "matrix"))
  type <- classify.corrmat(R)
  if ("nonsymmetric"==type){
    stop("corrmat must be symmetric")
  }
  if (!is.null(size)){
    stopifnot(NROW(R)==size)
  }
  type
}
