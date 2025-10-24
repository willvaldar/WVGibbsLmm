
gamma.prior <- function(shape, rate){
  # TODO: add scale option
  list(shape=shape, rate=rate)
}

normal.prior <- function(k=1, mean=0, var=1, corrmat=NULL){
  # Check mean and corrmat
  k <- length(mean)
  if (is.null(corrmat)){
    corrmat <- diag(rep(1, k))
  }
  verify.corrmat(corrmat, size=k)
  list(mean=mean, var=var, corrmat=corrmat, stz="off")
}

#'
#'
normal.ranef <- function(mean=0,
                         Z=NULL,
                         corrmat=NULL,
                         init.ranef=NULL,
                         init.var=1,
                         prior=NULL,
                         stz="off"){
  # TODO: fill this out
  list(Z=Z, corrmat=corrmat, stz=stz)
}


parse.varcomp.arg <- function(rname, rd, default.prior, nobs){
  cat("Parsing random effect ", rname, "\n")
  if (!is.list(rd)){
    if (is.factor(rd)){
      rd   <- list(factor=rd)
    } else {
      stop("Illegal object type of random effect in ", rname,
           "it should be a list or a factor.")
    }
  }
  if (!is.null(rd$factor)){
    rd$Z <- incidence.matrix(rd$factor)
  }
  if (is.null(rd$Z) && is.null(rd$R)) {
    stop("Must specify Z or R for random effect ", rname, "\n")
  }
  # use Z and R specifications
  if (!is.null(rd$Z) & is.null(rd$R)) {  # Z only
    if (nrow(rd$Z)!=nobs) {
      stop("Z should have n rows for random effect ", rname)
    }
    rd$R <- diag(rep(1, ncol(rd$Z))) # R=I
  } else if (is.null(rd$Z) && !is.null(rd$R)) { # R only
    if (ncol(rd$R)!=nobs){
      stop("R must be n x n if Z is unspecified for random effect ", rname, "\n")
    }
    rd$Z <- diag(rep(1, nobs))
  } else { # Z and R
    if (nrow(rd$Z)!=nobs) {
      stop("Z should have n rows for random effect ", rname)
    }
    if (ncol(rd$Z)!=ncol(rd$R)){
      stop("Z and R sizes are incompatible for random effect", rname, "\n")
    }
  }
  Rtype <- classify.corrmat(rd$R)
  if ("nonsymmetric"==Rtype){
    stop("R must be a symmmetric matrix for random effect", rname, "\n")
  }
  if (is.null(rd$prior)){
    rd$prior.tau2inv <- default.prior
  }
  return (rd)
}


# not currently used
parse.fixedeffects.arg <- function(X, prior.beta){
  k <- ncol(X)
  if (is.null(prior.beta$mean)){
    stop("Must specify prior mean")
  }
  prior.beta$mean <- matrix(nrow=k, rep(prior.beta$mean, length.out=k))
                            
  if (1==k) {
    prior.beta$corrmat <- matrix(prior.beta$corrmat)
  } else { 
    if (length(fixed$prior.beta$varcov)==1) {
      fixed$prior.beta$varcov <<- rep(fixed$prior.beta$varcov, length.out=fixed$k)
    }
    if (length(fixed$prior.beta$varcov)==fixed$k) {
      fixed$prior.beta$varcov <<- diag(fixed$prior.beta$varcov)
    }
    if (!all(dim(fixed$prior.beta$varcov)==rep(fixed$k,2))) {
      stop("Bad prior for beta\n")
    }
  }
  browser()
  fixed <- GibbsLmm_FixedEffects$new(
    X = X,
    R = prior.beta$varcov,
    tau2inv = prior.beta$tau2inv)
  fixed
}


lmmgibbs <- function(y,
                     X,
                     V = NULL,
                     prior.beta,
                     prior.sigma2inv,
                     prior.tau2inv = NULL,
                     random = NULL,
                     Z.updater = NULL,
                     ...){
  if (any(is.na(y))){
    stop("y cannot have missing values")
  }
  if (!is.matrix(X)){
    stopifnot(length(X)==length(y))
    X <- as.matrix(X, ncol=1)
  } else {
    stopifnot(nrow(X)==length(y))
  }
  if (is.null(V)){
    V <- diag(rep(1, length(y)))
  }
  if (length(random) > 0 & is.null(prior.tau2inv)){
    stop("Must specify prior for random effects model")
  }
  # check/reformat random effects
  for (r in seq_along(random)){
    rd <- random[[r]]
    rname <- names(random)[r]
    rd <- parse.varcomp.arg(rname, rd, default.prior=prior.tau2inv, nobs=length(y))
    random[[r]] <- rd
  }
  g <- GibbsLmm$new(y=y,
                    X=X,
                    V=V,
                    prior.beta = prior.beta,
                    prior.sigma2inv = prior.sigma2inv,
                    random = random,
                    Z.updater = Z.updater,
                    ...)
  g
}


