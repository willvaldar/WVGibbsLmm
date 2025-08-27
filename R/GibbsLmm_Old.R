

rinv.chisq <- function(n, scale, df){
  scale/rchisq(n=n, df=df)
}


lmmgibbs <- function(...){
  g <- GibbsLmm$new(...)
  g
}

GibbsLmm <- setRefClass("GibbsLmm",
  fields=list(
    y             = "numeric",
    prior.beta    = "list",
    prior.sigma2  = "list",
    nobs          = "integer",
    nrand         = "integer",
    fixed         = "list",
    random        = "list",
    block         = "list",
    state         = "list")
)


GibbsLmm$methods(initialize =
  function(
      y,
      X             = matrix(rep(1, length(y)), ncol=1), # default intercept
      random        = list(),
      init.beta     = rep(0,ncol(X)),
      init.sigma2   = sd(y, na.rm=TRUE),
      prior.beta    = list(
          mean=0,
          varcov=1e3),
      prior.sigma2  = list(
          scale=0.02,
          df=0.02),
      fixef.blocking=TRUE,
      ...)
  {
    y <<- c(y)
    nobs <<- length(y)

    stopifnot(nrow(X)==nobs)
    
    ## Fixed Effects
    # basics
    fixed$X <<- X
    fixed$k <<- ncol(X)
        
    # priors
    fixed$prior.beta <<- prior.beta
    fixed$prior.beta$mean <<- matrix(nrow=fixed$k, rep(prior.beta$mean, length.out=fixed$k))
    if (1==fixed$k) {
      fixed$prior.beta$varcov <<- matrix(prior.beta$varcov)
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
    # state
    state$fixed$beta <<- init.beta
    
    ## Noise    
    # priors
    prior.sigma2 <<- prior.sigma2
    state$sigma2 <<- init.sigma2
    
    ## Random effects
    nrand <<- length(random)
    state$random <<- list()
    for (r in seq_along(random)){
      rd <- random[[r]]
      rname <- names(random)[r]
      # simple factor-type specification
      if (!is.list(rd)){
        if (is.factor(rd)){
          rd   <- list(factor=rd)
        } else {
          stop("Illegal object type of random effect ", rname,
               "; it should be a list or a factor.")
        }
      }
      if (!is.null(rd$factor)){
        rd$Z <- incidence.matrix(rd$factor)
      }
      # use Z and R specifications
      if (is.null(rd$Z) && is.null(rd$R)) {
        stop("Must specify Z or R for random effect ", rname, "\n")
      }
      if (!is.null(rd$Z) && is.null(rd$R)) {  # Z only
        if (nrow(rd$Z)!=nobs) {
          stop("Z should have n rows for random effect ", rname)
        }
        rd$k <- ncol(rd$Z)
        rd$R <- diag(rep(1,rd$k))
      } else if (is.null(rd$Z) && !is.null(rd$R)) { # R only
        if (length(y)!=ncol(rd$R)){
          stop("R must be n x n if Z is unspecified for random effect ", rname, "\n")
        }
        rd$k <- length(y)
        rd$Z <- diag(length(y))
      } else { # Z and R
        if (ncol(rd$Z)!=ncol(rd$R)){
          stop("Z and R sizes are incompatible for random effect", rname, "\n")
        }
        rd$k <- ncol(rd$Z)
      }
      if (!isSymmetric(rd$R)){
        stop("R must be a symmmetric matrix for random effect", rname, "\n")
      }

      # dimension and correlation matrix      
      rd$Rinv = solve(rd$R)
      
      ## blocking options
      # allow to be in a block update by default
      if (is.null(rd$blocking)) {
        rd$blocking = TRUE
      }
      # optimizations for non-block updating
      if (!rd$blocking) {
        rd$Zt  = t(rd$Z)
        rd$ZtZ = t(rd$Z) %*% rd$Z
      }
    
      # priors
      if (is.null(rd$prior.tau2)) {
        rd$prior.tau2=prior.sigma2
      }
      random <<- c(.self$random, list(rd))
    
      # state
      state$random[[r]] <<- list(
          tau2=1,
          ranef=rep(0,rd$k))
    }
    names(state$random) <<- names(random)
    
    # optimizations
    if (!fixef.blocking) {
      fixed$Prec0     <<- solve(fixed$prior$varcov)
      fixed$Prec0.Mu0 <<- fixed$Prec0 %*% fixed$prior$mean
      fixed$Xt        <<- t(fixed$X)
      fixed$XtX       <<- fixed$Xt%*%fixed$X
    }

    ## Construct a block of effects, to be updated as one.
    # Members of the block, listed in block$contents, include
    #  0: fixed effects
    #  1: random effect 1
    #  2: random effect 2 ... etc ...
    
    # make a list of who is in the block
    block$contents <<- c(
        ifow(fixef.blocking, 0, NULL),
        which(sapply(.self$random, function(x){x$blocking})))

    if (fixef.blocking) {
      block$Mu0 <<- fixed$prior$mean
      block$X   <<- fixed$X
    }
    for (r in setdiff(block$contents, 0)) { # add each participating random effect to the block
      rd        <-  .self$random[[r]]
      block$Mu0 <<- c(block$Mu0, rep(0, rd$k))
      block$X   <<- cbind(block$X, .self$random[[r]]$Z)
    }
    block$Xt  <<- t(block$X)
    block$XtX <<- block$Xt %*% block$X
    block$SigmaR <<- as.matrix(
      Matrix::bdiag(
        c(
          ifow(fixef.blocking, fixed$prior$varcov, NULL),
          lapply(.self$random[block$contents], function(x){x$R})  # or is it Rinv?
        )
      )
    )

    callSuper(...)
  }
)

GibbsLmm$methods(sample =
  function(
    n.iter=1,
    burnin=0,
    thin=1,
    verbose.at=round(n.iter/10),
    wanted=NULL,
    extractor=NULL,
    as.coda.object=TRUE
    )
  {
    'Performs one round of Gibbs updating
    '
    # prepare extraction of state output
    nsave <- round((n.iter-burnin)/thin)
    mcmc.mat <- NULL
    if (!is.null(wanted)) {
      extractor <- .self$make.param.extractor(patterns=wanted)
    }
    if (!is.null(extractor)) {
      eg.output <- extractor(state)
      mcmc.mat <- matrix(NA, ncol=length(eg.output), nrow=nsave)
      colnames(mcmc.mat) <- names(eg.output)
    }

    isave <- 0
    for (i in 1:n.iter) {
      ## effects
      if (!(0 %in% block$contents)) {
        .self$update.fixef()
      }      
      .self$update.block()
      for (r in setdiff(seq_along(random), block$contents)) {
        .self$update.ranef(r)
      }
      
      ## variances
      .self$update.sigma2()
      for (r in seq_along(random)) {
        .self$update.tau2(r)
      }

      # recording values
      if (!is.null(extractor) && i > burnin && 0==i%%thin) {
        isave=isave+1
        mcmc.mat[isave,] = extractor(state)
      }
      if (0<verbose.at && 0==i%%verbose.at) {
        cat("[",i,"]", sep="")
      }
    }
    if (0<verbose.at) cat("\n")

    if (!is.null(extractor) && as.coda.object) {
        mcmc.mat=as.mcmc(mcmc.mat)
    }
    invisible(mcmc.mat)
  }
)

GibbsLmm$methods(make.param.extractor =
  function(patterns=c(), all=FALSE)
  {
    'Returns a function that will extract specified parameters
    from the state vector
    '
    params = names(unlist(state))
    wanted = integer(0)
    if (all) {
      wanted = 1:length(params)
    } else {
      for (p in patterns) {
        wanted = c(wanted, grep(p, params))
      }
    }
    function(s) {
      unlist(s)[wanted]
    }
  }
)    

GibbsLmm$methods(update.fixef =
  function()
  {
    'Updates the fixed effects vector
    using its full conditional distribution
    '
    yres = .self$yresid(exclude=0)
    Prec1 = fixed$XtX/state$sigma2 + fixed$Prec0

    ## slow step for large matrices
    # Sigma1 = solve(Prec1) # slower
    Sigma1 = chol2inv(chol(Prec1)) # faster
    
    Mu1.unscaled = fixed$Xt %*% yres/state$sigma2 + fixed$Prec0.Mu0
    Mu1 = Sigma1%*%Mu1.unscaled
    state$fixed$beta <<- c(rmnorm(mean=Mu1, varcov=Sigma1))
  }
)

GibbsLmm$methods(update.block =
  function()
  {
    'Updates one or more components that have been concatenated
    together as a block, from their block full conditional distribution.
    For example, if block$contents contains 0, 1 and 2, then a
    long vector containing fixed effects and random effects for 
    random components 1 and 2 are updated in a single draw.
    '
    yres = .self$yresid(exclude=block$contents)
    m <- NULL
    if (0 %in% block$contents){ # if fixed effects are in the block
      m <- rep(1, fixed$k)
    }
    for (r in setdiff(block$contents, 0)) {
      m <- c(m, rep(state$random[[r]]$tau2, times=random[[r]]$k))
    }
    Sigma <- block$SigmaR %*% diag(m)
    Prec0 <- solve(Sigma)
    Prec1 <- block$XtX/state$sigma2 + Prec0

    ## slow step for large matrices
    # Sigma1 = solve(Prec1) # slower
    Sigma1 = chol2inv(chol(Prec1)) # faster
    
    Mu1.unscaled = block$Xt %*% yres/state$sigma2 + Prec0 %*% block$Mu0
    Mu1 = Sigma1%*%Mu1.unscaled
    u = c(rmnorm(mean=Mu1, varcov=Sigma1))
    
    # allocate effects to state
    i = 1:fixed$k
    state$fixed$beta <<- u[i]
    u=u[-i]
    for (r in setdiff(block$contents, 0))
    {
      i = 1:random[[r]]$k
      state$random[[r]]$ranef <<- u[i]
      u = u[-i]
    }
  }
)

GibbsLmm$methods(update.ranef =
  function(r)
  {
    'Updates effects vector for specified random component
    using its full conditional distribution
    '
    yres = .self$yresid(exclude=r)
    rstate = state$random[[r]]
    rd = random[[r]]
    Prec1 = rd$ZtZ/state$sigma2 + rd$Rinv/rstate$tau2
    
    ## slowest step for large matrices
    # Sigma1 = solve(Prec1) # 20s/20
    Sigma1 = chol2inv(chol(Prec1))    # 11s/20
    
    Mu1.unscaled = rd$Zt %*% yres / state$sigma2
    Mu1= Sigma1 %*% Mu1.unscaled
    
    ## 2nd slowest step
    state$random[[r]]$ranef <<- c(rmnorm(mean=Mu1, varcov=Sigma1))    # 5s/20
  }
)

GibbsLmm$methods(update.sigma2 =
  function()
  {
    'Updates sigma2 noise variance using its full conditional distribution
    '
    yres = .self$yresid()
    m1 = c( sum(yres^2) + prior.sigma2$scale )
    nu1 = nobs + prior.sigma2$df
    state$sigma2 <<- rinv.chisq(n=1, scale=m1, df=nu1)
  }
)

GibbsLmm$methods(update.tau2 =
  function(r)
  {
    'Updates variance tau2 for the specified random component
    using its full conditional distribution
    ' 
    rstate = state$random[[r]]
    rd = random[[r]]
    m1 = c((t(rstate$ranef) %*% rd$Rinv %*% rstate$ranef) + rd$prior.tau2$scale)
    nu1 = rd$k + rd$prior.tau2$df
    state$random[[r]]$tau2 <<- rinv.chisq(n=1, scale=m1, df=nu1)
  }
)

GibbsLmm$methods(yhat =
  function(exclude=NULL)
  {
    'Returns the predicted value of y given all but the excluded components,
    where 0 denotes the fixed effects, and 1 or greater denotes the corresponding
    random component
    '
    yhat = numeric(nobs)
    if (!(0 %in% exclude))
    {
      yhat=fixed$X %*% state$fixed$beta
    }
    for (r in setdiff(seq_along(random),exclude))
    {
      yhat = yhat + random[[r]]$Z %*% state$random[[r]]$ranef
    }
    yhat
  }
)

GibbsLmm$methods(yresid =
  function(exclude=NULL)
  {
    'Returns the residual of y after prediction based on all components except
    those specified as excluded -- see method yhat()
    '
    y - .self$yhat(exclude)
  }
)
















