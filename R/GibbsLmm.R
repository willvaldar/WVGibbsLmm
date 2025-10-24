
  
## Variance component class
GibbsLmm_VarComp <- setRefClass(
  "GibbsLmm_VarComp",
  fields = list(
    # public
    Z              = "matrix",
    R              = "matrix",
    prior.tau2inv  = "list",
    Vinv           = "matrix",  # R_eps^{-1}
    stz            = "character",
    # added after construction
    ZtVinv         = "ANY",
    ZtVinvZ        = "ANY",
    Rinv           = "ANY"
  )
)

GibbsLmm_VarComp$methods(
  initialize = function(Z, R, prior.tau2inv, Vinv, stz, ...){
    ## assign arguments safely to fields
    stopifnot(is.matrix(Z), is.matrix(R), is.matrix(Vinv))
    .self$Z       <- Z
    .self$R       <- R
    .self$prior.tau2inv <- prior.tau2inv
    .self$Vinv    <- Vinv
    .self$stz     <- stz
    .self$Rinv    <- solve(R)
    .self$ZtVinv  <- NULL
    .self$ZtVinvZ <- NULL
  },
  getEffectNames = function(){
    colnames(.self$Z)
  },
  get.k = function(){
    ncol(.self$Z)
  },
  get.Rinv = function(){
    .self$Rinv
  },
  get.Z = function(){
    .self$Z
  },
  get.ZtVinv = function(){
    if (is.null(.self$ZtVinv)){
      .self$ZtVinv <- t(.self$Z) %*% .self$Vinv
    }
    .self$ZtVinv
  },
  get.ZtVinvZ = function(){
    if (is.null(.self$ZtVinvZ)){
      .self$ZtVinvZ <- .self$get.ZtVinv() %*% .self$Z
    }
    .self$ZtVinvZ
  },
  calc.init.tau2inv = function(){
    return (1)
  },
  calc.init.u = function(){
    u <- rep(0, times=.self$get.k())
    names(u) <- .self$getEffectNames()
    u
  },
  calc.post.tau2inv = function(u){
    k     <- ncol(.self$Z)
    shape <- 0.5*k + prior.tau2inv$shape
    rate  <- 0.5*crossprod(u, .self$Rinv %*% u) + prior.tau2inv$rate
    list(shape=shape, rate=rate)
  },
  calc.post.u.A = function(sigma2inv, tau2inv){
    stopifnot(!is.null(sigma2inv))
    stopifnot(!is.null(tau2inv))
    .self$get.ZtVinvZ()*sigma2inv + .self$Rinv*tau2inv
  },
  calc.post.u.B = function(sigma2inv, yres){
    .self$get.ZtVinv() %*% yres * sigma2inv
  },
  sample.post.u = function(sigma2inv, tau2inv, yres){
    A <- .self$calc.post.u.A(sigma2inv, tau2inv)
    B <- .self$calc.post.u.B(sigma2inv, yres)
    u <- c(rmvnormAB(A=A, B=B, stz=.self$stz))
    names(u) <- .self$getEffectNames()
    return (u)
  },
  set.Z = function(Z){
    if (ncol(Z) != ncol(.self$R)){
      stop("bad Z")
    }
    .self$Z       <- Z
    .self$ZtVinv  <- NULL
    .self$ZtVinvZ <- NULL
  }
)

## Fixed effects class - contains VarComp class
GibbsLmm_FixedEffects <- setRefClass(
  "GibbsLmm_FixedEffects",
  fields = list(
    vc      = "ANY", #GibbsLmm_VarComp", # use VarComp to hold X and R
    tau2inv = "numeric"
  )
)

GibbsLmm_FixedEffects$methods(
  initialize = function(X, R, tau2inv, Vinv, ...){
    .self$vc      <- GibbsLmm_VarComp$new(Z=X,
                                          R=R,
                                          Vinv=Vinv,
                                          prior.tau2inv=list(),
                                          stz="off") # stz doesn't make sense for mu
    .self$tau2inv <- tau2inv
  },
  getEffectNames = function(){
    .self$vc$getEffectNames()
  },
  get.X = function(){
    .self$vc$get.Z()
  },
  get.k = function(){
    .self$vc$get.k()
  },
  get.Rinv = function(){
    .self$vc$get.Rinv()
  },
  calc.post.beta.A = function(sigma2inv){
    .self$vc$get.ZtVinvZ()*sigma2inv + .self$vc$get.Rinv()*tau2inv
  },
  calc.post.beta.B = function(sigma2inv, yres){
    .self$vc$get.ZtVinv() %*% yres * sigma2inv
  },
  sample.post.beta = function(sigma2inv, yres){
    A <- .self$calc.post.beta.A(sigma2inv)
    B <- .self$calc.post.beta.B(sigma2inv, yres)
    beta <- c(rmvnormAB(A=A, B=B, stz="off"))
    names(beta) <- .self$getEffectNames()
    return (beta)
  }
)


GibbsLmm <- setRefClass("GibbsLmm",
  fields=list(
    y                = "numeric",
    prior.sigma2inv  = "list",
    Vinv             = "matrix",
    nobs             = "integer",
    nrand            = "integer",
    fixed            = "ANY", # GibbsLmm_FixedEffects object
    random           = "list", # list of GibbsLmm_VarComp objects
    state            = "list",
    Z.updater        = "ANY") # function to update Zs
)

GibbsLmm$methods(initialize =
  function(
      y,
      X,
      V,
      prior.beta,
      prior.sigma2inv,
      random  = list(),
      Z.updater = NULL,
      ...)
  {
    .self$y <- c(y)
    .self$nobs <- length(y)
    verify.corrmat(V, nobs)
    .self$Vinv <- solve(V)

    ## Fixed Effects
    # basics
    fixed <<- GibbsLmm_FixedEffects$new(
      X=X,
      R=prior.beta$corrmat,
      tau2inv=1/prior.beta$var,
      Vinv=Vinv)
    # state
    state$fixed$beta <<- rep(0, fixed$get.k())
    names(state$fixed$beta) <<- fixed$getEffectNames()
    
    ## Noise    
    # priors
    .self$prior.sigma2inv <<- prior.sigma2inv
    state$sigma2inv <<- 1/var(.self$y)

    ## Random effects
    nrand <<- length(random)
    state$random <<- list()
    for (r in seq_along(random)){
      rd <- random[[r]]
      rname <- names(random)[r]
      # structure and prior
      vc <- GibbsLmm_VarComp$new(
        Z = rd$Z,
        R = rd$R,
        prior.tau2inv = rd$prior.tau2inv,
        Vinv = Vinv,
        stz = rd$stz)
      .self$random[[rname]] <- vc
      # state  -- maybe should be generated by the vc?
      state$random[[r]] <<- list(
          tau2inv=vc$calc.init.tau2inv(),
          ranef=vc$calc.init.u())
    }
    names(state$random) <<- names(random)
    .self$Z.updater <<- Z.updater
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
    # Draw MCMC samples
    isave <- 0
    for (i in 1:n.iter) {
      ## sample effects
      .self$update.fixef()
      for (r in seq_along(random)) {
        .self$update.ranef(r)
      }
      ## sample variances
      .self$update.sigma2inv()
      for (r in seq_along(random)) {
        .self$update.tau2inv(r)
      }
      ## sample Z
      if (!is.null(.self$Z.updater)){
        .self$update.Z()
      }
      
      # recording values
      if (!is.null(extractor) && i > burnin && 0==i%%thin) {
        isave <- isave+1
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
    params <- names(unlist(state))
    wanted <- integer(0)
    if (all) {
      wanted <- 1:length(params)
    } else {
      for (p in patterns) {
        wanted <- c(wanted, grep(p, params))
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
    yres <- .self$yresid(exclude=0)
    beta <- fixed$sample.post.beta(state$sigma2inv, yres)
    state$fixed$beta <<- beta
  }
)

GibbsLmm$methods(update.ranef =
  function(r) {
    'Updates effects vector for specified random component
    using its full conditional distribution
    '
    yres <- .self$yresid(exclude=r)
    rstate <- state$random[[r]]
    vc <- random[[r]]
    u <- vc$sample.post.u(state$sigma2inv, rstate$tau2inv, yres)
    state$random[[r]]$ranef <<- u
  }
)

GibbsLmm$methods(update.sigma2inv =
  function()
  {
    'Updates sigma2 noise variance using its full conditional distribution
    '
    yres <- .self$yresid()
    shape <- .self$nobs/2 + prior.sigma2inv$shape
    rate <- crossprod(yres, .self$Vinv %*% yres) + prior.sigma2inv$rate
    state$sigma2inv <<- rgamma(1, shape=shape, rate=rate)
  }
)

GibbsLmm$methods(update.tau2inv =
  function(r)
  {
    'Updates variance tau2 for the specified random component
    using its full conditional distribution
    '
    rstate <- state$random[[r]]
    vc <- random[[r]]
    post <- vc$calc.post.tau2inv(state$random[[r]]$ranef)
    state$random[[r]]$tau2inv <<- rgamma(1, shape=post$shape, rate=post$rate)
  }
)

GibbsLmm$methods(yhat =
  function(exclude=NULL) {
    'Returns the predicted value of y given all but the excluded components,
    where 0 denotes the fixed effects, and 1 or greater denotes the corresponding
    random component
    '
    yhat = numeric(nobs)
    if (!(0 %in% exclude)) {
      yhat <- fixed$get.X() %*% state$fixed$beta
    }
    for (r in setdiff(seq_along(random), exclude)) {
      yhat <- yhat + random[[r]]$get.Z() %*% state$random[[r]]$ranef
    }
    yhat
  }
)

GibbsLmm$methods(yresid =
  function(exclude=NULL)  {
    'Returns the residual of y after prediction based on all components except
    those specified as excluded -- see method yhat()
    '
    y - .self$yhat(exclude)
  }
)

GibbsLmm$methods(update.Z =
  function() {                   
    which.r <- .self$Z.updater$whichComponents()
    yres <- .self$yresid(exclude=which.r)
    new.Z <- .self$Z.updater$sample.posterior.Z(yres=yres, Vinv=Vinv, state=state)
    for (r in which.r){
      .self$random[[r]]$set.Z( new.Z[[r]] )
    }
  }
)















