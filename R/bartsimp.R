#' Spatial BART Model
#'
#' Fits a Bayesian Additive Regression Trees (BART) model with
#' spatial variance components.
#'
#' This function is a simplified interface to the underlying C++ sampler.
#' The design is inspired by the original BART R package.
#'
#' @param x.train Predictor matrix for training data (n × p).
#' @param y.train Numeric response vector of length n.
#' @param s1 Spatial coordinate vector (length n).
#' @param s2 Spatial coordinate vector (length n).
#' @param x.test Optional predictor matrix for test data.
#'
#' @param size Cluster sizes.
#'
#' @param ntree Number of trees in the ensemble.
#' @param numcut Number of cut points per predictor.
#'
#' @param ndpost Number of posterior draws to keep.
#' @param nskip Number of burn-in iterations.
#' @param keepevery Thinning interval.
#'
#' @param sigest Initial estimate of residual standard deviation.
#' @param sigdf Degrees of freedom for sigma prior.
#' @param sigquant Quantile used to calibrate sigma prior.
#'
#' @param k Tree prior scale parameter.
#' @param power Tree depth prior parameter.
#' @param base Tree depth prior parameter.
#'
#' @param nwarmup Number of warmup iterations.
#' @param seed Random seed.
#'
#' @param usecoords Logical; whether to use spatial coordinates.
#'
#' @return
#' An object containing:
#'
#' \describe{
#'   \item{yhat.train}{Posterior samples for training data}
#'   \item{yhat.test}{Posterior samples for test data}
#'   \item{yhat.train.mean}{Posterior mean predictions for training data}
#'   \item{yhat.test.mean}{Posterior mean predictions for test data}
#'   \item{sigma}{Posterior draws of the observation-level residual standard deviation}
#'   \item{kappas}{Posterior draws of the Matérn inverse range parameter
#'   \eqn{\kappa = \sqrt{8\nu}/\rho}}
#'   \item{sigmams}{Posterior draws of the spatial random-effect standard deviation}
#'   \item{varcount}{Variable inclusion counts}
#'   \item{varprob}{Variable inclusion probabilities}
#'   \item{proc.time}{Execution time}
#' }
#'
#' }
#'
#' @references
#' Jiang, A. Z., & Wakefield, J. (2025).
#' BARTSIMP: Flexible spatial covariate modeling and prediction using
#' Bayesian Additive Regression Trees.
#' *Spatial and Spatio-temporal Epidemiology*, 55, 100757.
#' https://doi.org/10.1016/j.sste.2025.100757
#'
#' The interface design is inspired by the BART R package.
#' @export
bartsimp=function(
    x.train, y.train,s1,s2, x.test=matrix(0.0,0,0),
    size = rep(1,length(s1)),
    theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0.0,0,0), usequants=FALSE,
    cont=FALSE, rm.const=TRUE,
    sigest=0.1, sigdf=3, sigquant=.90,
    k=2.0, power=2.0, base=.95,
    sigmaf=NA, lambda=NA,
    treeprev = NA,
    fmean=mean(y.train),
    w=rep(1,length(y.train)),
    ntree=3L, numcut=100L,
    ndpost=500L, nskip=1L, keepevery=1L,
    nkeeptrain=ndpost, nkeeptest=ndpost,
    nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
    printevery=100L, transposed=FALSE,
    seed = 99L, nwarmup=500,
    usecoords = TRUE,
    sigest_new = 1, range_new = 5, sig_m_new = 1,
    rho_0 = 2.4, sigmam_0 = 0.55,
    alpha_1 = 0.05, alpha_2 = 0.05
)
{
  # --------------------------------------------------
  # DEBUG MODE CONTROL
  debug_mode <- isTRUE(getOption("BARTSIMP.debug", FALSE))
  iftest  <- debug_mode
  isexact <- debug_mode
  doBART  <- debug_mode

  # --------------------------------------------------
  # INTERNAL FLAGS (not user-facing)
  weighted <- FALSE
  sparse <- FALSE

  # --------------------------------------------------
  #data
  colvars <- colnames(x.train)
  n = length(y.train)
  set.seed(seed)

  if(!transposed) {
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    if(length(x.test)>0) {
      x.test = bartModelMatrix(x.test)
      x.test = t(x.test[ , temp$rm.const])
    }
    rm.const <- temp$rm.const
    grp <- temp$grp
    rm(temp)
  } else {
    rm.const <- NULL
    grp <- NULL
  }

  if(n!=ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')

  p = nrow(x.train)
  np = ncol(x.test)
  if(length(rho)==0) rho=p
  if(length(rm.const)==0) rm.const <- 1:p
  if(length(grp)==0) grp <- 1:p

  y.train = y.train-fmean

  if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
    nkeeptrain=ndpost
    cat('*****nkeeptrain set to ndpost\n')
  }
  if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
    nkeeptest=ndpost
    cat('*****nkeeptest set to ndpost\n')
  }
  if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
    nkeeptestmean=ndpost
    cat('*****nkeeptestmean set to ndpost\n')
  }
  if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
    nkeeptreedraws=ndpost
    cat('*****nkeeptreedraws set to ndpost\n')
  }

  nu=sigdf
  if(is.na(lambda)) {
    if(is.na(sigest)) {
      if(p < n) {
        df = data.frame(t(x.train),y.train)
        lmf = lm(y.train~.,df)
        sigest = summary(lmf)$sigma
      } else {
        sigest = sd(y.train)
      }
    }
    qchi = qchisq(1.0-sigquant,nu)
    lambda = (sigest*sigest*qchi)/nu
  } else {
    sigest=sqrt(lambda)
  }

  if(is.na(sigmaf)) {
    tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
  } else {
    tau = sigmaf/sqrt(ntree)
  }

  x_train <- as.data.frame(t(x.train))
  colnames(x_train) <- colvars

  sigest_new <- sigest_new
  range_new <- range_new
  sig_m_new <- sig_m_new
  rho_0 <- rho_0
  sigmam_0 <- sigmam_0
  alpha_1 <- alpha_1
  alpha_2 <- alpha_2

  x.unique <- x.train[,cumsum(size)]
  n.unique <- length(size)

  if (is.na(treeprev)) {
    tree.update <- FALSE
    treeprev_last <- NA
  } else {
    print("use previous updates")
    tree.update <- TRUE
    treeprev_last <- treeprev
  }

  ptm <- proc.time()

  res <- invisible(
    capture.output({
      tmp <- .Call("cwbart",
                   n, p, n.unique, np,
                   x.train, x.unique, y.train, x.test,
                   s1, s2, weighted, size,
                   ntree, numcut,
                   ndpost*keepevery, nskip,
                   power, base, tau, nu, lambda,
                   sigest_new, rho_0, sigmam_0,
                   alpha_1, alpha_2,
                   range_new, sig_m_new, sigest_new,
                   w, sparse, theta, omega, grp,
                   a, b, rho, augment,
                   nkeeptrain, nkeeptest, nkeeptestmean,
                   nkeeptreedraws, printevery,
                   xinfo, iftest, isexact, nwarmup,
                   usecoords, tree.update, treeprev_last,
                   doBART
      )
    })
  )

  res <- tmp

  res$proc.time <- proc.time()-ptm

  res$mu = fmean
  res$yhat.train.mean = res$yhat.train.mean+fmean
  res$yhat.train = res$yhat.train+fmean
  res$yhat.test.mean = res$yhat.test.mean+fmean
  res$yhat.test = res$yhat.test+fmean
  if(nkeeptreedraws>0)
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
  dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
  dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
  res$varcount.mean <- apply(res$varcount, 2, mean)
  res$varprob.mean <- apply(res$varprob, 2, mean)
  res$rm.const <- rm.const

  res$sigmams <- res$sigmams[(nwarmup+1):(nwarmup+ndpost)]
  res$kappas <- res$kappas[(nwarmup+1):(nwarmup+ndpost)]
  res$sigma <- res$sigma[(nwarmup+1):(nwarmup+ndpost)]

  attr(res, 'class') <- 'wbart'
  return(res)
}
