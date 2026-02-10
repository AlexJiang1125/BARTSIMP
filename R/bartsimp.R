#' @references
#' Chipman, H., George, E., McCulloch, R. (2010).
#' BART: Bayesian Additive Regression Trees.
#'
#' The interface design is inspired by the BART R package.
#' @export
bartsimp=function(
  x.train, y.train,s1,s2, x.test=matrix(0.0,0,0),
  size = rep(1,length(s1)),
  weighted = FALSE,
  sparse=FALSE, theta=0, omega=1,
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
  seed = 99L, iftest=FALSE, isexact=FALSE, nwarmup=500,
  usecoords = TRUE,
  sigest_new = 1, range_new = 5, sig_m_new = 1,
  rho_0 = 2.4, sigmam_0 = 0.55,
  alpha_1 = 0.05, alpha_2 = 0.05,
  doBART = FALSE
)
{
  #--------------------------------------------------
  #data
  colvars <- colnames(x.train)
  n = length(y.train)
  #print(numcut)
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

  ##if(p>1 & length(numcut)==1) numcut=rep(numcut, p)

  y.train = y.train-fmean
  #--------------------------------------------------
  #set nkeeps for thinning
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
  #--------------------------------------------------
  #prior
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
    lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
  } else {
    sigest=sqrt(lambda)
  }

  if(is.na(sigmaf)) {
    tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
  } else {
    tau = sigmaf/sqrt(ntree)
  }
  ## 20220918 update
  #--- set spatial hyperpars
  x_train <- as.data.frame(t(x.train))
  colnames(x_train) <- colvars
  # hyperpar <- tunehyperpriors(x_train = x_train,
  #                 y_train = y.train,
  #                 s1 = s1,
  #                 s2 = s2)
  sigest_new <- sigest_new#1.1#hyperpar$sigest_new
  range_new <- range_new#hyperpar$range_new
  sig_m_new <- sig_m_new#hyperpar$sig_m_new
  rho_0 <- rho_0#2.4#hyperpar$rho_0
  sigmam_0 <- sigmam_0#0.55#hyperpar$sigmam_0
  alpha_1 <- alpha_1
  alpha_2 <- alpha_2
  # dist_vec <- as.vector(dist(s_mat))
  # rho_0 <- quantile(dist_vec, 0.5)/2.5
  # sigmam_0 <- sigest * 2.5
  # ## fit INLA model



  #--- set unique
  x.unique <- x.train[,cumsum(size)]
  #print(dim(x.unique))
  n.unique <- length(size)

  #--- check if we need to load from a previous tree object
  if (is.na(treeprev)) {
    tree.update <- FALSE
    treeprev_last <- NA
  } else {
    print("use previous updates")
    tree.update <- TRUE
    # only need the last one!
    treeprev_last <- treeprev#[[length(treeprev)]]
  }

  #--------------------------------------------------
  ptm <- proc.time()
  #call
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
  ##res$nkeeptreedraws=nkeeptreedraws
  res$varcount.mean <- apply(res$varcount, 2, mean)
  res$varprob.mean <- apply(res$varprob, 2, mean)
  res$rm.const <- rm.const

  res$sigmams <- res$sigmams[(nwarmup+1):(nwarmup+ndpost)]
  res$kappas <- res$kappas[(nwarmup+1):(nwarmup+ndpost)]
  res$sigma <- res$sigma[(nwarmup+1):(nwarmup+ndpost)]

  attr(res, 'class') <- 'wbart'
  return(res)
}


# Backward-compatibility wrapper
#' @rdname wbart_spatial
#' @export
mywbart <- function(...) {
  .Deprecated("bartsimp")
  bartsimp(...)
}
