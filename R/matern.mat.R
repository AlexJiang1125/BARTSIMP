#' @export
matern.mat <- function(nu, kappa, coords, sigma2e, sigma2x) {
  n <- dim(coords)[1]
  dmat <- dist(coords)
  mcor <- as.matrix(2^(1-nu)*(kappa*dmat)^nu*
                      besselK(dmat*kappa,nu)/gamma(nu))
  diag(mcor) <- 1
  mcov <- sigma2e*diag(n) + sigma2x*mcor
  return(mcov)
}

