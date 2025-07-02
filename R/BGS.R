#' @title Bayesian Gibbs Sampling method
#'
#' @description
#' A Bayesian Gibbs Sampling method to obtain GEBVs for GBLUP model.
#' BGS for additive spatial GBLUP model:
#' yt ~ mu*1 + gt + e0, where
#' yt is phenotypic values for a given training population;
#' mu is fixed parameter of population mean;
#' 1 denotes column vector with length nt;
#' gt is genotypic values of corresponding individuals for training set;                       # gt ~ MVN(0, sigma_g^2*Kt)
#' e0 is random error for a specific field with known spatially correlated experimental units. # e0 ~ MVN(0, sigma_e^2*R0)
#'
#' @param yt A numeric vector of phenotypic values corresponding to individuals of allocation of a training set.
#' @param Kt The corresponding Kinship matrix of the training set.
#' @param R0 The known spatial correlation matrix for a specific field. Consider it independent for units if NULL.
#' @param mu.ini Initial values for the constant term.
#' @param ga.ini Initial values for the additive genotypic values.
#' @param vE.ini Initial values for the error variance.
#' @param va.ini Initial values for the additive variance.
#' @param nu_star The degree of freedom of the variance components for the prior distribution of a scaled inverse Chi square.
#' @param nIters Number of iterations.
#' @param tolParInv A tolerance parameter of matrix inversion.
#' @param verbose Show details in console if TRUE and vice versa.
#'
#' @return This function will return a list containing estimates of ga, mu, Va, and VE for each iteration.
#'
#' @import MASS stats
#' @importFrom compiler cmpfun
#' @export

BGS <- function(yt, Kt, R0 = NULL, mu.ini = NULL, ga.ini = 0, vE.ini = 1, vA.ini = 0.5, nu_star = 5,
                nIters = 5000, tolParInv = 10^-4, verbose = TRUE){
  n <- length(yt)
  n_fix <- 1
  Kinv <- INV(Kt, tolParInv)

  if (is.null(R0)) {
    R0 <- diag(n)
  }
  Rinv <- INV(R0, tolParInv)

  if (is.null(mu.ini)) {
    mu.ini <- mean(yt, na.rm = T)
  }
  vA <- vA.ini
  vE <- vE.ini
  lamA <- vE / vA

  X <- rep(1, n)
  C11 <- t(X) %*% Rinv %*% X
  C12 <- t(X) %*% Rinv
  C21 <- Rinv %*% X
  C22 <- Rinv + Kinv * lamA

  r1 <- t(X) %*% Rinv %*% yt
  r2 <- Rinv %*% yt
  g1 <- rep(mu.ini, n_fix)
  g2 <- rep(ga.ini, n)

  ga_dat <- matrix(NA, n, iter)
  mu_dat <- rep(NA, iter)
  vA_dat <- rep(NA, iter)
  vE_dat <- rep(NA, iter)

  for (i in 0:nIters) {
    C11inv <- INV(C11, tolParInv)
    Sigma_g1 <- vE*C11inv
    g1_star <- C11inv %*% (r1 - C12 %*% g2)
    g1 <- MASS::mvrnorm(n = 1, mu = g1_star, Sigma = Sigma_g1)

    C22inv <- INV(C22, tolParInv)
    Sigma_g2 <- vE*C22inv
    g2_star <- C22inv %*% (r2 - C21 %*% g1)
    g2 <- MASS::mvrnorm(n = 1, mu = g2_star, Sigma = Sigma_g2)

    e0 <- yt-g1-g2

    S_star <- 0.5*var(yt)

    vE <- (t(e0) %*% e0 + S_star*nu_star) / rchisq(1, n+nu_star)
    vA <- (t(g2) %*% Kinv %*% g2 + S_star*nu_star) / rchisq(1, n+nu_star)
    vE <- as.numeric(vE)
    vA <- as.numeric(vA)
    lamA <- vE/vA
    C22 <- Rinv + Kinv*lamA

    ga_dat[, i] <- g2
    mu_dat[i] <- g1
    vA_dat[i] <- vA
    vE_dat[i] <- vE
  }
  result <- list(ga = ga_dat, mu = mu_dat, vA = vA_dat, vE = vE_dat)
  return(result)
}
BGS <- compiler::cmpfun(BGS)
