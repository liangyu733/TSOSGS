#' @title CDmean(v2) criterion
#'
#' @description
#' Compute the CDmean.v2 value for a given training set.
#'
#' @param Kb Kinship matrix of breeding population.
#' @param train An integer vector indicating corresponding allocation of a training population.
#' @param test An integer vector of targeted population.
#' @param R A nt*nt covariance matrix of random error for spatial correlation.
#' @param lambda A variance parameter representing the ratio of error variance to additive genetic variance for the computation of CDmean. Defaulted to 1.
#' @param tolParInv A tolerance parameter of matrix inversion.
#'
#' @export

CDmean.v2 <- function(Kb, train, test, R = NULL, lambda = NULL, tolParInv = 10^-4) {
  if (is.null(lambda)) {lambda <- 1}
  n1 <- length(train)
  n2 <- length(test)
  I <- diag(n1)
  J <- matrix(1, n1, n1)
  if (is.null(R)) {R <- I}
  Rinv <- INV(R, tolParInv)
  M <- Rinv %*% (I - (J %*% Rinv)/(sum(Rinv[])))
  A <- Kb[test, train] %*% INV(M %*% Kb[train, train] + lambda*I, tolParInv) %*% M %*% Kb[train, test]
  B <- Kb[test, test]
  CDmean <- sum(diag(A)/diag(B)) / n2
  return(CDmean)
}
