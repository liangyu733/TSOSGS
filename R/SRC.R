#' @title Spearman's Rank Correlation
#'
#' @description
#' An evaluation metric describing the linear relationship between the
#' performance of true and predicted top k individuals within breeding
#' population.
#'
#' @param v A numeric vector representing true values for breeding population.
#' @param vhat A numeric vector representing predicted values for breeding
#'    population.
#' @param k An integer representing the number of individuals that will be
#'    selected in GS procedure.
#'
#' @import stats
#' @export

SRC <- function(v, vhat, k) {
  if (length(v) != length(vhat)) {
    stop("Dimension of two vectors are NOT equal.")
  }
  N <- length(v)
  if (k > N) {
    k <- N
    warning("Force k <- length(v); k must <= N.")
  }
  pi0 <- rank(-vhat,ties.method="random")
  pi <- rank(-v,ties.method="random")
  pi0_sort <- sort(pi0)
  pi_order_pi0 <- pi[order(pi0)]
  SRC <- cor(pi0_sort[1:k], pi_order_pi0[1:k])
  return(SRC)
}
