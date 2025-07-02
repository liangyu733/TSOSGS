#' @title Normalized Discounted Cumulative Gain
#'
#' @description
#' An evaluation metric that measures the performance of predicted top k
#' individuals by giving discounted weights to their corresponding true values.
#'
#' @param v A numeric vector representing true values for breeding population.
#' @param vhat A numeric vector representing predicted values for breeding
#'    population.
#' @param k k An integer representing the number of individuals that will be
#'    selected in GS procedure.
#'
#' @export

NDCG <- function(v, vhat, k, gain = "linear") {
  if (length(v) != length(vhat)) {
    stop("Dimension of two vectors are NOT equal.")
  }
  N <- length(v)
  if (k > N) {
    k <- N
    warning("Force k <- length(v); k must <= N.")
  }
  if (gain=="exp") v <- exp(v)
  pi0 <- rank(-vhat,ties.method="random")
  pi <- rank(-v,ties.method="random")

  d <- 1/log2(seq(k) +1)
  DCG <- as.numeric(v[order(pi0)[seq(k)]] %*% d)
  IDCG <- as.numeric(v[order(pi)[seq(k)]] %*% d)
  NDCG <- DCG/IDCG
  return(NDCG)
}
