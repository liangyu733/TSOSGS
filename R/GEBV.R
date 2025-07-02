#' @title GEBV prediction
#'
#' @description
#' Obtain the GEBVs for a given training set by REML or Bayesian estimation.
#'
#' @param Kb A nb*nb Kinship matrix of breeding population.
#' @param train An integer vector indicating corresponding allocation of a training population.
#' @param yt A numeric vector of phenotypic values corresponding to individuals of allocation of a training set.
#' @param R0 The known spatial correlation matrix for a specific field. Consider it independent for units if NULL.
#' @param nIters Number of iterations.
#' @param tolParInv A tolerance parameter of matrix inversion.
#' @param method 2 estimation methods "REML" and "BGS" are provided.
#'
#' @return This function will return a numeric vector represents the correponding GEBVs for entire breeding population.
#'
#' @import sommer
#' @export

GEBV <- function(Kb, train, yt, R0 = NULL, nIters = 5000, tolParInv = 10^-4, method = "BGS") {
  if (!method %in% c("REML", "BGS")) {
    method <- "BGS"
    warning("Default BGS was applied.")
  }
  if (is.null(R0)) R0 <- diag(length(train))
  rownames(R0) <- colnames(R0) <- levels(factor(paste("u", seq(ncol(R0)), sep = "")))
  train0 <- sort(train)
  id <- factor(colnames(Kb)[train])
  pheno <- data.frame(yt = yt, id = id)
  Kt <- Kb[train, train]
  if(method == "BGS") {
    fit <- BGS(yt, Kt, R0, nIters = nIters, tolParInv = tolParInv)
    gthat <- apply(fit[[1]], 3, function(x) mean(x[-seq(nIters*0.9)]))
    muhat <- mean(fit[[2]][-seq(nIters*0.9)])
    gebv <- muhat + Kb[, train] %*% INV(Kt, tolParInv) %*% gthat
  }
  if(method == "REML") {
    fit <- sommer::mmes(fixed = yt~1,
                        random = ~vsm(ism(id), Gu = Kt),
                        rcov = ~vsm(ism(units), Gu = R0),
                        data = pheno,
                        nIters = nIters,
                        tolParInv = tolParInv)
    gebv <- fit$fitted[1] + Kb[, train0] %*% INV(Kb[train0, train0], tolParInv) %*% fit$u
  }
  return(gebv)
}
