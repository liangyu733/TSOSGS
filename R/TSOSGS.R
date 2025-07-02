#' @title Training Set Optimization for Spatial Genomic Selection
#'
#' @description
#' This function provides Two-Stage or One-Stage Genetic Search algorithm for
#' Training Set Optimization for Spatial Genomic Selection under an additive
#' spatial GBLUP model.
#' GBLUP model: yt ~ mu*1 + gt + e0, where
#' yt is phenotypic values for a given training population;
#' mu is fixed parameter of population mean;
#' 1 denotes column vector with length nt;
#' gt is genotypic values of corresponding individuals for training set;
#' e0 is random error for a specific field with known spatially correlated
#' experimental units.
#' gt ~ MVN(0, sigma_g^2*Kt),
#' e0 ~ MVN(0, sigma_e^2*R0).
#'
#' @param Kb A nb*nb Kinship matrix of breeding population.
#' @param nt A integer of training set size as desired.
#' @param test A vector indicating target population.
#' @param R A nt*nt covariance matrix of random error for spatial correlation.
#' @param lambda A single value for computing CDmean. Defaulted to 1.
#' @param gia.rf A integer vector of length 2 indicating nrow and ncol of a
#' retangular field for geometrically isomorphic allocation. See function `GIA`.
#' @param method Two-Stage "TS" or One-Stage "OS" method are provided.
#' @param ns Number of solutions for one generation in GA.
#' @param threshold An acceptable threshold of improvement of nearest 100
#' iterations in GA stop criterion.
#' @param nIters Number of iterations that must be reached in GA.
#' @param tolParInv A tolerance parameter of matrix inversion.
#' @param mc.cores Number of used cores. See funtion `mclapply` in package
#' `parallel`.
#' @param verbose Show details in console if TRUE and vice versa.
#'
#' @return This function will return the optimal training population from
#' candidate with the best allocation in field.
#'
#' @examples
#' # Using built-in dataset as demonstration.
#' Kb <- data(tropical.rice)
#' nt <- 50
#' gia.rf <- c(5, 10)
#' R0 <- mkR(c(5, 10), rho = 0.8)
#' # Optimization without spatial correlation.
#' TSOSGS(Kb, nt)
#' # Two-Stage method considering geometrically isomorphic allocation.
#' TSOSGS(Kb, nt, gia.rf = gia.rf)
#' # One-Stage method with specific spatial correlation matrix R.
#' TSOSGS(Kb, nt, R = R0, method = "OS")
#' @importFrom compiler cmpfun
#' @import parallel
#' @export

TSOSGS <- function(Kb, nt, test = NULL, R = NULL, lambda = NULL, gia.rf = NULL,
                   method = "TS", ns = NULL, threshold = 10^-4, nIters = 10^4,
                   tolParInv = 10^-4, mc.cores = 1, verbose = TRUE) {

  ## initial condition
  if (!is.null(gia.rf)) {
    if (!is.vector(gia.rf) ||
        length(gia.rf) != 2 ||
        !is.numeric(gia.rf) ||
        any(gia.rf < 0 | gia.rf %% 1 != 0) ||
        prod(gia.rf) != nt
    ) {
      gia.rf <- NULL
      warning("gia.rf <- c(nrow(Your_field), ncol(Your_field))")
    }
  }

  I <- diag(nt)
  if (is.null(R)) {
    if (is.null(gia.rf)) {
      R <- I
      method <- "OS"
    } else {
      R <- mkR(gia.rf)
    }
  }
  if (ncol(R) == nt && is.pd(R)) {
    colnames(R) <- levels(factor(paste("u", seq(ncol(R)), sep = "")))
    rownames(R) <- colnames(R)
  } else {
    stop("R must be nt*nt and positive definite.")
  }

  if (!method %in% c("TS","OS")) {
    method <- "TS"
    warning("Default Two-Stage method was applied.")
  }

  if (is.null(lambda)) lambda <- 1

  ## candidate population
  if (is.null(test)) {
    nc <- ncol(Kb)
    cand <- test <- seq(nc)       # untargeted
  } else {
    nc <- ncol(Kb) - length(test)
    cand <- seq(nc)[-test]        #   targeted
  }

  stage <- 1
  if (method == "TS") {
    stage <- 0
    nIters <- nIters/2
  }

  ## set ns ##
  ## size of candidate solutions for one generation in GA ##
  if (is.null(ns)) {
    if ( nt >= 10) {
      ns <- nt
    } else {
      ns <- 10
    }
  }

  ## Setup mc.cores for parallel processing
  if (Sys.info()["sysname"] == "Windows") mc.cores <- 1

  while (stage < 2) {

    ## solution matrix
    sol <- mclapply(seq(ns), function(i) {
      GIA(sample(cand, nt), gia.rf)
    }, mc.cores = mc.cores)
    sol <- do.call(cbind, sol)

    ## start GA
    iter <- 0
    stop <- 0
    while(stop == 0) {
      iter <- iter +1
      if (verbose) cat(iter, "...\n", sep = "")

      ## solution score
      score <- c()
      if (iter == 1) {
        max.score <- c()
      }

      ## solution score
      score <- mclapply(seq(ns), function(i) {
        if (stage == 0) {
          return(CDmean.v2(Kb, sol[,i], test, I, lambda, tolParInv))
        }
        if (stage == 1) {
          return(CDmean.v2(Kb, sol[,i], test, R, lambda, tolParInv))
        }
      }, mc.cores = mc.cores)
      score <- do.call(c, score)

      max.score <- c(max.score, max(score)[1])
      if (verbose) cat("max(CDmean):", round(max.score[iter], 5), "\n")

      ## stop criteria
      if (iter >= nIters) {
        if (max(score) - max.score[iter - 100] < threshold) {
          stop <- stop +1
          opt.sol <- which(rank(-score, ties.method = "random") == 1)
        }
      }

      ## elite solution
      elite <- which(rank(-score, ties.method = "random") <= floor(log2(ns)))
      iteropt <- which(rank(-score, ties.method = "random") == 1)
      if (verbose) {
        cat("iter.opt:\n")
        if (is.null(gia.rf)) {
          iteropt <- sol[, iteropt]
        } else {
          iteropt <- matrix(sol[, iteropt], gia.rf[1], gia.rf[2], byrow = TRUE)
        }
        print(iteropt)
      }
      ## Delete solution
      del <- sample(seq(ns)[-elite], floor(log2(ns)),
                    prob = (1/score[-elite]) / sum(1/score[-elite]))

      ## crossover
      cross <- mclapply(del, function(i) {
        chr <- sample(seq(ns)[-del], 2)
        tpl <- sol[, chr[1]]
        cdn <- sol[, chr[2]]
        pos <- sample(nt - 1, 1)
        cdn <- cdn[-seq(pos)]
        upd <- which(cdn %in% tpl[seq(pos)])
        pool <- cand[!cand %in% c(tpl[seq(pos)], cdn[-upd])]
        cdn[upd] <- pool[sample.int(length(pool), length(upd))]
        tpl[-seq(pos)] <- cdn
        return(GIA(tpl, gia.rf))
      }, mc.cores = mc.cores)
      sol[, del] <- do.call(cbind, cross)

      ## mutation
      mut <- mclapply(seq(ns), function(i) {
        sol.new <- sol[, i]
        pos <- sample(nt, 2)
        pool <- c(cand[!cand %in% sol.new], sol[pos, i])
        sol.new[pos] <- sample(pool, 2)
        sol.new <- GIA(sol.new, gia.rf)

        if (i %in% elite) {
          if (stage == 0) {
            old <- CDmean.v2(Kb, sol[, i], test, I, lambda, tolParInv)
            new <- CDmean.v2(Kb, sol.new, test, I, lambda, tolParInv)
          }
          if (stage == 1) {
            old <- CDmean.v2(Kb, sol[, i], test, R, lambda, tolParInv)
            new <- CDmean.v2(Kb, sol.new, test, R, lambda, tolParInv)
          }
          if (new < old) {
            sol.new <- sol[, i]
          }
        }
        return(sol.new)
      }, mc.cores = mc.cores)
      sol <- do.call(cbind, mut)

      if (stop != 0) {
        if (verbose) cat("\n stage-",stage," Genatic Algorithm ended!\n",
                         sep = "")
      }
    } ## END GA

    ## OPT Training set
    if (stage == 0) {
      cand <- sol[, opt.sol]
    } else {
      opt <- sol[, opt.sol]
    }
    stage <- stage +1
  }

  return(opt)
}
TSOSGS <- compiler::cmpfun(TSOSGS)
