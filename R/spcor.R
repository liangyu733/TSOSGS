#' @title mkR
#'
#' @description
#' A function to generate spatial correlation matrix R for a given retangular field.
#'
#' @param i A integer indicates the plot number corresponding to the row number of spatial correlation matrix R.
#' @param j A integer indicates the plot number corresponding to the column number of spatial correlation matrix R..
#' @param rho A numeric value for spatial correlation parameter ranges from 0 to 1.
#' @param fun 3 spatial correlation functions "power", "exponential", and "Gaussian" are provided.
#'
#' @noRd

spcor <- function(i, j, rho = 0.5, fun = "power") {
  dij <- sqrt((i[1]-j[1])^2+(i[2]-j[2])^2)
  if (fun == "power") {value <- rho^dij}
  if (fun == "exponential") {value <- exp(-dij/rho)}
  if (fun == "Gaussian") {value <- exp(-dij^2/rho^2)}
  return(value)
}
