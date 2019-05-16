#'Power for covaroances matrices
#'
#'Compute the power of a positive definite symmetric
#'
#'@param x a positive definite symmetric matrix
#'
#'@param n a real number
#'
#'@return a matrix of the same dimensions as \code{x}
#'
#'@keywords internal

"%^%" <- function(x, n) {
    x_pow <- eigen(x)
    x_pow$vectors %*% (x_pow$values^n * t(x_pow$vectors))
}
