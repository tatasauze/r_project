vecss <- function(t, d, nr, nc, sym = FALSE) {
  # vecss - Constructs a nrxnc matrix from a sparse vector where d is
  # a vector with the indices of the theta parameters of M to be estimated.
  # sym shows whether it is a symmetric matrix or not.
  #
  # Args:
  #   t: theta parameters vector
  #   d: indices of the theta parameters
  #   nr: number of rows
  #   nc: number of columns
  #   sym: boolean, TRUE if the matrix is symmetric
  #
  # Returns:
  #   The constructed matrix M.

  t <- as.vector(t)
  M <- matrix(0, nrow = nr, ncol = nc)

  if (length(t) > 0) {
    M[d] <- t
    if (sym) {
      # The parameters 't' at indices 'd' are assumed to fill the lower triangle.
      # This makes the matrix symmetric by copying the lower triangle to the upper triangle.
      M[upper.tri(M)] <- t(M)[upper.tri(M)]
    }
  }
  return(M)
}