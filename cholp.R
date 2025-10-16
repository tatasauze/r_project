cholp <- function(H) {
  # cholp - Computes the perturbed Cholesky factor R for a symmetric matrix H.
  #
  # This function uses the Matrix::nearPD function to find the nearest
  # positive definite matrix to H and then computes the Cholesky factor.
  # This is an alternative to the complex perturbation logic in the original
  # MATLAB choldc/cholp functions.
  #
  # Args:
  #   H: A symmetric matrix.
  #
  # Returns:
  #   The Cholesky factor R, such that R'R is a positive definite matrix
  #   close to H.

  # Ensure the Matrix package is available
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 'Matrix' package is required for the cholp function. Please install it.", call. = FALSE)
  }

  # Check if the matrix is already positive definite
  is_pd <- tryCatch({
    chol(H)
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (is_pd) {
    return(chol(H))
  }

  # If not positive definite, find the nearest PD matrix
  pd_matrix_info <- Matrix::nearPD(H, corr = FALSE, keepDiag = TRUE, do2eigen = TRUE, trace = FALSE)

  # Get the resulting positive definite matrix
  pd_matrix <- as.matrix(pd_matrix_info$mat)

  # Compute the Cholesky factor of the new matrix
  R <- chol(pd_matrix)

  return(R)
}