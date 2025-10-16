# Source the necessary helper functions
source("thd2arm2.R")
source("e4gthead.R")

thd2arma <- function(theta, din, e4option = rep(0, 28)) {
  # thd2arma - Converts a THD model to reduced form VARMAX notation.
  #
  # Args:
  #   theta: A vector of parameter values.
  #   din: A vector describing the model's dynamic structure.
  #   e4option: A vector for E4 options.
  #
  # Returns:
  #   A list containing the reduced form VARMAX matrices: F, A, V, G.

  # Get factored matrices from thd2arm2
  factored_model <- thd2arm2(theta, din, e4option = e4option)
  FR_in <- factored_model$FR
  FS_in <- factored_model$FS
  AR_in <- factored_model$AR
  AS_in <- factored_model$AS
  V <- factored_model$V
  G0 <- factored_model$G

  # Get model header information
  header <- e4gthead(din)
  m <- header$m
  r <- header$r
  s <- header$s
  type <- header$type
  H_D <- header$H_D

  if (type > 1) {
    stop("Model type > 1 not supported by thd2arma.")
  }

  p <- din[H_D + 2]
  P <- din[H_D + 3]
  q <- din[H_D + 4]
  Q <- din[H_D + 5]
  g <- din[H_D + 6]

  k <- max(p + P * s, q + Q * s, g)

  # Initialize output matrices
  F_mat <- cbind(diag(m), matrix(0, nrow = m, ncol = k * m))
  A_mat <- cbind(diag(m), matrix(0, nrow = m, ncol = k * m))
  G_mat <- matrix(0, nrow = m, ncol = (k + 1) * r)
  if (r > 0) {
    G_mat[, 1:(r * (g + 1))] <- G0
  }

  # --- Polynomial multiplication for F ---
  # F(B) = FR(B) * FS(B^s) = (I + FR1*B + ...)(I + FS1*B^s + ...)

  # Place regular AR coefficients (FR)
  if (p > 0) {
    F_mat[, (m + 1):(m * (p + 1))] <- FR_in
  }

  # Add seasonal and mixed AR coefficients
  if (P > 0) {
    for (i in 1:P) {
      # Get the i-th seasonal coefficient matrix
      FS_i <- FS_in[, ((i - 1) * m + 1):(i * m), drop = FALSE]

      # Pure seasonal term: F_{i*s} += FS_i
      f_cols_seasonal <- (i * s * m + 1):((i * s + 1) * m)
      F_mat[, f_cols_seasonal] <- F_mat[, f_cols_seasonal] + FS_i

      # Mixed terms: F_{i*s + j} += FR_j * FS_i
      if (p > 0) {
        for (j in 1:p) {
          # Get the j-th regular coefficient matrix
          FR_j <- FR_in[, ((j - 1) * m + 1):(j * m), drop = FALSE]

          f_cols_mixed <- ((i * s + j) * m + 1):((i * s + j + 1) * m)
          F_mat[, f_cols_mixed] <- F_mat[, f_cols_mixed] + FR_j %*% FS_i
        }
      }
    }
  }

  # --- Polynomial multiplication for A ---
  # A(B) = AR(B) * AS(B^s) = (I + AR1*B + ...)(I + AS1*B^s + ...)

  # Place regular MA coefficients (AR)
  if (q > 0) {
    A_mat[, (m + 1):(m * (q + 1))] <- AR_in
  }

  # Add seasonal and mixed MA coefficients
  if (Q > 0) {
    for (i in 1:Q) {
      # Get the i-th seasonal coefficient matrix
      AS_i <- AS_in[, ((i - 1) * m + 1):(i * m), drop = FALSE]

      # Pure seasonal term: A_{i*s} += AS_i
      a_cols_seasonal <- (i * s * m + 1):((i * s + 1) * m)
      A_mat[, a_cols_seasonal] <- A_mat[, a_cols_seasonal] + AS_i

      # Mixed terms: A_{i*s + j} += AR_j * AS_i
      if (q > 0) {
        for (j in 1:q) {
          # Get the j-th regular coefficient matrix
          AR_j <- AR_in[, ((j - 1) * m + 1):(j * m), drop = FALSE]

          a_cols_mixed <- ((i * s + j) * m + 1):((i * s + j + 1) * m)
          A_mat[, a_cols_mixed] <- A_mat[, a_cols_mixed] + AR_j %*% AS_i
        }
      }
    }
  }

  return(list(F = F_mat, A = A_mat, V = V, G = G_mat))
}