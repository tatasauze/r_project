# Source the necessary helper functions
source("e4gthead.R")
source("vecss.R")
source("cholp.R")

thd2arm2 <- function(theta, din, fromgrad = FALSE, e4option = rep(0, 28)) {
  # thd2arm2 - Converts a THD model to factored VARMAX notation.
  #
  # Args:
  #   theta: A vector of parameter values.
  #   din: A vector describing the model's dynamic structure.
  #   fromgrad: A boolean flag (not used in this simplified version).
  #   e4option: A vector for E4 options, e4option[5] is relevant for V matrix.
  #
  # Returns:
  #   A list containing the factored VARMAX matrices: FR, FS, AR, AS, V, G.

  header <- e4gthead(din)
  m <- header$m
  r <- header$r
  type <- header$type
  H_D <- header$H_D

  if (type > 2) {
    stop("Model type > 2 not supported by thd2arm2.")
  }

  p <- din[H_D + 2]
  P <- din[H_D + 3]
  q <- din[H_D + 4]
  Q <- din[H_D + 5]
  g <- din[H_D + 6]

  if (type == 2) {
    p <- p + 1
    P <- P + 1
    q <- q + 1
    Q <- Q + 1
  }

  ptr1 <- 1
  ptr2 <- H_D + 7

  # --- FR ---
  len <- din[ptr2]
  if (len == 0) {
    FR <- matrix(0, nrow = m, ncol = p * m)
  } else {
    indices <- din[(ptr2 + 1):(ptr2 + len)]
    params <- theta[ptr1:(ptr1 + len - 1)]
    FR <- vecss(params, indices, m, p * m)
  }
  ptr1 <- ptr1 + len
  ptr2 <- ptr2 + len + 1

  # --- FS ---
  len <- din[ptr2]
  if (len == 0) {
    FS <- matrix(0, nrow = m, ncol = P * m)
  } else {
    indices <- din[(ptr2 + 1):(ptr2 + len)]
    params <- theta[ptr1:(ptr1 + len - 1)]
    FS <- vecss(params, indices, m, P * m)
  }
  ptr1 <- ptr1 + len
  ptr2 <- ptr2 + len + 1

  # --- AR ---
  len <- din[ptr2]
  if (len == 0) {
    AR <- matrix(0, nrow = m, ncol = q * m)
  } else {
    indices <- din[(ptr2 + 1):(ptr2 + len)]
    params <- theta[ptr1:(ptr1 + len - 1)]
    AR <- vecss(params, indices, m, q * m)
  }
  ptr1 <- ptr1 + len
  ptr2 <- ptr2 + len + 1

  # --- AS ---
  len <- din[ptr2]
  if (len == 0) {
    AS <- matrix(0, nrow = m, ncol = Q * m)
  } else {
    indices <- din[(ptr2 + 1):(ptr2 + len)]
    params <- theta[ptr1:(ptr1 + len - 1)]
    AS <- vecss(params, indices, m, Q * m)
  }
  ptr1 <- ptr1 + len
  ptr2 <- ptr2 + len + 1

  # --- G ---
  len <- din[ptr2]
  if (r > 0 && len > 0) {
    indices <- din[(ptr2 + 1):(ptr2 + len)]
    params <- theta[ptr1:(ptr1 + len - 1)]
    G <- vecss(params, indices, m, (g + 1) * r)
  } else {
    G <- matrix(0, nrow = m, ncol = (g + 1) * r)
  }
  ptr1 <- ptr1 + len
  ptr2 <- ptr2 + len + 1

  # --- V ---
  len <- din[ptr2]
  sym <- if (e4option[5] == 1) FALSE else TRUE
  if (len == 0) {
    V <- matrix(0, nrow = m, ncol = m)
  } else {
    indices <- din[(ptr2 + 1):(ptr2 + len)]
    params <- theta[ptr1:(ptr1 + len - 1)]
    V <- vecss(params, indices, m, m, sym = sym)
  }

  # --- Force V positive definite ---
  if (!fromgrad) {
    if (e4option[5] != 1) {
      # Check if V is positive definite, if not, use cholp to find nearest PD matrix
      is_pd <- tryCatch({
        chol(V)
        TRUE
      }, error = function(e) FALSE)

      if (!is_pd) {
        Vu <- cholp(V)
        V  <- t(Vu) %*% Vu
      }
    } else {
      # This corresponds to the MATLAB V = V*V'
      V <- V %*% t(V)
    }
  }

  return(list(FR = FR, FS = FS, AR = AR, AS = AS, V = V, G = G))
}