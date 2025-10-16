# Test script for the thd2arma R conversion

# Source the main function, which in turn sources all dependencies
source("thd2arma.R")

# --- 1. Define Model Parameters ---
# We will create a simple bivariate VARMA(1,1) model with a seasonal MA component (s=4, Q=1)
# m = 2 (two endogenous variables)
# r = 1 (one exogenous variable)
# s = 4 (quarterly seasonality)
# p = 1 (regular AR order)
# P = 0 (seasonal AR order)
# q = 1 (regular MA order)
# Q = 1 (seasonal MA order)
# g = 0 (exogenous variable lag order)

m <- 2
r <- 1
s <- 4
p <- 1
P <- 0
q <- 1
Q <- 1
g <- 0

# The maximum lag k is the maximum of the AR and MA polynomial degrees.
# AR degree = p + P*s
# MA degree = q + Q*s
# k = max(p + P*s, q + Q*s, g)
k <- max(p + P * s, q + Q * s, g)
n_states <- k * m

# --- 2. Construct `theta` vector ---
# This vector contains all model parameters in a flat structure.
# Order: FR, FS, AR, AS, G, V
# FR (p=1): 4 parameters for a 2x2 matrix
# FS (P=0): 0 parameters
# AR (q=1): 4 parameters for a 2x2 matrix
# AS (Q=1): 4 parameters for a 2x2 seasonal MA matrix
# G (g=0): 2 parameters for a 2x1 G0 matrix
# V: 3 parameters for a symmetric 2x2 covariance matrix (V11, V21, V22)

theta <- c(
  # FR1 (4 params)
  0.5, 0.1, 0.2, 0.6,
  # FS (is empty)
  # AR1 (4 params)
  0.3, -0.1, 0.1, 0.4,
  # AS1 (4 params)
  0.0, 0.2, 0.2, 0.0,
  # G0 (2 params)
  1.1, 1.2,
  # V (3 params for symmetric)
  1.0, 0.5, 2.0
)

# --- 3. Construct `din` vector ---
# This vector defines the structure of the model.

# Header part (first 28 elements)
din <- numeric(100) # Initialize with zeros
din[1] <- 1  # type = ARMA
din[2] <- m
din[3] <- r
din[4] <- s
din[5] <- n_states
din[6] <- length(theta) # np (number of parameters)

# Body part (parameter structure information)
H_D <- 28
din[H_D + 2] <- p
din[H_D + 3] <- P
din[H_D + 4] <- q
din[H_D + 5] <- Q
din[H_D + 6] <- g

ptr <- H_D + 7

# FR block
din[ptr] <- 4 # 4 params in FR
din[(ptr + 1):(ptr + 4)] <- 1:4 # dense matrix
ptr <- ptr + 4 + 1

# FS block
din[ptr] <- 0 # 0 params in FS
ptr <- ptr + 1

# AR block
din[ptr] <- 4 # 4 params in AR
din[(ptr + 1):(ptr + 4)] <- 1:4 # dense matrix
ptr <- ptr + 4 + 1

# AS block
din[ptr] <- 4 # 4 params in AS
din[(ptr + 1):(ptr + 4)] <- 1:4 # dense matrix
ptr <- ptr + 4 + 1

# G block
din[ptr] <- 2 # 2 params in G
din[(ptr + 1):(ptr + 2)] <- 1:2 # dense matrix
ptr <- ptr + 2 + 1

# V block
din[ptr] <- 3 # 3 params in V for symmetric
din[(ptr + 1):(ptr + 3)] <- c(1, 2, 4) # Indices for lower triangle of V
ptr <- ptr + 3 + 1

# --- 4. Run the conversion ---
# Use default e4option (e4option[5] = 0, so V is treated as symmetric)
e4option <- rep(0, 28)
result <- thd2arma(theta, din, e4option)

# --- 5. Print the results ---
cat("--- Converted VARMAX Matrices ---\n\n")

cat("F Matrix (AR part):\n")
print(result$F)
cat("\nDimension of F: ", dim(result$F), "\n\n")

cat("A Matrix (MA part):\n")
print(result$A)
cat("\nDimension of A: ", dim(result$A), "\n\n")

cat("G Matrix (Exogenous part):\n")
print(result$G)
cat("\nDimension of G: ", dim(result$G), "\n\n")

cat("V Matrix (Covariance):\n")
print(result$V)
cat("\nDimension of V: ", dim(result$V), "\n")

cat("\nTest script finished.\n")