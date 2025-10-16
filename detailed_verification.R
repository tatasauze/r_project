# Detailed verification to ensure R implementation matches MATLAB logic exactly
source("e4init.R")
source("e4gthead.R") 
source("vecss.R")
source("cholp.R")
source("thd2arm2.R")
source("thd2arma.R")

# Create a simple test case that allows us to manually verify the algorithm steps
create_simple_test_data <- function() {
  # Simple univariate ARMA(1,1) model - no seasonal components
  # This should be straightforward to verify step by step
  
  m <- 1  # one endogenous variable
  r <- 0  # no exogenous variables
  s <- 4  # seasonal period (not used in this simple case)
  p <- 1  # regular AR order
  P <- 0  # seasonal AR order
  q <- 1  # regular MA order
  Q <- 0  # seasonal MA order
  g <- 0  # exogenous variable lag order
  
  k <- max(p + P * s, q + Q * s, g)  # Should be max(1, 1, 0) = 1
  n_states <- k * m  # 1
  
  # Parameters: just simple AR and MA coefficients
  theta <- c(
    # FR1 (1 param) - AR coefficient
    0.5,
    # FS (0 params) - no seasonal AR
    # AR1 (1 param) - MA coefficient
    0.3,
    # AS (0 params) - no seasonal MA
    # G (0 params) - no exogenous
    # V (1 param) - variance
    1.0
  )
  
  # Build din vector
  din <- numeric(50)
  din[1] <- 1  # type = ARMA
  din[2] <- m
  din[3] <- r
  din[4] <- s
  din[5] <- n_states
  din[6] <- length(theta) # np
  
  H_D <- 28
  din[H_D + 2] <- p  # p = 1
  din[H_D + 3] <- P  # P = 0
  din[H_D + 4] <- q  # q = 1
  din[H_D + 5] <- Q  # Q = 0
  din[H_D + 6] <- g  # g = 0
  
  ptr <- H_D + 7
  
  # FR block (1 param)
  din[ptr] <- 1
  din[(ptr + 1):(ptr + 1)] <- 1  # dense element
  ptr <- ptr + 1 + 1
  
  # FS block (0 params)
  din[ptr] <- 0
  ptr <- ptr + 1
  
  # AR block (1 param)
  din[ptr] <- 1
  din[(ptr + 1):(ptr + 1)] <- 1  # dense element
  ptr <- ptr + 1 + 1
  
  # AS block (0 params)
  din[ptr] <- 0
  ptr <- ptr + 1
  
  # G block (0 params)
  din[ptr] <- 0
  ptr <- ptr + 1
  
  # V block (1 param)
  din[ptr] <- 1
  din[(ptr + 1):(ptr + 1)] <- 1  # single diagonal element
  ptr <- ptr + 1 + 1
  
  return(list(theta = theta, din = din, m = m, r = r, s = s, p = p, P = P, q = q, Q = Q, g = g, k = k))
}

cat("Testing simple univariate ARMA(1,1) case:\n")
simple_data <- create_simple_test_data()
theta_simple <- simple_data$theta
din_simple <- simple_data$din

e4option <- rep(0, 28)
simple_result <- thd2arma(theta_simple, din_simple, e4option)

cat("Simple ARMA(1,1) Results:\n")
cat("F Matrix (should be [1, 0.5]):\n")
print(simple_result$F)
cat("A Matrix (should be [1, 0.3]):\n")
print(simple_result$A)
cat("G Matrix (should be [0]):\n")
print(simple_result$G)
cat("V Matrix (should be [1]):\n")
print(simple_result$V)

# For a simple ARMA(1,1) without seasonal components:
# F(B)y(t) = A(B)e(t) where:
# F(B) = I + F1*B = [1] + [0.5]*B  -> F = [1, 0.5]
# A(B) = I + A1*B = [1] + [0.3]*B  -> A = [1, 0.3]

expected_F <- c(1, 0.5)
expected_A <- c(1, 0.3)

F_match <- all(abs(simple_result$F[1,] - expected_F) < 1e-10)
A_match <- all(abs(simple_result$A[1,] - expected_A) < 1e-10)

cat("\nSimple test validation:\n")
cat("F matrix matches expected [1, 0.5]:", F_match, "\n")
cat("A matrix matches expected [1, 0.3]:", A_match, "\n")

if(F_match && A_match) {
  cat("PASS: Simple ARMA(1,1) test matches expected results\n")
} else {
  cat("FAIL: Simple ARMA(1,1) test does not match expected results\n")
}

# Now let's test the original complex case but verify the logic
cat("\n\nTesting the complex case with manual verification of algorithm:\n")

# Create the same complex test as before
create_complex_test <- function() {
  m <- 2  # two endogenous variables
  r <- 1  # one exogenous variable
  s <- 4  # quarterly seasonality
  p <- 1  # regular AR order
  P <- 0  # seasonal AR order
  q <- 1  # regular MA order
  Q <- 1  # seasonal MA order
  g <- 0  # exogenous variable lag order
  
  k <- max(p + P * s, q + Q * s, g)  # max(1, 1+4, 0) = 5
  n_states <- k * m  # 10
  
  # Parameters
  theta <- c(
    # FR1 (4 params) - AR coefficients
    0.5, 0.1,  # First row of FR1
    0.2, 0.6,  # Second row of FR1
    # FS (is empty for this test - 0 params)
    # AR1 (4 params) - MA coefficients
    0.3, -0.1, # First row of AR1
    0.1, 0.4,  # Second row of AR1
    # AS1 (4 params) - Seasonal MA coefficients
    0.0, 0.2,  # First row of AS1
    0.2, 0.0,  # Second row of AS1
    # G0 (2 params)
    1.1, 1.2,
    # V (3 params for symmetric covariance)
    1.0,   # V[1,1]
    0.5,   # V[2,1] 
    2.0    # V[2,2]
  )
  
  # Build din vector
  din <- numeric(100)
  din[1] <- 1  # type = ARMA
  din[2] <- m
  din[3] <- r
  din[4] <- s
  din[5] <- n_states
  din[6] <- length(theta) # np
  
  H_D <- 28
  din[H_D + 2] <- p  # p
  din[H_D + 3] <- P  # P 
  din[H_D + 4] <- q  # q
  din[H_D + 5] <- Q  # Q
  din[H_D + 6] <- g  # g
  
  ptr <- H_D + 7
  
  # FR block (4 params)
  din[ptr] <- 4
  din[(ptr + 1):(ptr + 4)] <- 1:4
  ptr <- ptr + 4 + 1
  
  # FS block (0 params)
  din[ptr] <- 0
  ptr <- ptr + 1
  
  # AR block (4 params)
  din[ptr] <- 4
  din[(ptr + 1):(ptr + 4)] <- 1:4
  ptr <- ptr + 4 + 1
  
  # AS block (4 params)
  din[ptr] <- 4
  din[(ptr + 1):(ptr + 4)] <- 1:4
  ptr <- ptr + 4 + 1
  
  # G block (2 params)
  din[ptr] <- 2
  din[(ptr + 1):(ptr + 2)] <- 1:2
  ptr <- ptr + 2 + 1
  
  # V block (3 params for symmetric)
  din[ptr] <- 3
  din[(ptr + 1):(ptr + 3)] <- c(1, 2, 4)  # Lower triangular indices
  ptr <- ptr + 3 + 1
  
  return(list(theta = theta, din = din, m = m, r = r, s = s, p = p, P = P, q = q, Q = Q, g = g, k = k))
}

complex_data <- create_complex_test()
theta_complex <- complex_data$theta
din_complex <- complex_data$din

complex_result <- thd2arma(theta_complex, din_complex, e4option)

cat("Complex model results:\n")
cat("F Matrix:\n")
print(complex_result$F)
cat("A Matrix:\n")
print(complex_result$A)

# Extract factored matrices for verification
factored_result <- thd2arm2(theta_complex, din_complex, fromgrad = FALSE, e4option = e4option)

cat("\nFactored components:\n")
cat("FR:", dim(factored_result$FR), "\n")
print(factored_result$FR)
cat("FS:", dim(factored_result$FS), " (should be 2x0)\n")
print(factored_result$FS)
cat("AR:", dim(factored_result$AR), "\n")
print(factored_result$AR)
cat("AS:", dim(factored_result$AS), "\n")
print(factored_result$AS)

# For the complex case: 
# F(B) = FR(B) * FS(B^s) = (I + FR1*B) * (I + 0) = I + FR1*B
# So F should be [I, FR1] which is [eye(2), FR1]
expected_F_complex <- cbind(diag(2), factored_result$FR)
cat("\nExpected F based on algorithm [I, FR1]:\n")
print(expected_F_complex)
cat("Actual F:\n")
print(complex_result$F[,1:4])  # First 4 columns (I and FR1)

F_match_complex <- all(abs(complex_result$F[,1:4] - expected_F_complex) < 1e-10)
cat("F matrix matches expected [I, FR1]:", F_match_complex, "\n")

# For A(B) = AR(B) * AS(B^s) = (I + AR1*B) * (I + AS1*B^4)
# When expanded: I + AR1*B + AS1*B^4 + AR1*AS1*B^5
# This means A should contain: [I, AR1, 0, 0, AS1, AR1%*%AS1]
AR1 <- factored_result$AR
AS1 <- factored_result$AS

# Build expected A matrix based on the algorithm
expected_A_parts <- array(0, c(2, 12))  # 2x12 for m=2, k=5
expected_A_parts[, 1:2] <- diag(2)      # Identity part
expected_A_parts[, 3:4] <- AR1          # Regular MA part (B^1)
expected_A_parts[, 9:10] <- AS1         # Seasonal MA part (B^4, index 4*s*m+1 to (4*s+1)*m = 9-10)
expected_A_parts[, 11:12] <- AR1 %*% AS1 # Mixed part (B^5, index 5*m+1 to 6*m = 11-12)

cat("\nExpected A based on polynomial multiplication:\n")
print(expected_A_parts)
cat("Actual A:\n")
print(complex_result$A)

A_match_complex <- all(abs(complex_result$A - expected_A_parts) < 1e-10)
cat("A matrix matches expected polynomial multiplication result:", A_match_complex, "\n")

if(F_match_complex && A_match_complex) {
  cat("PASS: Complex test matches expected results based on algorithm\n")
} else {
  cat("FAIL: Complex test does not match expected algorithm results\n")
}

cat("\nOverall validation complete.\n")