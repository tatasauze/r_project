# Validate the R implementation by checking for logical correctness
# Since we can't run MATLAB directly, we'll verify the R implementation follows the same logic
# as outlined in the MATLAB code and produces mathematically sensible results

source("e4init.R")
source("e4gthead.R") 
source("vecss.R")
source("cholp.R")
source("thd2arm2.R")
source("thd2arma.R")

# Function to create test data that will be used in testing
create_test_data <- function() {
  m <- 2  # two endogenous variables
  r <- 1  # one exogenous variable
  s <- 4  # quarterly seasonality
  p <- 1  # regular AR order
  P <- 0  # seasonal AR order
  q <- 1  # regular MA order
  Q <- 1  # seasonal MA order
  g <- 0  # exogenous variable lag order
  
  k <- max(p + P * s, q + Q * s, g)
  n_states <- k * m
  
  # Parameters for the model
  theta <- c(
    # FR1 (4 params) - AR coefficients
    0.5, 0.1,  # First row of FR1
    0.2, 0.6,  # Second row of FR1
    # FS (is empty for this test)
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

# Run the test
cat("Validating R implementation of thd2arma...\n")

test_data <- create_test_data()
theta <- test_data$theta
din <- test_data$din
m <- test_data$m
r <- test_data$r
s <- test_data$s
p <- test_data$p
P <- test_data$P
q <- test_data$q
Q <- test_data$Q
g <- test_data$g
k <- test_data$k

e4option <- rep(0, 28)
result <- thd2arma(theta, din, e4option)

cat("R Implementation Results:\n")
cat("F Matrix (dimensions:", dim(result$F)[1], "x", dim(result$F)[2], "):\n")
print(result$F)
cat("\nA Matrix (dimensions:", dim(result$A)[1], "x", dim(result$A)[2], "):\n")
print(result$A)
cat("\nG Matrix (dimensions:", dim(result$G)[1], "x", dim(result$G)[2], "):\n")
print(result$G)
cat("\nV Matrix (dimensions:", dim(result$V)[1], "x", dim(result$V)[2], "):\n")
print(result$V)

# Validate matrix dimensions
expected_F_cols <- m * (1 + k)  # m columns for identity + k*m columns for lags
expected_A_cols <- m * (1 + k)
expected_G_cols <- r * (1 + k) 
expected_G_rows <- m

cat("\nValidation checks:\n")
cat("Expected F matrix dimensions:", m, "x", expected_F_cols, "\n")
cat("Actual F matrix dimensions:", dim(result$F)[1], "x", dim(result$F)[2], "\n")
cat("F matrix dimension check:", if(dim(result$F)[1] == m && dim(result$F)[2] == expected_F_cols) "PASS" else "FAIL", "\n")

cat("\nExpected A matrix dimensions:", m, "x", expected_A_cols, "\n")
cat("Actual A matrix dimensions:", dim(result$A)[1], "x", dim(result$A)[2], "\n")
cat("A matrix dimension check:", if(dim(result$A)[1] == m && dim(result$A)[2] == expected_A_cols) "PASS" else "FAIL", "\n")

cat("\nExpected G matrix dimensions:", expected_G_rows, "x", expected_G_cols, "\n")
cat("Actual G matrix dimensions:", dim(result$G)[1], "x", dim(result$G)[2], "\n")
cat("G matrix dimension check:", if(dim(result$G)[1] == expected_G_rows && dim(result$G)[2] == expected_G_cols) "PASS" else "FAIL", "\n")

cat("\nV matrix is symmetric check:", if(all(result$V == t(result$V))) "PASS" else "FAIL", "\n")

# Check that F and A start with identity matrices
F_identity_check <- all(diag(result$F[, 1:m, drop=FALSE]) == 1) && all(result$F[upper.tri(result$F[, 1:m, drop=FALSE])] == 0) && all(result$F[lower.tri(result$F[, 1:m, drop=FALSE])] == 0)
A_identity_check <- all(diag(result$A[, 1:m, drop=FALSE]) == 1) && all(result$A[upper.tri(result$A[, 1:m, drop=FALSE])] == 0) && all(result$A[lower.tri(result$A[, 1:m, drop=FALSE])] == 0)

cat("F starts with identity matrix:", if(F_identity_check) "PASS" else "FAIL", "\n")
cat("A starts with identity matrix:", if(A_identity_check) "PASS" else "FAIL", "\n")

# Additional validation based on the expected algorithm behavior

# For the first test, let's run thd2arm2 separately to see the intermediate results
factored_result <- thd2arm2(theta, din, fromgrad = FALSE, e4option = e4option)
cat("\nIntermediate results from thd2arm2:\n")
cat("FR dimensions:", dim(factored_result$FR), "\n")
cat("FS dimensions:", dim(factored_result$FS), "\n")
cat("AR dimensions:", dim(factored_result$AR), "\n")
cat("AS dimensions:", dim(factored_result$AS), "\n")
cat("V is positive definite check:", if(min(eigen(result$V)$values) > 0) "PASS" else "FAIL (not PD)", "\n")

cat("\nValidation complete.\n")

# Run the same test with the existing test script to verify consistency
cat("\nRunning existing test script for consistency check:\n")
source("test_thd2arma.R")