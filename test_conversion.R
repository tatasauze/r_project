# Test script to compare MATLAB and R implementations of thd2arma
# This script will create test data and compare outputs between MATLAB and R

# Load required libraries
library(R.matlab)

# Source all required R files
source("e4init.R")
source("e4gthead.R") 
source("vecss.R")
source("cholp.R")
source("thd2arm2.R")
source("thd2arma.R")

# Function to create test data that will be used in both MATLAB and R
create_test_data <- function() {
  # --- Define Model Parameters ---
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
  
  # --- Construct `theta` vector ---
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
  
  # --- Construct `din` vector ---
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
  din[(ptr + 1):(ptr + 3)] <- c(1, 3, 4) # Indices for lower triangle of V [1,2,4 would be (1,1), (2,1), (2,2); using 1,3,4 for (1,1), (1,2), (2,2)]
  # Actually corrected: for a 2x2 matrix in column major order: (1,1)=1, (2,1)=2, (1,2)=3, (2,2)=4
  # So for lower triangular we need: (1,1)=1, (2,1)=2, (2,2)=4
  din[(ptr + 1):(ptr + 3)] <- c(1, 2, 4) # Indices for lower triangle of V
  ptr <- ptr + 3 + 1
  
  return(list(theta = theta, din = din))
}

# Function to run the R implementation
run_r_implementation <- function(theta, din) {
  e4option <- rep(0, 28)  # Use default e4option
  result <- thd2arma(theta, din, e4option)
  return(result)
}

# Function to save test data to a MATLAB-compatible format
save_test_data_for_matlab <- function(theta, din, filename = "test_data.mat") {
  # Save test data in MATLAB format
  R.matlab::writeMat(filename, theta = theta, din = din)
  cat("Test data saved to", filename, "\n")
}

# Function to print results for comparison
print_results <- function(r_result, matlab_result = NULL) {
  cat("R Implementation Results:\n")
  cat("F Matrix:\n")
  print(r_result$F)
  cat("A Matrix:\n")
  print(r_result$A)
  cat("G Matrix:\n")
  print(r_result$G)
  cat("V Matrix:\n")
  print(r_result$V)
  
  if (!is.null(matlab_result)) {
    cat("\nMATLAB Implementation Results:\n")
    cat("F Matrix:\n")
    print(matlab_result$F)
    cat("A Matrix:\n")
    print(matlab_result$A)
    cat("G Matrix:\n")
    print(matlab_result$G)
    cat("V Matrix:\n")
    print(matlab_result$V)
    
    # Compare results
    cat("\nComparison Results:\n")
    F_diff <- max(abs(r_result$F - matlab_result$F))
    A_diff <- max(abs(r_result$A - matlab_result$A))
    G_diff <- max(abs(r_result$G - matlab_result$G))
    V_diff <- max(abs(r_result$V - matlab_result$V))
    
    cat("Max difference in F matrix:", F_diff, "\n")
    cat("Max difference in A matrix:", A_diff, "\n")
    cat("Max difference in G matrix:", G_diff, "\n")
    cat("Max difference in V matrix:", V_diff, "\n")
    
    tolerance <- 1e-10
    if (F_diff < tolerance && A_diff < tolerance && G_diff < tolerance && V_diff < tolerance) {
      cat("SUCCESS: R and MATLAB implementations produce equivalent results!\n")
    } else {
      cat("FAILURE: R and MATLAB implementations produce different results!\n")
    }
  }
}

cat("Creating test data...\n")
test_data <- create_test_data()
theta <- test_data$theta
din <- test_data$din

cat("Running R implementation...\n")
r_result <- run_r_implementation(theta, din)

cat("Saving test data for MATLAB...\n")
save_test_data_for_matlab(theta, din, "test_data.mat")

cat("Test data prepared. R results computed.\n")
cat("To complete the comparison, run the MATLAB version with the same inputs and compare the outputs.\n")

print_results(r_result)