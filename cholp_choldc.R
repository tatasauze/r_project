cholp <- function(H, S = NULL) {
  # cholp - Computes the perturbed Cholesky factor R for the symmetric
  # matrix H, such that R'*R = H + mu*I.
  # [R, mu] = cholp(H, S)
  # This function is identical to HESSP, but their tolerances may be adapted to
  # work with other matrices (mainly variance-covariance matrices).
  #
  # Based on MATLAB code by Jaime Terceiro (11/3/97)
  # MODELHESS. Dennis & Schnabel, pp. 315
  
  scaling <- FALSE
  
  if (!is.null(S)) {
    scaling <- TRUE
    if (is.vector(S) || any(dim(S) == 1)) {
      S <- as.vector(S)
    } else {
      S <- diag(S)
    }
  }
  
  dims <- dim(H)
  n <- dims[1]
  m <- dims[2]
  
  if (n != m) {
    warning("Matrix H is not square")
    n <- min(n, m)
    H <- H[1:n, 1:n]
  }
  
  if (scaling) {
    if (length(S) >= n) {
      S <- S[1:n]
      d <- diag(1/S, nrow = n, ncol = n)  # Explicit dimension to avoid vector issues
      H <- d %*% H %*% t(d)
    } else {
      scaling <- FALSE
    }
  }
  
  sqrteps <- sqrt(.Machine$double.eps)
  maxdiag <- max(diag(H))
  mindiag <- min(diag(H))
  maxpos <- max(0.0, maxdiag)
  
  if (mindiag <= sqrteps * maxpos) {
    mu <- 2 * (maxpos - mindiag) * sqrteps - mindiag
    maxdiag <- maxdiag + mu
  } else {
    mu <- 0.0
  }
  
  maxoff <- 0.0
  if (n > 1) {
    maxoff1 <- numeric(n - 1)
    for (i in 1:(n-1)) {
      maxoff1[i] <- max(abs(H[i, (i+1):n]))
    }
    maxoff <- max(maxoff1)
  }
  
  if ((maxoff * (1 + 2*sqrteps)) > maxdiag) {
    mu <- mu + (maxoff - maxdiag) + (2 * sqrteps * maxoff)
    maxdiag <- maxoff * (1 + 2*sqrteps)
  }
  
  if (maxdiag == 0) {
    mu <- 1.0
    maxdiag <- 1.0
  }
  
  if (mu > 0) {
    H <- H + mu * diag(n)
  }
  
  maxoffl <- sqrt(max(maxdiag, maxoff/n))
  
  result1 <- choldc(H, maxoffl)
  L <- result1$L
  maxadd <- result1$maxadd
  
  if (maxadd > 0) {
    maxev <- H[1, 1]
    minev <- H[1, 1]
    
    for (i in 1:n) {
      if (i < n) {
        offrow_right <- sum(abs(H[i, (i+1):n]))
      } else {
        offrow_right <- 0
      }
      
      if (i > 1) {
        offrow_left <- sum(abs(H[1:(i-1), i]))
      } else {
        offrow_left <- 0
      }
      
      offrow <- offrow_right + offrow_left
      maxev <- max(maxev, H[i, i] + offrow)
      minev <- min(minev, H[i, i] - offrow)
    }
    
    sdd <- (maxev - minev) * sqrteps - minev
    sdd <- max(sdd, 0.0)
    mu <- min(maxadd, sdd)
    H <- H + mu * diag(n)
    
    result2 <- choldc(H, 0.0)
    L <- result2$L
    maxadd <- result2$maxadd
  }
  
  if (scaling) {
    d <- diag(S, nrow = n, ncol = n)  # Explicit dimension
    L <- L %*% d
  }
  
  R <- t(L)  # Same return as MATLAB's chol function
  
  return(list(R = R, mu = mu))
}


choldc <- function(H, maxoffl) {
  # choldc - Computes the perturbed Cholesky decomposition LL'=H+D.
  # [L, maxadd] = choldc(H, maxoffl)
  # D is a diagonal non-negative matrix which is computed when it is
  # necessary to grant that a) the elements in the diagonal of L are
  # greater than a tolerance and b) that the elements in the lower triangle
  # are less than maxoffl. L is the perturbed Cholesky factor and maxadd
  # is the bigger element of D.
  #
  # Based on MATLAB code by Jaime Terceiro (11/3/97)
  # Based on CHOLDECOMP. Dennis & Schnabel (1983), pp. 315
  
  n <- nrow(H)
  minl <- sqrt(sqrt(.Machine$double.eps)) * maxoffl
  minl2 <- sqrt(.Machine$double.eps) * maxoffl
  
  if (maxoffl == 0.0) {
    maxoffl <- sqrt(max(abs(diag(H))))
    minl2 <- sqrt(.Machine$double.eps) * maxoffl
  }
  
  maxadd <- 0.0
  L <- matrix(0, nrow = n, ncol = n)
  
  if (n > 1) {
    L[1, 1] <- H[1, 1]
    minljj <- 0.0
    L[2:n, 1] <- H[1, 2:n]
    minljj <- max(max(abs(L[2:n, 1])), minljj)
    minljj <- max(minljj / maxoffl, minl)
    
    if (L[1, 1] > minljj^2) {
      L[1, 1] <- sqrt(L[1, 1])
    } else {
      if (minljj < minl2) {
        minljj <- minl2
      }
      maxadd <- max(maxadd, (minljj^2) - L[1, 1])
      L[1, 1] <- minljj
    }
    L[2:n, 1] <- L[2:n, 1] / L[1, 1]
  }
  
  if (n > 2) {
    for (j in 2:(n-1)) {
      L[j, j] <- H[j, j] - sum(L[j, 1:(j-1)]^2)
      minljj <- 0.0
      L[(j+1):n, j] <- H[j, (j+1):n] - L[(j+1):n, 1:(j-1), drop = FALSE] %*% L[j, 1:(j-1)]
      minljj <- max(max(abs(L[(j+1):n, j])), minljj)
      minljj <- max(minljj / maxoffl, minl)
      
      if (L[j, j] > minljj^2) {
        L[j, j] <- sqrt(L[j, j])
      } else {
        if (minljj < minl2) {
          minljj <- minl2
        }
        maxadd <- max(maxadd, (minljj^2) - L[j, j])
        L[j, j] <- minljj
      }
      L[(j+1):n, j] <- L[(j+1):n, j] / L[j, j]
    }
  }
  
  L[n, n] <- H[n, n] - sum(L[n, 1:(n-1)]^2)
  
  if (L[n, n] > minl^2) {
    L[n, n] <- sqrt(L[n, n])
  } else {
    if (minl < minl2) {
      minl <- minl2
    }
    maxadd <- max(maxadd, (minl^2) - L[n, n])
    L[n, n] <- minl
  }
  
  return(list(L = L, maxadd = maxadd))
}


#===================================
#example
#===================================
if(0>1)
{
 H <- matrix(c(4, 2, 2, 3), nrow = 2)  #should be symmetric
 result <- choldc(H, maxoffl = 1e-8)

 L <- result$L
 maxadd <- result$maxadd

 print(L)
 print(maxadd)
 
 
 H <- matrix(c(4, 2, 2, 3), nrow = 2)
 S <- c(1, 1)  # Scaling vector

 result <- cholp(H, S)

 R <- result$R
 mu <- result$mu

 print(R)
 print(mu)

 # Check: R' * R == H + mu * I (within tolerance)
 print(t(R) %*% R - (H + mu * diag(2)))
}