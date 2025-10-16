e4init <- function() {
  # 定義全域變數
  assign("E4OPTION", if (!exists("E4OPTION", envir = .GlobalEnv)) rep(0, 20) else get("E4OPTION", envir = .GlobalEnv), envir = .GlobalEnv)

  # 印出 ASCII 標誌與資訊
  cat("\n")
  cat("                      EEEEEEEEE   444  444\n")
  cat("                     EEEEEEEEEEE  444  444\n")
  cat("                     EEE          44444444\n")          
  cat("                     EEE           4444444\n")
  cat("                     EEEEEEE           444\n")
  cat("                     EEEEEEE           444\n")
  cat("                     EEE\n")
  cat("                     EEE\n")
  cat("                     EEEEEEEEEE\n") 
  cat("                      EEEEEEEE\n")
  cat("\n")
  cat("Toolbox for State Space Estimation of Econometric Models\n")
  cat("                   Version  APR-2013\n")
  cat("\n")
  cat("Web: www.ucm.es/e-4/\n\n")

  # 呼叫設定函數
  sete4opt()

  # 設定錯誤訊息
  ERRORSTR <- c(
    "1.  THETA and DIN do not fit",
    "2.  i inconsistent with THETA (out of range)",
    "3.  Incorrect number of arguments",
    "4.  Badly conditioned covariance matrix",
    "5.  Incorrect model specification",
    "6.  Only one series allowed",
    "7.  Should be more than 1 observation",
    "8.  Model %1d inconsistent",
    "9.  Endogenous variables model should be simple",
    "10. Model not identified",
    "11. Inconsistent input arguments",
    "12. Inconsistent error model",
    "13. User function should be passed as argument in user models",
    "14. Incorrect model",
    "15. Inconsistent system matrix dimension",
    "16. Impossible to compute with missing data",
    "17. Invalid number of lags",
    "18. File not found: %s",
    "19. The equation has no solution",
    "20. SETE4OPT. Unrecognized option %s",
    "21. SETE4OPT. Unrecognized value %s",
    "22. SETE4OPT. Invalid value for %s",
    "23. Run E4INIT before using E4",
    "24. Use ARMA2THD for ARMA models",
    "25. Initial conditions are meaningless",
    "26. Non-stationary system. Initial conditions not compatible with Chandrasekhar",
    "27. E4MIN. No decision variables; check second column of THETA",
    "28. E4LNSRCH. THETA vector is meaningless",
    "29. Multivariate time-varying parameters models are not supported",
    "30. The sample size should be an integer multiple of the seasonal period",
    "31. SETE4OPT. If vcond=De Jong, filter must be Kalman",
    "32. Argument should be scalar",
    "33. E4MIN. Objective function not found",
    "34. For this type of model ARMA2THD or STR2THD should be used",
    "35. Not enough data for using e4preest()",
    "36. User function not defined for the analytic gradient",
    "37. Function only compatible with new versions of matlab (uses cell arrays)"
  )
  assign("ERRORSTR", ERRORSTR, envir = .GlobalEnv)

  # 設定警告訊息
  WARNSTR <- c(
    "1.  Should be one title per series",
    "2.  Invalid number of lags",
    "3.  Invalid %s option",
    "4.  PLOTSERS. A maximum of seven series can be represented in mode 2",
    "5.  RMEDSERS. Invalid group length",
    "6.  LFMODINI. Roots within the circle of radius 1",
    "7.  Approximate computation of information matrix",
    "8.  Information matrix sd+ o d-. Pseudo-inverse computed",
    "9.  E4MIN. Surpassed the maximum number of iterations",
    "10. - - - - - empty - - - - -",
    "11. E4MIN. Hessian reinitialized",
    "12. - - - - - empty - - - - -",
    "13. E4LNSRCH. Precision problem",
    "14. CHOLP. Matrix not square",
    "15. Kalman filter will be used",
    "16. E4MIN. Analytic gradient function not found. Numeric approximation used",
    "17. ARE has no solution because system is non detectable",
    "18. Newton algorithm applied in ARE. NOT convergence after %2d iterations",
    "19. Newton algorithm applied in ARE. Convergence achieved after %2d iter."
  )
  assign("WARNSTR", WARNSTR, envir = .GlobalEnv)

  # 設定顯示訊息
  DISPSTR <- c(
    "\nIteration: %4d | Objf: %8.4f",
    "   x =",
    "   g =",
    "  x0 =",
    "  g0 =",
    "Maximum step length: %8.4f",
    "Step =",
    "E4MIN. Convergence reached in x",
    "E4MIN. Convergence reached in gradient"
  )
  assign("DISPSTR", DISPSTR, envir = .GlobalEnv)
}


sete4opt <- function(...) {
  args <- list(...)
  nargs <- length(args)

  if (!exists("E4OPTION", envir = .GlobalEnv)) {
    e4error(23)
  }

  E4OPTION <- get("E4OPTION", envir = .GlobalEnv)

  if (nargs > 1 && nargs %% 2 != 0) {
    e4error(3)
  }

  e4file <- if (length(E4OPTION) >= 20) E4OPTION[20] else 1

  if (nargs == 0) {
    # Default options
    E4OPTION <- numeric(20)
    E4OPTION[1] <- 1       # Kalman filter
    E4OPTION[2] <- 0       # No scaling
    E4OPTION[3] <- 4       # IDEJ (Inverse De Jong)
    E4OPTION[4] <- 5       # auto
    E4OPTION[5] <- 0       # variance
    E4OPTION[6] <- 1       # BFGS
    E4OPTION[7] <- 0.1     # step
    E4OPTION[8] <- 1e-5    # tolerance
    E4OPTION[9] <- 75      # maxiter
    E4OPTION[10] <- 1      # verbose

    # Internal-use options
    E4OPTION[11] <- 1e-5
    E4OPTION[12] <- 1 - 1e-5
    E4OPTION[13] <- 1e-10
    E4OPTION[14] <- 1e-4
    E4OPTION[15] <- 1e-10
    E4OPTION[16] <- 1e5
    E4OPTION[17] <- 1e-6
    E4OPTION[18] <- if (exists("OCTAVE_VERSION")) 1 else 0
    E4OPTION[19] <- as.numeric(substr(getRversion(), 1, 3))
    E4OPTION[20] <- 1

    assign("E4OPTION", E4OPTION, envir = .GlobalEnv)

    cat("\n")
    return(invisible(E4OPTION))
  }

  if (nargs == 1 && is.character(args[[1]]) && tolower(substr(args[[1]], 1, 3)) == "sho") {
    cat("*********************** Options set by user ********************\n")
    s1 <- if (E4OPTION[1] == 1) "KALMAN" else "CHANDRASEKHAR"
    s2 <- if (E4OPTION[2] == 1) "YES" else "NO"
    s3 <- switch(E4OPTION[3], "LYAPUNOV", "ZERO", "DJONG", "IDEJONG", "UNKNOWN")
    s4 <- switch(E4OPTION[4],
                 "u0 = EXOGENOUS MEAN",
                 "MAXIMUM LIKELIHOOD",
                 "ZERO",
                 "u0 = EXOGENOUS FIRST VALUE (u1)",
                 "AUTOMATIC SELECTION")
    s5 <- if (E4OPTION[5] == 1) "FACTOR" else "VARIANCE"
    s6 <- if (E4OPTION[6] == 1) "BFGS" else "NEWTON"
    s10 <- if (E4OPTION[10] == 1) "YES" else "NO"

    cat(sprintf("Filter. . . . . . . . . . . . . : %s\n", s1))
    cat(sprintf("Scaled B and M matrices . . . . : %s\n", s2))
    cat(sprintf("Initial state vector. . . . . . : %s\n", s4))
    cat(sprintf("Initial covariance of state v.  : %s\n", s3))
    cat(sprintf("Variance or Cholesky factor? .  : %s\n", s5))
    cat(sprintf("Optimization algorithm. . . . . : %s\n", s6))
    cat(sprintf("Maximum step length . . . . . . : %8.6f\n", E4OPTION[7]))
    cat(sprintf("Stop tolerance. . . . . . . . . : %8.6f\n", E4OPTION[8]))
    cat(sprintf("Max. number of iterations . . . : %8i\n", E4OPTION[9]))
    cat(sprintf("Verbose iterations. . . . . . . : %s\n", s10))
    cat("****************************************************************\n\n\n")
    return(invisible(E4OPTION))
  }

  # Modify user-supplied options
  cat("\n\n************** The following options are modified **************\n")

  for (i in seq(1, nargs, by = 2)) {
    optstr <- tolower(args[[i]])
    optval <- args[[i + 1]]

    if (!is.character(optval) && !is.numeric(optval)) e4error(22, optstr)

    if (startsWith(optstr, "fil")) {
      if (startsWith(tolower(optval), "k")) {
        E4OPTION[1] <- 1
      } else if (startsWith(tolower(optval), "c")) {
        if (E4OPTION[3] == 3) e4error(31)
        E4OPTION[1] <- 2
      } else e4error(21, optval)
      cat(sprintf("Filter. . . . . . . . . . . . . : %s\n", toupper(optval)))

    } else if (startsWith(optstr, "sca") || startsWith(optstr, "esc")) {
      if (startsWith(tolower(optval), "y") || startsWith(tolower(optval), "s")) {
        E4OPTION[2] <- 1
      } else if (startsWith(tolower(optval), "n")) {
        E4OPTION[2] <- 0
      } else e4error(21, optval)
      cat(sprintf("Scaled B and M matrices . . . . : %s\n", toupper(optval)))

    } else if (startsWith(optstr, "vco")) {
      val <- substr(tolower(optval), 1, 1)
      E4OPTION[3] <- switch(val,
                            "l" = 1,
                            "z" = 2,
                            "c" = 2,
                            "i" = 4,
                            e4error(21, optval))
      cat(sprintf("Initial covariance of state v.  : %s\n", toupper(optval)))

    } else if (startsWith(optstr, "eco")) {
      val <- tolower(optval)
      E4OPTION[4] <- if (startsWith(val, "a") && nchar(val) < 3) 1 else
                     if (startsWith(val, "m")) 2 else
                     if (startsWith(val, "z") || startsWith(val, "c")) 3 else
                     if (startsWith(val, "i")) 4 else
                     if (startsWith(val, "aut")) 5 else e4error(21, optval)
      cat(sprintf("Initial state vector. . . . . . : %s\n", toupper(optval)))

    } else if (startsWith(optstr, "var")) {
      E4OPTION[5] <- if (startsWith(tolower(optval), "f")) 1 else
                     if (startsWith(tolower(optval), "v")) 0 else e4error(21, optval)
      cat(sprintf("Variance or Cholesky factor? .  : %s\n", toupper(optval)))

    } else if (startsWith(optstr, "alg")) {
      E4OPTION[6] <- if (startsWith(tolower(optval), "b")) 1 else
                     if (startsWith(tolower(optval), "n")) 2 else e4error(21, optval)
      cat(sprintf("Optimization algorithm. . . . . : %s\n", toupper(optval)))

    } else if (startsWith(optstr, "ste") || startsWith(optstr, "pas")) {
      optval <- as.numeric(optval)
      if (optval > 0) {
        E4OPTION[7] <- optval
        cat(sprintf("Maximum step length . . . . . . : %8.6f\n", optval))
      } else e4error(22, optval)

    } else if (startsWith(optstr, "tol")) {
      optval <- as.numeric(optval)
      if (optval > 0) {
        E4OPTION[8] <- optval
        cat(sprintf("Stop tolerance. . . . . . . . . : %8.6f\n", optval))
      } else e4error(22, optval)

    } else if (startsWith(optstr, "max")) {
      optval <- as.integer(optval)
      if (optval > 0) {
        E4OPTION[9] <- optval
        cat(sprintf("Max. number of iterations . . . : %8i\n", optval))
      } else e4error(22, optval)

    } else if (startsWith(optstr, "ver")) {
      if (startsWith(tolower(optval), "y")) {
        E4OPTION[10] <- 1
      } else if (startsWith(tolower(optval), "n")) {
        E4OPTION[10] <- 0
      } else e4error(21, optval)
      cat(sprintf("Verbose iterations. . . . . . . : %s\n", toupper(optval)))

    } else {
      e4error(20, optstr)
    }
  }

  cat("****************************************************************\n\n\n")
  assign("E4OPTION", E4OPTION, envir = .GlobalEnv)
  invisible(E4OPTION)
}
