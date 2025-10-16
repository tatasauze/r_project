e4gthead <- function(din) {
  # e4gthead - Gets the common header of din.
  #
  # Args:
  #   din: din of THD format
  #
  # Returns:
  #   A list containing:
  #     H_D: common head dimension
  #     type: model type
  #     m: number of endogenous variables
  #     r: number of exogenous variables
  #     s: seasonal period
  #     n: number of states
  #     np: number of parameters
  #     userflag: flag for user models
  #     userf: user function name
  #     innov: innovation model flags
  #     szpriv: private data sizes

  if (missing(din)) {
    stop("Input 'din' is required.")
  }
  if (length(din) < 17) {
    stop("Input 'din' has an invalid size.")
  }

  H_D <- 28
  type <- din[1]
  m <- din[2]
  r <- din[3]
  s_val <- din[4]
  n_val <- din[5]
  np <- din[6]
  userflag <- din[7]
  innov <- din[24:26]
  szpriv <- din[27:28]

  userf <- NULL
  if (userflag == 1) {
    userf <- intToUtf8(din[8:15])
  } else if (userflag == 2) {
    userf <- c(intToUtf8(din[8:15]), intToUtf8(din[16:23]))
  }

  return(list(
    H_D = H_D,
    type = type,
    m = m,
    r = r,
    s = s_val,
    n = n_val,
    np = np,
    userflag = userflag,
    userf = userf,
    innov = innov,
    szpriv = szpriv
  ))
}