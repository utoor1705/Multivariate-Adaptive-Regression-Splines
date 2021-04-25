print.mars <- function(object) {
  cat("mars object\n\n")
  cat("Call:\n",deparse(object$call),"\n")
  # or just print(object$call)
  cat("\nCoefficients: \n")
  print(object$coefficients) 
}
