anova.mars <- function(object) {
  cat("ANOVA decomposition, basis functions of one or two variables.\n")
  xn <- object$x_names
  bn <- names(object$coefficients)
  cc <- object$coefficients
  ss <- object$splits
  nhinge <- sapply(ss,function(x) nrow(x) - 1)
  i1 <- which(nhinge==1)
  i2 <- which(nhinge==2)
  imore <- setdiff(2:length(nhinge),union(i1,i2))
  cat("Basis functions of one variable\n")
  for(i in i1) {
    cat(paste0(bn[i],":\n"))
    cat(paste0("coefficient=",cc[i],",\n"))
    cat(paste0(" variable ",xn[ss[[i]][2,"v"]],";"))
    cat(paste0(" sign ",ss[[i]][2,"s"],";"))
    cat(paste0(" split at value ",ss[[i]][2,"t"],"\n"))
  }
  if(length(i2)>0) {
    cat("\nBasis functions of two variables\n")
    for(i in i2) {
      cat(paste0(bn[i],":\n"))
      cat(paste0("coefficient=",cc[i],",\n"))
      for(j in 2:nrow(ss[[i]])) {
        cat(paste0("  Component ",j-1,":  variable ",xn[ss[[i]][j,"v"]],";"))
        cat(paste0(" sign ",ss[[i]][j,"s"],";"))
        cat(paste0(" split at value ",ss[[i]][j,"t"],"\n"))
      }
    } 
  } else {
    cat("\nNo basis functions of 2 variables\n\n")
  }
  if(length(imore)>0) {
    cat(paste(length(imore),"basis functions of 2 or more variables\n"))
  } else {
    cat("No basis functions of 2 or more variables\n")
  }
}
