plot.mars <- function(object,which=NULL) {
  if(is.null(which)) {
    plot.mars.default(object)
  } else {
    stats:::plot.lm(object,which=which)
  }
}
plot.mars.default <- function(object) {
  yvar <- all.vars(object$formula)[1]
  plot(object$fitted,object$y,xlab="fitted values",ylab=yvar,
       main=deparse(object$call))
}
