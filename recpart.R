# Implementation of recursive partitioning using a binary tree
# to store the partition. 
#-------------------------------------------------------------#
# Binary trees can be implemented as a linked list of nodes. 
# Constructor for the node data structure:
new_node <- function(data,childl=NULL,childr=NULL){
  nn <- list(data=data,childl=childl,childr=childr)
  class(nn) <- "node"
  return(nn)
}
# The data stored in the node are a partition, or region of the 
# covariate space. Constructor for region data structure:
new_region <- function(coords=NULL,x,y){
  if(is.null(coords)) {
    coords <- sapply(x,range)
  }
  out <- list(coords=coords,x=x,y=y)
  class(out) <- "region"
  out
}
#---------------------------------------------------#
# Recursive partitioning function.
recpart <- function(x,y,debug=FALSE){
  init <- new_node(new_region(x=x,y=y))
  tree <- recpart_recursive(init,debug)
  class(tree) <- c("tree",class(tree))
  return(tree)
}
recpart_recursive <- function(node,debug=FALSE) {
  R <- node$data
  # stop recursion if region has a single data point
  if(length(R$y) == 1) { return(NULL) }
  # else find a split that minimizes LOF criterion
  # Initialize
  lof_best  <- Inf
  # Loop over variables and splits
  for(v in 1:ncol(R$x)){ 
    tt <- split_points(R$x[,v]) # Exercise: write split_points()
    for(t in tt) { 
      gdat <- data.frame(y=R$y,x=as.numeric(R$x[,v] <= t))
      lof <- LOF(y~.,gdat) # Exercise: write LOF()
      if(lof < lof_best) { 
        lof_best <- lof
        if(debug) best_split <- c(v=v,t=t)
        childRs <- split(R,xvar=v,spt=t) # Exercises: write split.region()
      }
    }
  } 
  if(debug) {
    cat("best split on variable",best_split["v"], "at", best_split["t"],"\n")
  }
  # Call self on best split
  node$childl <- recpart_recursive(new_node(childRs$Rl),debug)
  node$childr <- recpart_recursive(new_node(childRs$Rr),debug)
  return(node)
}
#------------------------------------------------------------#
split_points <- function(x) {
  # Take a vector of covariate values and return the sorted 
  # unique values. Trim off the max unique value, as this can't
  # be used as a split point.
  x <- sort(unique(x))
  x <- x[-length(x)] # can't split on last value
  return(x)
}
LOF <- function(form,data) {
  ff <- lm(form,data)
  return(sum(residuals(ff)^2))
}
split.region <- function(R,xvar,spt){
  r1_ind <- (R$x[,xvar] <= spt)
  c1 <- c2 <- R$coords
  c1[2,xvar] <- spt; c2[1,xvar] <- spt 
  Rl <- new_region(c1,R$x[r1_ind,,drop=FALSE],R$y[r1_ind])
  Rr <- new_region(c2,R$x[!r1_ind,,drop=FALSE],R$y[!r1_ind])
  return(list(Rl=Rl,Rr=Rr))
}
#---------------------------------------------------#
print.region <- function(R,print.data=FALSE){
  cat("coordinates:\n")
  print(R$coords)
  if(print.data) {
    cat("y:\n")
    print(R$y)
    cat("x:\n")
    print(R$x)
  }
}
# plot partitions on first two variables
plot_regions <- function(x,...)   UseMethod("plot_regions")
plot_regions.tree <- function(tree){
  plot(tree$data$x[,1],tree$data$x[,2],xlab="X1",ylab="X2") 
  plot_regions.node(tree$childl)
  plot_regions.node(tree$childr)
}
plot_regions.node<- function(node) {
  if(is.null(node)) return(NULL)
  x <- node$data$coords[,1]
  y <- node$data$coords[,2]
  lines(c(x[1],x[2],x[2],x[1],x[1]),c(y[1],y[1],y[2],y[2],y[1]),
        col="red")
  plot_regions.node(node$childl)
  plot_regions.node(node$childr)
}
