recpart_fwd <- function(y,x,Mmax){
  if(Mmax<2) {
    warning("Input Mmax must be >= 2; setting to 2")
    Mmax <- 2
  } 
  N <- length(y) 
  n <- ncol(x) 
  B <- init_B(N,Mmax) 
  splits <- vector(mode="list",length=Mmax) # **2
  #---------------------------------------------------
  for(M in 1:Mmax) { 
    lof_best <- Inf
    for(m in 1:M) { 
      for(v in 1:n){ 
        tt <- split_points(x[,v],B[,m])
        for(t in tt) { 
          Bnew <- data.frame(B[,(1:M)[-m]],
                               #Btem1=B[,m]*(x[,v]>t),Btem2=B[,m]*(x[,v]<=t)) 
                               Btem1=B[,m]*H(x[,v]-t),Btem2=B[,m]*H(-(x[,v]-t))) # ** 1
          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(y~.,gdat)
          if(lof < lof_best) { 
            lof_best <- lof
            split_best <- c(m=m,v=v,t=t) # **2
          } 
        } 
      } 
    } 
    cat("best split",split_best,"\n")
    m <- split_best["m"]; v <- split_best["v"]; t <- split_best["t"]
    # B[,M+1] <- B[,m]*(x[,v]<=t) 
    B[,M+1] <- B[,m]*H(-(x[,v]-t)) # **1
    splits[[M+1]] <- rbind(splits[[m]], # **2
                           c(s=-1,v,t)) 
    #B[,m] <- B[,m]*(x[,v]>t)
    B[,m] <- B[,m]*H(x[,v]-t) # **1
    splits[[m]] <- rbind(splits[[m]], # **2
                         c(s=1,v,t))
  } # end loop over M
  return(list(B=B,splits=splits))
}
H <- function(x) { 
  return(as.numeric(x>=0)) 
}
LOF <- function(form,data) {
  ff <- lm(form,data)
  return(sum(residuals(ff)^2))
}
#-------------------------------------
init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}
split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}

set.seed(123); n <- 10
x <- data.frame(x1=rnorm(n),x2=rnorm(n))
y <- rnorm(n)
rp_fwd <- recpart_fwd(y,x,Mmax=9)
rp_fwd$splits
