mars <- function(formula,data,control=NULL,...) {
    cc <- match.call() # save the call
    mf <- model.frame(formula,data)
    y <- model.response(mf)
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    x_names <- colnames(x)
    if(is.null(control)) control <- mars.control()
    fwd <- fwd_selection(y,x,control)
    bwd <- bwd_selection(fwd,control)
    fit <- lm(y~.-1,data=data.frame(y=y,bwd$B)) # notice -1 added
    out <- c(list(call=cc,formula=formula,y=y,B=bwd$B,splits=bwd$splits,
                  x_names=x_names),fit)
    class(out) <- c("mars",class(fit))
    out
  }

fwd_selection <- function(y,x,control=mars.control()){
  N <- length(y) # sample size
  n <- ncol(x) # number of predictors
  B <- matrix(1,nrow=N,ncol=1) 
  splits <- list(data.frame(m=0,v=0,s=NA,t=NA))
  M <- 1
  #---------------------------------------------------
  while(!(M>control$Mmax)) { 
    if(control$trace) cat("M",M,"\n")
    lof_best <- Inf 
    for(m in 1:M) { # choose a basis function to split
      svars <- setdiff(1:n,splits[[m]]$v)
      if(control$trace) cat("M, m, svars",M,m,svars,"\n")
      for(v in svars){ # select a variable to split on
        tt <- split_points(x[,v],B[,m]) 
        for(t in tt) { 
          Bnew <- data.frame(B, 
                             Btem1=B[,m]*h(x[,v],+1,t), 
                             Btem2=B[,m]*h(x[,v],-1,t))
          
          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(y~.,gdat,control) 
          if(lof < lof_best) { 
            lof_best <- lof; m_best <- m; v_best <- v; t_best <- t 
          } 
        } 
      }
    } 
    
    left_split<- rbind(splits[[m_best]], c(m_best,v_best,+1,t_best)) # had order bwds
    right_split<- rbind(splits[[m_best]], c(m_best,v_best,-1,t_best)) # order bwds
    
    splits <- c(splits,list(left_split),list(right_split)) 
    
    B <- cbind(B, 
               B[,m_best]*h(x[,v_best],+1,t_best), 
               B[,m_best]*h(x[,v_best],-1,t_best))
    M <- M+2 # ** 2
    
  } # end while loop over M
  
  colnames(B) <- paste0("B",(0:(ncol(B)-1))) # optional
  return(list(y=y,B=B,splits=splits)) # ** 7
}

bwd_selection <- function(fwd,control) {
  #fwd is a list with elements y, B and splits
  Mmax <- ncol(fwd$B)
  #Jstar <- 1:Mmax
  Jstar <- 2:Mmax
  Kstar <- Jstar
  dat <- data.frame(y=fwd$y,fwd$B)
  lofstar <- LOF(y~.,dat,control)
  for(M in Mmax:2) {
    b <- Inf
    L <- Kstar
    if(control$trace) cat("L:",L,"\n")
    for(m in L){
      K <- setdiff(L,m)
      dat <- data.frame(y=fwd$y,fwd$B[,K])
      lof <- LOF(y~.,dat,control)
      if(control$trace) cat("M:K:lof",M,":",K,":",lof,"\n")
      if(lof < b) {
        b <- lof
        Kstar <- K
      }
      if(lof < lofstar) {
        lofstar <- lof
        Jstar <- K
      }
    }
  }
  Jstar <- c(1,Jstar)
  return(list(y=fwd$y,B=fwd$B[,Jstar],splits=fwd$splits[Jstar]))
}
LOF <- function(form,data,control) {
  ff <- lm(form,data) # should be replaced with lm.fit()
  RSS <- sum(residuals(ff)^2)
  N <- nrow(data)
  M <- length(coef(ff))-1
  Ctilde <- sum(diag(hatvalues(ff))) + control$d*M # ???
  return(RSS * N/(N-Ctilde)^2)
}
h <- function(x,s,t) { 
  return(pmax(0,s*(x-t)))
}
split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}
#------------------------------------------------------------------------
# constructor, validator and helper for class mars.control
new_mars.control <- function(Mmax,d,trace) {
  structure(list(Mmax=Mmax,d=d,trace=trace),class="mars.control")
}
validate_mars.control <- function(control) {
  if(control$Mmax < 2) {
    warning("Mmax must be >= 2; setting to 2")
    control$Mmax <- 2
  }
  control
}
mars.control <- function(Mmax=2,d=3,trace=FALSE) {
  control <- new_mars.control(Mmax,d,trace)
  validate_mars.control(control)
}
