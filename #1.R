####

f1<-function(x){
  log(x)/(1+x)
}

##Bisec function

myBisec <- function(g,a,b,e=10^(-8)){
  maxiter<-1000
  err<-1
  niter<-0
  x0<-(a+b)/2
  
  while ( niter<=maxiter && err >= e){
    if (genD(g,a)$D[1]*genD(g,x0)$D[1] <=0)   {b<-x0}
    else   {a<-x0}
    
    oldx0<-x0
    x0<-(a+b)/2
    
    err<-abs(oldx0-x0)
    niter<-niter+1
  }
  return(x0)
}


##Newton function

myNewton <- function(g,x,e=10^(-8)){
  maxiter<-1000
  err<-1
  niter<-0
  x0<-x
  
  while ( niter<=maxiter && err >= e ){
    oldx0<-x0
    x0<-oldx0-genD(g,oldx0)$D[1]/genD(g,oldx0)$D[2]
    
    err<-abs(oldx0-x0)
    niter<-niter+1
  }
  return(x0)
}

##Secant function

mySecant <- function(g,a,b,e=10^(-8)){
  maxiter<-1000
  err<-1
  niter<-0
  x0<-a
  x1<-b
  
  while ( niter<=maxiter && err >= e ){
    x2<-x0-genD(g,x0)$D[1]*(x0-x1)/(genD(g,x0)$D[1]-genD(g,x1)$D[1])
    x0<-x1
    x1<-x2
    
    err<-abs(x0-x1)
    niter<-niter+1
  }
  return(x1)
}
