# HW1
#initial values
nC <- 85; nI <- 196; nT <- 341; N <- nC+nI+nT
p <- 1/3; q <- 1/3; r <- 1/3

iter <- 1; maxiter <- 1000; err <- 1
while(err > 10^(-10) & iter < maxiter){
  #E-step
  En_CC <- nC*p^2/(p^2+2*p*q+2*p*r)
  En_CI <- nC*2*p*q/(p^2+2*p*q+2*p*r)
  En_CT <- nC*2*p*r/(p^2+2*p*q+2*p*r)
  En_II <- nI*q^2/(q^2+2*q*r)
  En_IT <- nI*2*q*r/(q^2+2*q*r)
  
  
  #M-step
  new_p <- (2*En_CC+En_CI+En_CT)/(2*N)
  new_q <- (2*En_II+En_IT+En_CI)/(2*N)
  new_r <- (2*nT+En_CT+En_IT)/(2*N)
  
  
  #Iterate until convergence
  err <- (abs(p - new_p) + abs(q - new_q) + abs(r - new_r))/sqrt(p^2+q^2+r^2)
  p <- new_p; q <- new_q; r <- new_r
  iter <- iter + 1
}
cbind(p, q, r, err, iter)
cbind(p, q, r, iter, En_CC, En_CI, En_CT, En_II, En_IT)

MyRiemann <- function(f, a, b, e=10^(-5)){
  k <- 5; iter <- 1; err <- 1; old_R <- (b-a)*f(a)/32
  while(err > e & iter < 23){
    h <- (b-a)/2^k; i <- c(0:2^k-1)
    R <- h*sum(f(a+i*h))
    err <- abs(R-old_R)
    old_R <- R
    iter <- iter+1; k <- k+1
  }
  return(list(R=R,iter=iter, err=err))
}

MyTrap <- function(f, a, b, e=10^(-5)){
  k <- 5; iter <- 1; err <- 1; old_Tr <- (b-a)*f(a)/32
  while(err > e & iter < 23){
    h <- (b-a)/2^k; i <- c(1:2^k-1)
    Tr <- h*f(a)/2+h*sum(f(a+i*h))+h*f(b)/2
    err <- abs(Tr-old_Tr)
    old_Tr <- Tr
    iter <- iter+1; k <- k+1
  }
  return(list(Tr=Tr,iter=iter, err=err))
}

MySimpson <- function(f, a, b, e=10^(-5)){
  k <- 5; iter <- 1; err <- 1; old_S <- (b-a)*f(a)/32
  while(err > e & iter < 25){
    h <- (b-a)/2^k; i <- c(1:2^(k-1))
    S <- h*(f(a)-f(b)+4*sum(f(a+(2*i-1)*h))+2*sum(f(a+(2*i)*h)))/3
    err <- abs(S-old_S)
    old_S <- S
    iter <- iter+1; k <- k+1
  }
  return(list(S=S,iter=iter, err=err))
}

##eg1
f1 <- function(x){x^2}
plot(seq(-10,10,0.01),f1(seq(-10,10,0.01)),type="l",main="Graph of f1(x)",xlab="x",ylab="f1")

MyRiemann(f1, 0,1); system.time(MyRiemann(f1, 0,1))
MyTrap(f1,0,1); system.time(MyTrap(f1, 0,1))
MySimpson(f1,0,1); system.time(MySimpson(f1,0,1))
integrate(f1, lower = 0, upper = 1); system.time(integrate(f1, lower = 0, upper = 1))

##eg2
f2 <- function(x){(x^2+5*x+4*x^3-1)/sqrt(abs(x)+4*x^5)}

plot(seq(-10,10,0.01),f2(seq(-10,10,0.01)),type="l",main="Graph of f2(x)",xlab="x",ylab="f2")
plot(seq(-1,10,0.0001),f2(seq(-1,10,0.0001)),type="l",main="Graph of f2(x)",xlab="x",ylab="f2")
f2(0)

MyRiemann(f2,0.001,10)
system.time(MyRiemann(f2,0.001,10))
MyTrap(f2,0.001,10)
system.time(MyTrap(f2,0.001,10))
MySimpson(f2,0.001,10)
system.time(MySimpson(f2,0.001,10))
integrate(f2, lower = 0, upper = 10)$subdivisions
system.time(integrate(f2, lower = 0, upper = 10))

