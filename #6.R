#HW2
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
  while(err > e & iter < 24){
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

#a

f1 <- function(mu){
  data1 <- c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
  xbar <- mean(data1)
  dnorm(xbar, mu, 3/sqrt(7))*dcauchy(mu,5,2)
}
plot(seq(-10,10,0.01),f1(seq(-10,10,0.01)), type="l", xlab="mu", ylab="prior * likelihood")
(result <- 1/7.84654)
MyRiemann(f1, 0, 10); system.time(MyRiemann(f1, 0, 10))
MyTrap(f1, 0, 10); system.time(MyTrap(f1, 0, 10))
MySimpson(f1, 0, 10); system.time(MySimpson(f1, 0, 10))
integrate(f1, lower = 0, upper = 10); system.time(integrate(f1, lower = 0, upper = 10))

#b
f2 <- function(mu){
  data1 <- c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
  xbar <- mean(data1)
  7.84654*dnorm(xbar, mu, 3/sqrt(7))*dcauchy(mu,5,2)
}
MyRiemann(f2, 2, 8); system.time(MyRiemann(f2, 2, 8))
MyTrap(f2, 2, 8); system.time(MyTrap(f2, 2, 8))
MySimpson(f2, 2, 8); system.time(MySimpson(f2, 2, 8))
integrate(f2, lower = 2, upper = 8); system.time(integrate(f2, lower = 2, upper = 8))

#c
f3 <- function(u){
  data1 <- c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
  xbar <- mean(data1)
  mu <- log(u/(1-u))
  7.84654*dnorm(xbar, mu, 3/sqrt(7))*dcauchy(mu,5,2)/(u*(1-u))
}

integrate(f2, -Inf, Inf)
integrate(f3, exp(3)/(1+exp(3)), 1)
MyRiemann(f3, exp(3)/(1+exp(3)), 1-10^(-7)); system.time(MyRiemann(f3, exp(3)/(1+exp(3)), 1-10^(-7)))
MyTrap(f3, exp(3)/(1+exp(3)), 1-10^(-7)); system.time(MyTrap(f3, exp(3)/(1+exp(3)), 1-10^(-7)))
MySimpson(f3, exp(3)/(1+exp(3)), 1-10^(-7)); system.time(MySimpson(f3, exp(3)/(1+exp(3)), 1-10^(-7)))

#d
f4 <- function(u){
  data1 <- c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
  xbar <- mean(data1)
  mu <- 1/u
  -7.84654*dnorm(xbar, mu, 3/sqrt(7))*dcauchy(mu,5,2)/(u^2)
}

integrate(f4, 1/3, 0)
MyRiemann(f4, 1/3, 0); system.time(MyRiemann(f4, 1/3, 0))
MyTrap(f4, 1/3, 0+10^(-10)); system.time(MyTrap(f4, 1/3,0+10^(-10)))
MySimpson(f4, 1/3, 0+10^(-10)); system.time(MySimpson(f4, 1/3, 0+10^(-10)))

#HW3
Romberg <- function(f, a, b, n){
  T0 <- (f(a)+f(b))/2; h <- (b-a)/n
  TT <- matrix(nrow=(n+1),ncol=(n+1))
  TT[1,1] <- T0
  for(i in 1:n){
    TT[i+1,1] <- (TT[i,1] + h*sum(f(a+(1:i-0.5)*h)))/2
    for(j in 1:n){
     TT[i+1,j+1] <- (4^(j+1)*TT[i+1,j]-TT[i,j])/(4^(j+1)-1) 
    }
  }
  return(TT)
}
f5 <- function(x){1/x}
a <- 2:10
result <- list()
for(i in 1:length(a)){
  print(Romberg(f5, 1, a[i], 6))
  print(cat(a[i],":",log(a[i]),'\n'))
}

#HW4
##Rejection Sampling of Gamma deviates
MyRS <- function(r){
  a <- r/3; b<- 1/sqrt(9*a)
  z <- rnorm(20000); u <- runif(20000)
  t <- function(y){a*(1+b*y)^3}
  ge <- exp(z^2/2+a*(log(t(z))-log(a))-t(z)+a)
  result <- ifelse(u<ge, 1, 0)
  sum(result, na.rm=T)/length(u)
}

MyRS(1)
MyRS(4)
MyRS(10)
MyRS(20)

qqplot()












