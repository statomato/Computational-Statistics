# HW 1
S0 <- 100; K <- 102; sigma <- 0.3; N <- 50; r <- 0.05
n <- 100000; m <- 100
mu_mc <- c(); theta_mc <- c(); A <- c(); theta <- c()
lambda <- -1


for(j in 1:m){
  for(i in 1:n){
    temp <- exp((r-sigma^2/2)/365+sigma*rnorm(N-1)/sqrt(365))
    S <- S0*cumprod(temp)
    A[i] <- exp(-r*N/365)*max(0,mean(S)-K)
    theta[i] <- exp(-r*N/365)*max(c(0,exp(mean(log(S))) - K))
  }
  mu_mc[j] <- mean(A)
  theta_mc[j] <- mean(theta)
}

c3 <- 1+1/N
c2 <- sigma*(c3*N*(1+1/(2*N))/1095)^(1/2)
c1 <- (log(S0/K)+(c3*N/730*(r-sigma^2/2)+c3*sigma^2*N/1095*(1+1/(2*N))))/c2

theta0 <- S0*pnorm(c1)*exp(-N*(r+c3*sigma^2/6)*(1-1/N)/730)-K*pnorm(c1-c2)*exp(-r*N/365)

mu_cv <- mu_mc + lambda*(theta_mc-theta0)
cbind(mean(mu_mc),mean(mu_cv))
cbind(sd(mu_mc),sd(mu_cv))
# cor(mu_cv,theta_mc)

#HW 2

## a
q <- function(x){exp(-abs(x)^3/3)}
g <- function(x){dnorm(x)}
h <- function(x){x^2}

plot(seq(-4,4,0.0001),q(seq(-4,4,0.0001)),type="l",main="plot of q(x), g(x)",xlab="x",ylab="q(x)")
lines(seq(-4,4,0.0001),g(seq(-4,4,0.0001)),type="l",lwd=2,lty=2)
legend(2.5,1,c("q(x)","g(x)"),lty=1:2,lwd=2)

plot(seq(-10,10,0.0001),q(seq(-10,10,0.0001))/g(seq(-10,10,0.0001)),type="l",main="plot of q(x)/g(x)",xlab="x",ylab="q(x)/g(x)")

n <- 10^5
X <- rnorm(n)
w_star <- q(X)/g(X)
w <- w_star/sum(w_star)
mu_IS <- sum(h(X)*w)
mu_IS

## b
n <- 10^5
X <- rnorm(n); U <-runif(n)
xx <- q(X)/(sqrt(2*pi)*g(X))
Y <-X[U<=xx]
length(Y)/n; 1-length(Y)/n
mu_b <- mean(h(Y))
mu_b

plot(seq(-4,4,0.0001),q(seq(-4,4,0.0001)),type="l",main="Rejection Sampling",xlab="x",ylab="y")
lines(seq(-4,4,0.0001),g(seq(-4,4,0.0001))*sqrt(2*pi),type="l",lwd=2,lty=2)
legend(1.8,1,c("target","env"),lty=1:2,lwd=2)

## c
X_sort <- sort(X)
mu_c <- sum((X_sort[2:n]-X_sort[1:n-1])*h(X_sort[1:n-1])*q(X_sort[1:n-1]))/sum((X_sort[2:n]-X_sort[1:n-1])*q(X_sort[1:n-1]))
mu_c

## d
n <- 10^5
result <- matrix(nrow=100, ncol=2)
for(j in 1:100){
  X <- rnorm(n); U <-runif(n)
  xx <- q(X)/(sqrt(2*pi)*g(X))
  Y <-X[U<=xx]
  length(Y)/n; 1-length(Y)/n
  mu_b <- mean(h(Y))
  result[j,1] <- mu_b
  
  X_sort <- sort(X)
  mu_c <- sum((X_sort[2:n]-X_sort[1:n-1])*h(X_sort[1:n-1])*q(X_sort[1:n-1]))/sum((X_sort[2:n]-X_sort[1:n-1])*q(X_sort[1:n-1]))
  result[j,2] <- mu_c
}
mu_result <- apply(result,2,mean)
var_result <- apply(result,2,var)

mu_result; var_result

#HW 3

## a

### Naive MC
ld <- 2; n <- 1000
Z <- c(); mu_MC <- c()
# X <- matrix(rpois(m*n, ld), ncol = m)
for(j in 1:1000){
  for(i in 1:1000){
    Z[i] <- (mean(rpois(25,ld))-2)/sqrt(2/25)
  }
  hZ <- ifelse(Z>1.645,1,0)
  mu_MC[j] <- sum(hZ)/n
}
# mu_MC
mean(mu_MC); sd(mu_MC)
quantile(mu_MC,0.05); quantile(mu_MC,0.95)

### antithetic MC
mu_AS <- c()
zz <- function(x){ifelse((mean(x)-2)/sqrt(2/25)>1.645,1,0)}
for(i in 1:1000){
  u1 <- matrix(runif(25*(n/2)), ncol=25) ; u2 <- 1 - u1
  X1 <- qpois(u1,ld) ; X2 <- qpois(u2,ld)
  hZ1 <- apply(X1,1,zz) ; hZ2 <- apply(X2,1,zz)
  mu_AS[i] <- mean(hZ1 + hZ2)/2
}
mean(mu_AS); sd(mu_AS)
quantile(mu_AS,0.05); quantile(mu_AS,0.95)

### IS unstandardized & standardized
p <- function(x){dpois(x,ld*25)}
g <- function(x){dpois(x,2.4653*25)}
n <- 1000
mu_unIS <- c(); mu_IS <- c(); mu_ISCV <- c()

for(j in 1:1000){
  X <- matrix(rpois(25*n,2.4653),ncol=25)
  Y <- apply(X,1,sum)
  w_star <- p(Y)/g(Y)
  w <- w_star/sum(w_star)
  h <- ifelse(Y>=25*(1.645*sqrt(2/25)+2),1,0)
  mu_unIS[j] <- mean(h*w_star)
  mu_IS[j] <- sum(h*w)
  
  lam <- -lm(h*w_star~w_star)$coef[2]
  mu_ISCV[j] <- mu_IS[j] + lam*(mean(w_star)-1)
}
mean(mu_unIS); sd(mu_unIS)
mean(mu_IS); sd(mu_IS)
mean(mu_ISCV); sd(mu_ISCV)

quantile(mu_unIS,0.05); quantile(mu_unIS,0.95)
quantile(mu_IS,0.05); quantile(mu_IS,0.95)
quantile(mu_ISCV,0.05); quantile(mu_ISCV,0.95)


## b
lam1<-seq(2.2,4,length=5)
power<-function(lam){  
  m<-matrix(nrow=length(lam),ncol=5)  ; colnames(m)<-c("standard","AS","IS_un","IS_s","IS_c")  
  rownames(m)<-lam ;s<-m  
  for(j in 1:length(lam)){    
    n<-5000    
    res<-matrix(ncol=100,nrow=5)    
    result<-matrix(ncol=2,nrow=5)    
    for(i in 1:100){      
      ## Naive MC      
      z<-function(x){ifelse((mean(x)-2)/sqrt(2/25)>1.645,1,0)}      
      x<-matrix(rpois(25*n,lam[j]),ncol=25)      
      h<-apply(x,1,z)      
      res[1,i]<-mean(h)            
      
      ## AS      
      u1<-matrix(runif(25*(n/2),0,1),ncol=25) ; u2<-1-u1      
      x1<-qpois(u1,lam[j]) ; x2<-qpois(u2,lam[j])      
      h1<-apply(x1,1,z) ; h2<-apply(x2,1,z)      
      res[2,i]<-mean(h1+h2)/2            
      
      ##IS      
      f.target<-function(x){dpois(x,lambda=lam[j]*25)}      
      g.env<-function(x){dpois(x,lambda=(lam[j]+0.4653)*25)}            
      
      x.s<-matrix(rpois(25*n,lam[j]+0.4653),ncol=25)      
      y<-apply(x.s,1,sum)      
      w.s<-f.target(y)/g.env(y)      
      w<-w.s/sum(w.s)      
      h3<-ifelse(y>=25*(1.645*sqrt(2/25)+2),1,0)            
      
      #unstd'd IW      
      res[3,i]<-mean(h3*w.s)      
      
      #std'd IW      
      res[4,i]<-sum(h3*w)            
      
      #control      
      lambda<--lm(h3*w.s~w.s)$coef[2]      
      res[5,i]<-mean(h3*w.s)+lambda*(mean(w.s)-1)    
      
    }    
    result[,1:2]<-c(apply(res,1,mean),apply(res,1,sd))   
    m[j,]<-result[,1] ; s[j,]<-result[,2]  
    }  
  result2<-list(mean=m,sd=s) 
  return(result2) 
  } 
result_6.6_b<-power(lam1) 


library(gplots) 
plotCI(lam1,result_6.6_b$mean[,1],uiw=1.96*result_6.6_b$sd[,1],liw=1.96*result_6.6_b$sd[,1],pch=20,       
       gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of Naive MC") 
plotCI(lam1,result_6.6_b$mean[,2],uiw=1.96*result_6.6_b$sd[,2],liw=1.96*result_6.6_b$sd[,2],pch=20,       
       gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of Antithetic Sampling") 
plotCI(lam1,result_6.6_b$mean[,3],uiw=1.96*result_6.6_b$sd[,3],liw=1.96*result_6.6_b$sd[,3],pch=20,       
       gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of IS with unstd'd") 
plotCI(lam1,result_6.6_b$mean[,4],uiw=1.96*result_6.6_b$sd[,4],liw=1.96*result_6.6_b$sd[,4],pch=20,       
       gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of IS with std'd") 
plotCI(lam1,result_6.6_b$mean[,5],uiw=1.96*result_6.6_b$sd[,5],liw=1.96*result_6.6_b$sd[,5],pch=20,       
       gap=0,type="o",cex=0.95,xlab=expression(lambda),ylab="power",main="Power Curve of IS with control variate")

      

#HW 4

S0 <- 50; K <- 52; sigma <- 0.5; N <- 30; r <- 0.05

## a
European_call <- function(S0, K, sigma, N, r){
  Z <- rnorm(10^5)
  St <- S0*exp((r-sigma^2/2)*N/365 + sigma*Z*sqrt(N/365))
  C <- exp(-r*N/365)*sapply(St,function(x){max(0,x-K)})
  
  mu_mc <- mean(C)
  mu_sd <- sd(C)/sqrt(10^5)
  return(list(mu_mc=mu_mc, mu_sd=mu_sd))
}
European_call(S0,K,sigma,N,r)

## b
Asian_call <- function(S0, K, sigma, N, r){
  mu_mc <- c(); mu_as <- c(); theta <- c(); theta_mc <- c(); A <- c(); A_as1 <- c(); A_as2 <- c()
  # lambda <- -1
  for(j in 1:100){
    for(i in 1:1000){
      ## MC
      temp1 <- exp((r-sigma^2/2)/365+sigma*rnorm(N-1)/sqrt(365))
      S <- S0*cumprod(temp1)
      A[i] <- exp(-r*N/365)*max(0,mean(S)-K)
      
      ## AS
      temp2 <- exp((r-sigma^2/2)/365+sigma*(-rnorm(N-1))/sqrt(365))
      S_as1 <- S0*cumprod(temp1); S_as2 <- S0*cumprod(temp2)
      A_as1[i] <- exp(-r*N/365)*max(0,mean(S_as1)-K);  A_as2[i] <- exp(-r*N/365)*max(0,mean(S_as2)-K)
      
      ## CV
      theta[i] <- exp(-r*N/365)*max(c(0,exp(mean(log(S))) - K))
    }
    mu_mc[j] <- mean(A)
    theta_mc[j] <- mean(theta)
    mu_as[j] <- (mean(A_as1)+mean(A_as2))/2
  }
  
  c3 <- 1+1/N
  c2 <- sigma*(c3*N*(1+1/(2*N))/1095)^(1/2)
  c1 <- (log(S0/K)+(c3*N/730*(r-sigma^2/2)+c3*sigma^2*N/1095*(1+1/(2*N))))/c2
  
  theta0 <- S0*pnorm(c1)*exp(-N*(r+c3*sigma^2/6)*(1-1/N)/730)-K*pnorm(c1-c2)*exp(-r*N/365)
  lambda <- -cov(mu_mc,theta_mc)/var(theta_mc)
  
  mu_cv <- mu_mc + lambda*(theta_mc-theta0)
  sd_cv <- sd(mu_cv); sd_mc <- sd(mu_mc); sd_as <- sd(mu_as)
  return(list(mu_mc=mean(mu_mc),mu_cv=mean(mu_cv),mu_as=mean(mu_as),sd_mc=sd_mc,sd_cv=sd_cv,sd_as=sd_as,lambda=lambda))
}
Asian_call(S0,K,sigma,N,r)


