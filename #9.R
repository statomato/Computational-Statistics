library(ggplot2)
#HW 1.1

# data <- read.csv("C:/대학원/2018-2/1. 전공/통계계산특론1/coal.csv")
# y <- data$disasters 

n <- 100; d<-0.7
y <- c(rnorm(d*n, 7, 0.5), rnorm((1-d)*n, 10, 0.5))

##(a) Beta(1,1)
f <- function(x){prod(x*dnorm(y,7,0.5)+(1-x)*dnorm(y,10,0.5))}
g <- function(x){dbeta(x,1,1)}
R <- function(xt,x_star){f(x_star)*g(xt)/(f(xt)*g(x_star))}


N <- 10000


d1 <- c(); d1[1] <- runif(1)
for(i in 1:N){
  xt <- d1[i]
  x <- rbeta(1,1,1)
  p <- min(R(xt,x),1)
  d <- rbinom(1,1,p)
  d1[i+1] <- x*d + xt*(1-d)
}
plot(d1,ylim=c(0,1),type="l",ylab=expression(delta^(t)),xlab="t",xlim=c(0,10000),main="Sample paths with Beta(1,1)")
ggplot() + geom_histogram(aes(x=d1)) + labs(title="Histogram of Independence chains with Beta(1,1)",x=expression(delta),y="Frequency")

summary(d1[501:(n+1)])





##(b) Beta(2,10)
f <- function(x){prod(x*dnorm(y,7,0.5)+(1-x)*dnorm(y,10,0.5))}
g <- function(x){dbeta(x,2,10)}
R <- function(xt,x_star){f(x_star)*g(xt)/(f(xt)*g(x_star))}


d2 <- c(); d2[1] <- rbeta(1,2,10)
for(i in 1:N){
  xt <- d2[i]
  x <- rbeta(1,2,10)
  p <- min(R(xt,x),1)
  d <- rbinom(1,1,p)
  d2[i+1] <- x*d + xt*(1-d)
}
plot(d2,ylim=c(0,1),type="l",ylab=expression(delta^(t)),xlab="t",xlim=c(0,10000),main="Sample paths with Beta(2,10)")
ggplot() + geom_histogram(aes(x=d2)) + labs(title="Histogram of Independence chains with Beta(2,10)",x=expression(delta),y="Frequency")

summary(d2[501:(n+1)])
 



#HW 1.2
likelihood <- function(x){prod(x*dnorm(y,7,0.5)+(1-x)*dnorm(y,10,0.5))}
J <- function(u){(1+exp(u))^2/exp(u)}

n <- 100000
inv_logit <- function(u){exp(u)/(1+exp(u))}



##(a) b=1
u1 <- c(); d1 <- c(); b <- 1
u1[1] <- runif(1); d1[1] <- inv_logit(u1[1])
for(i in 1:n){
  u_star <- u1[i] + runif(1,-b,b)
  d_star <- inv_logit(u_star)
  MH_R <- likelihood(d_star)*J(u_star)/(likelihood(d1[i])*J(u1[i]))
  p <- min(MH_R,1)
  d <- rbinom(1,1,p)
  u1[i+1] <- u_star*d + u1[i]*(1-d)
  d1[i+1] <- exp(u1[i+1])/(1+exp(u1[i+1]))
}
plot(d1,ylim=c(0,1),type="l",ylab=expression(delta^(t)),xlab="t",xlim=c(0,10000),main="Sample paths with Unif(-1,1)")
ggplot() + geom_histogram(aes(x=d1)) + labs(title="Histogram of Random walks chains with Unif(-1,1)",x=expression(delta),y="Frequency")

summary(d1[501:(n+1)])

##(b) b=0.01
u2 <- c(); d2 <- c(); b <- 0.01
u2[1] <- runif(1); d2[1] <- inv_logit(u1[1])
for(i in 1:n){
  u_star <- u2[i] + runif(1,-b,b)
  d_star <- inv_logit(u_star)
  MH_R <- likelihood(d_star)*J(u_star)/(likelihood(d2[i])*J(u2[i]))
  p <- min(MH_R,1)
  d <- rbinom(1,1,p)
  u2[i+1] <- u_star*d + u2[i]*(1-d)
  d2[i+1] <- exp(u2[i+1])/(1+exp(u2[i+1]))
}
plot(d2,ylim=c(0,1),type="l",ylab=expression(delta^(t)),xlab="t",xlim=c(0,10000),main="Sample paths with Unif(-0.01,0.01)")
ggplot() + geom_histogram(aes(x=d2)) + labs(title="Histogram of Random walks chains with Unif(-0.01,0.01)",x=expression(delta),y="Frequency")

summary(d2[501:(n+1)])


#HW 2
data2 <- read.table("C:/대학원/2018-2/1. 전공/통계계산특론1/Datasets/Datasets/breastcancer.txt",head=T)

#a
summary(data2)
summary(data2$recurtime[data2$treatment==1])
summary(data2$recurtime[data2$treatment==0])
ggplot(data2,aes(x=factor(treatment),y=recurtime))+geom_boxplot()+stat_summary(fun.y="mean",geom="point",shape=22,size=3,fill="blue")
ggplot(data2,aes(x=factor(treatment),y=recurtime,fill=factor(censored)))+geom_boxplot()+stat_summary(fun.y="mean",geom="point",shape=22,size=3,fill="blue")

#c
a <- 3; b <- 1; c <- 60; d <- 120
n <- 100000
time_c <- data2$recurtime[data2$treatment==0]
time_h <- data2$recurtime[data2$treatment==1]

data2$delta <- abs(data2$censored-1)

## theta & tau
theta <- c(); theta[1] <- 0
tau <- c(); tau[1] <- 1

for(i in 1:n){
  theta[i+1] <- rgamma(1,sum(data2$censored)+a+1,sum(time_c)+c+tau[i]*(sum(time_h)+d))
  tau[i+1] <- rgamma(1,sum(data2$censored[data2$treatment==1])+b+1,theta[i+1]*(sum(time_h)+d))
}
result1 <- list(theta=theta,tau=tau)
theta1.1 <- result1$theta[-(1:500)]; tau1.1 <- result1$tau[-(1:500)]

plot(result1$theta,main="sample path of theta",type="l",xlab="t",ylab=expression(theta^(t)))
acf(theta1.1,main="ACF of theta")

plot(result1$tau,main="sample path of tau",type="l",xlab="t",ylab=expression(tau^(t)))
acf(tau1.1,main="ACF of tau")


#d
post <- theta1.1^(sum(delta)+a)*tau1.1^(sum(delta_h)+b)*exp(-theta1.1*(sum(x_c)+c)-tau1.1*theta1.1*(sum(x_h)+d))
summary(post)
summary(theta1.1)
summary(tau1.1)

mean(theta1.1); sd(theta1.1)
quantile(theta1.1,0.05); quantile(theta1.1,0.95)

mean(tau1.1); sd(tau1.1)
quantile(tau1.1,0.05); quantile(tau1.1, 0.95)

#e
tau_prior <- function(tau,a,b,c,d){
  tau^b*gamma(a+1)/(c+d*tau)^(a+1)
}

x<-seq(0,6,by=0.01) 
plot(x,tau_prior(x,a,b,c,d)/integrate(tau_prior,0,Inf,a=a,b=b,c=c,d=d)$value,type="l",xlim=c(0,6), ylim=c(0,3),ylab="dens ity",xlab=expression(tau),main="density of tau") 
lines(density(result1$tau),lwd=2,lty=2,col=2) 
legend(4,3,c("Prior","Posterior"),lwd=1:2,lty=1:2,col=1:2)

#f
quantile(tau1.1,0.05); quantile(tau1.1,0.95)

#g
##(1)
a <- 3/2; b <- 1/2; c <- 60/2; d <- 120/2

theta <- c(); theta[1] <- 0
tau <- c(); tau[1] <- 1

for(i in 1:n){
  theta[i+1] <- rgamma(1,sum(data2$censored)+a+1,sum(time_c)+c+tau[i]*(sum(time_h)+d))
  tau[i+1] <- rgamma(1,sum(data2$censored[data2$treatment==1])+b+1,theta[i+1]*(sum(time_h)+d))
}
result1 <- list(theta=theta,tau=tau)
theta1.1 <- result1$theta[-(1:500)]; tau1.1 <- result1$tau[-(1:500)]
quantile(tau1.1,0.05); quantile(tau1.1,0.95)

##(2)
a <- 3*2; b <- 1*2; c <- 60*2; d <- 120*2

theta <- c(); theta[1] <- 0
tau <- c(); tau[1] <- 1

for(i in 1:n){
  theta[i+1] <- rgamma(1,sum(data2$censored)+a+1,sum(time_c)+c+tau[i]*(sum(time_h)+d))
  tau[i+1] <- rgamma(1,sum(data2$censored[data2$treatment==1])+b+1,theta[i+1]*(sum(time_h)+d))
}
result1 <- list(theta=theta,tau=tau)
theta1.1 <- result1$theta[-(1:500)]; tau1.1 <- result1$tau[-(1:500)]
quantile(tau1.1,0.05); quantile(tau1.1,0.95)
