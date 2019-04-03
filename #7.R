#HW1
library(extraDistr)
m <- 100000; n<- 5000
y_norm <- rnorm(m); y_slash <- rslash(m) 

w1 <- function(x){
  r1 <- dslash(x)/dnorm(x)
  r2 <- sum(r1)
  return(r1/r2)
}
w2 <- function(x){
  r1 <- dnorm(x)/dslash(x)
  r2 <- sum(r1)
  return(r1/r2)
}

## weights & resample
weights = w1(y_norm)
SIR1 <- sample(y_norm, n, replace=TRUE, prob=weights)
weights2 = w2(y_slash)
SIR2 <- sample(y_slash, n, replace=TRUE, prob=weights2)

## graph
par(mfrow=c(1,2))
hist(SIR2,freq=FALSE,breaks=seq(-7,7,by=.25),main="Histogram of draws",
     ylab="Normal density")
points(seq(-10,10,by=.01),dnorm(seq(-10,10,by=.01)),type="l")

hist(SIR1,freq=FALSE,breaks=seq(-7,7,by=.25),main="Histogram of draws",
     ylab="Slash density")
points(seq(-10,10,by=.01),dslash(seq(-10,10,by=.01)),type="l")


#HW2
##Rejection Sampling of Poisson
MyRS_poi <- function(obs){
  n <- 10^5
  z <- rlnorm(n,log(4),0.5); u<-runif(length(z),0,1)
  e <- prod(dpois(obs,mean(obs)))
  q <- sapply(z,function(x){prod(dpois(obs,x))})
  z<-z[u<q/e]
  xx<-seq(0,20,by=0.01)
  target_y<-dlnorm(xx,log(4),0.5)*sapply(xx,function(x){prod(dpois(obs,x))})*10^11
  env_y<-dlnorm(xx,log(4),0.5)*prod(dpois(obs,mean(obs)))*10^11
  par(mfrow=c(1,1))
  plot(xx,env_y,type="l",lwd=2,xlab=expression(lambda),ylab="Unnormalized Density x 10^11")
  lines(xx,target_y,type="l",lwd=2,lty=2)
  legend(12,3.1,c("Envelope","Target"),lty=1:2,lwd=2)
  reject_rate <- length(z)/n; mean <- mean(z); var <- var(z)
  return(list(reject_rate=reject_rate, mean=mean, var=var, z=z))
}

obs1 <- c(8,3,4,3,1,7,2,6,2,7)
RS1 <- MyRS_poi(obs1)
RS1$reject_rate; RS1$mean; RS1$var
RS1$z

n <- 10000; m <- 100000
lambda <- rlnorm(m, log(4), 0.5)

## weight & resample
w3 <- sapply(lambda, function(x){prod(dpois(obs1,x))})
w3 <- w3/sum(w3)
SIR3 <- sample(lambda, n , replace=T, prob=w3)
hist(SIR3, prob=T,breaks=50,yaxs="i",xlab=expression(lambda),main="Histogram of SIR")
mean(SIR3); var(SIR3)

## qqplot
par(mfrow=c(1,2))
qqnorm(RS1$z, main = "Reject Sampling Q-Q Plot"); qqline(RS1$z)
qqnorm(SIR3, main = "SIR Q-Q Plot"); qqline(SIR3)

#HW3
myroute<-function(x){
  m<-matrix(c(	1,x[1],x[2],x[3],x[4],0,0,0,0,0,
               x[1],1,x[5],0,0,x[8],x[9],0,0,0,
               x[2],x[5],1,x[6],0,0,x[10],0,0,0,
               x[3],0,x[6],1,x[7],0,0,x[11],0,0,
               x[4],0,0,x[7],1,0,0,x[12],x[13],0,
               0,x[8],0,0,0,1,x[14],0,0,x[17],
               0,x[9],x[10],0,0,x[14],1,x[15],0,x[18],
               0,0,0,x[11],x[12],0,x[15],1,x[16],x[19],
               0,0,0,0,x[13],0,0,x[16],1,x[20],
               0,0,0,0,0,x[17],x[18],x[19],x[20],1),ncol=10,byrow=T)
  res<-ifelse((m%*%m%*%m%*%m%*%m%*%m%*%m%*%m%*%m)[1,10]>0,0,1)
  return(res)
}
# myplot<-function(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20){
#   plot(c(0,1,1,1,1,2,2,2,2,3),c(0,1.5,0.5,-0.5,-1.5,1.5,0.5,-0.5,-1.5,0),axes=FALSE,xlab="",ylab="",pch=16)
#   text(0,-0.1,"a")
#   text(3,-0.1,"b")
#   if(x1==1){lines(c(0,1),c(0,1.5))}
#   if(x2==1){lines(c(0,1),c(0,0.5))}
#   if(x3==1){lines(c(0,1),c(0,-0.5))}
#   if(x4==1){lines(c(0,1),c(0,-1.5))}
#   if(x5==1){lines(c(1,1),c(1.5,0.5))}
#   if(x6==1){lines(c(1,1),c(0.5,-0.5))}
#   if(x7==1){lines(c(1,1),c(-0.5,-1.5))}
#   if(x8==1){lines(c(1,2),c(1.5,1.5))}
#   if(x9==1){lines(c(1,2),c(1.5,0.5))}
#   if(x10==1){lines(c(1,2),c(0.5,0.5))}
#   if(x11==1){lines(c(1,2),c(-0.5,-0.5))}
#   if(x12==1){lines(c(1,2),c(-1.5,-0.5))}
#   if(x13==1){lines(c(1,2),c(-1.5,-1.5))}
#   if(x14==1){lines(c(2,2),c(1.5,0.5))}
#   if(x15==1){lines(c(2,2),c(0.5,-0.5))}
#   if(x16==1){lines(c(2,2),c(-0.5,-1.5))}
#   if(x17==1){lines(c(2,3),c(1.5,0))}
#   if(x18==1){lines(c(2,3),c(0.5,0))}
#   if(x19==1){lines(c(2,3),c(-0.5,0))}
#   if(x20==1){lines(c(2,3),c(-1.5,0))}
# }

# myroute(c(1,1,0,0,0,0,1,0,1,0,1,1,1,0,1,0,1,0,0,1))
# myplot(1,1,0,0,0,0,1,0,1,0,1,1,1,0,1,0,1,0,0,1)
# myroute(c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0))
# myplot(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0)
# myroute(c(1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,0,0,0,1))
# myplot(1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,0,0,0,1)



## MC estimating of mu
p <- 0.05; p_star <- 0.2; n<- 10^6
x_MC <- matrix(rbinom(20*n,1,1-p),ncol=20)
hx_MC <- apply(x_MC,1,myroute)
mu_MC <- mean(hx_MC); var_MC <- mu_MC*(1-mu_MC)/n


## IS_star estimating of mu
x_star <- matrix(rbinom(20*n,1,1-p_star),ncol=20)
hx_star <- apply(x_star,1,myroute)
bx_star <- 20-apply(x_star,1,sum)  
w_star<-((1-p)/(1-p_star))^20 * (p*(1-p_star)/(p_star*(1-p)))^bx_star
mu_star <- mean(hx_star*w_star) 
var_star <- (sum(w_star[hx_star==1]*p^bx_star[hx_star==1]*(1-p)^(20-bx_star[hx_star==1]))-mu_star^2)/n

w<-w_star/sum(w_star)
mu_ <- sum(hx_star*w)  
var_<-(var(hx_star*w_star)+mu_star^2*var(w_star)-2*mu_star*cov(hx_star*w_star,w_star))/n

list(mu_MC=mu_MC,var_MC=var_MC,mu_star=mu_star,var_star=var_star,mu_=mu_,var_=var_)

######################################################
## MC estimating of mu
p <- 0.05; p_star <- 0.4; n<- 10^6

## IS_star estimating of mu
x_star <- matrix(rbinom(20*n,1,1-p_star),ncol=20)
hx_star <- apply(x_star,1,myroute)
bx_star <- 20-apply(x_star,1,sum)  
w_star<-((1-p)/(1-p_star))^20 * (p*(1-p_star)/(p_star*(1-p)))^bx_star
mu_star <- mean(hx_star*w_star) 
var_star <- (sum(w_star[hx_star==1]*p^bx_star[hx_star==1]*(1-p)^(20-bx_star[hx_star==1]))-mu_star^2)/n

w<-w_star/sum(w_star)
mu_ <- sum(hx_star*w)  
var_<-(var(hx_star*w_star)+mu_star^2*var(w_star)-2*mu_star*cov(hx_star*w_star,w_star))/n

list(mu_MC=mu_MC,var_MC=var_MC,mu_star=mu_star,var_star=var_star,mu_=mu_,var_=var_)

######################################################
## MC estimating of mu
p <- 0.05; p_star <- 0.5; n<- 10^6

## IS_star estimating of mu
x_star <- matrix(rbinom(20*n,1,1-p_star),ncol=20)
hx_star <- apply(x_star,1,myroute)
bx_star <- 20-apply(x_star,1,sum)  
w_star<-((1-p)/(1-p_star))^20 * (p*(1-p_star)/(p_star*(1-p)))^bx_star
mu_star <- mean(hx_star*w_star) 
var_star <- (sum(w_star[hx_star==1]*p^bx_star[hx_star==1]*(1-p)^(20-bx_star[hx_star==1]))-mu_star^2)/n

w<-w_star/sum(w_star)
mu_ <- sum(hx_star*w)  
var_<-(var(hx_star*w_star)+mu_star^2*var(w_star)-2*mu_star*cov(hx_star*w_star,w_star))/n

list(mu_MC=mu_MC,var_MC=var_MC,mu_star=mu_star,var_star=var_star,mu_=mu_,var_=var_)
