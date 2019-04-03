data <- c(379,299,222,145,109,95,73,59,45,30,24,12,4,2,0,1,1)
N <- sum(data)

#(b)
MyEM <- function(alpha, beta, mu, lambda, maxiter=1000){
  iter <- 1; j <- c(0:16)
  err1 <- 1; err2 <- 1; err3 <- 1; err4 <- 1
  Pi <- c(); t <- c(); p <- c()
  while(err1 > 10^(-10) &err2 > 10^(-10)& err3 > 10^(-10) & err4 > 10^(-10)& iter < maxiter){
    
    #E-step
    for(i in 1:length(data)){
      Pi[i] <- alpha*ifelse(j[i]==0,1,0)+beta*mu^j[i]*exp(-mu)+(1-alpha-beta)*lambda^j[i]*exp(-lambda)
      t[i] <- beta*mu^j[i]*exp(-mu)/Pi[i]
      p[i] <- (1-alpha-beta)*lambda^j[i]*exp(-lambda)/Pi[i]
    }
    z0 <- alpha/Pi[1]
    
    
    #M-step
    new_alpha <- data[1]*z0/N
    new_beta <- sum(data*t)/N
    new_mu <- sum(j*data*t)/sum(data*t)
    new_lambda <- sum(j*data*p)/sum(data*p)
    
    rbind(new_alpha,new_beta,new_mu,new_lambda,iter)
    #Iterate until convergence
    err1 <- abs(alpha - new_alpha)
    err2 <- abs(beta - new_beta)
    err3 <- abs(mu - new_mu)
    err4 <- abs(lambda - new_lambda)
    alpha <- new_alpha; beta <- new_beta; mu <- new_mu; lambda <- new_lambda
    iter <- iter + 1
  }
  rbind(alpha, beta, mu, lambda, iter)
}
alpha <- sum(data[1:3])/N; beta <- sum(data[4:12])/N; mu <- 1; lambda <- 8
MyEM(alpha,beta,mu,lambda)
alpha <- 0.4; beta <- 0.1; mu <- 2; lambda <- 6
MyEM(alpha,beta,mu,lambda)

system.time(MyEM(alpha,beta,mu,lambda))



#------------------------------------------------------------------------------------------------------------#
#(c)
prob <- data/N           

set.seed(2)
result_bs <- c()

for(i in 1:1000) {        
  freq <- rmultinom(1, size=N, prob=prob);
  #Initial Values
  alpha <- 0.6; beta <- 0.38; mu <- 1; lambda <- 8
  j <- c(0:16)
  data_bs <- freq
  
  iter <- 1; maxiter <- 1000; j <- c(0:16)
  err1 <- 1; err2 <- 1; err3 <- 1; err4 <- 1
  Pi <- c(); t <- c(); p <- c()
  
  while(err1 > 10^(-10) &err2 > 10^(-10)& err3 > 10^(-10) & err4 > 10^(-10)& iter < maxiter){
    
    #E-step
    for(i in 1:length(data_bs)){
      Pi[i] <- alpha*ifelse(j[i]==0,1,0)+beta*mu^j[i]*exp(-mu)+(1-alpha-beta)*lambda^j[i]*exp(-lambda)
      t[i] <- beta*mu^j[i]*exp(-mu)/Pi[i]
      p[i] <- (1-alpha-beta)*lambda^j[i]*exp(-lambda)/Pi[i]
    }
    z0 <- alpha/Pi[1]
    
    
    #M-step
    new_alpha <- data_bs[1]*z0/N
    new_beta <- sum(data_bs*t)/N
    new_mu <- sum(j*data_bs*t)/sum(data_bs*t)
    new_lambda <- sum(j*data_bs*p)/sum(data_bs*p)
    
    
    #Iterate until convergence
    err1 <- abs(alpha - new_alpha)
    err2 <- abs(beta - new_beta)
    err3 <- abs(mu - new_mu)
    err4 <- abs(lambda - new_lambda)
    alpha <- new_alpha; beta <- new_beta; mu <- new_mu; lambda <- new_lambda
    iter <- iter + 1
  }
  result_bs <- rbind(result_bs, c(alpha,beta,mu,lambda))
}
cov(result_bs)
sqrt(cov(result_bs))
cor(result_bs)

#------------------------------------------------------------------------------------------------------------#
x <- rep(j, data)
plot(table(x), ylab="Freq",main="Freq of respondents reporting numbers")

lam1 <- mean(x)
dist <- dpois(0:16, lambda = lam1)
dist <- dist*sum(data)
dist <- as.data.frame(dist)
lines(dist, lwd=2)

library(vcd)
df <- data.frame(0:16,data)
colnames(df) <- c("num","freq")
gf <- goodfit(df[, 2:1],"poisson")
plot(gf, type = "standing", scale = "raw")
#------------------------------------------------------------------------------------------------------------#


log_likelihood <- function(param, data){
  p <- c(); result <- c()
  facto <- factorial(0:(length(data)-1))
  
  for(i in 1:length(data)){
    p[i] <- param[1]*ifelse(i==1,1,0)+param[2]*param[3]^(i-1)*exp(-param[3])+(1-param[1]-param[2])*param[4]^(i-1)*exp(-param[4])
  }
  return(-sum(data*log(p/facto)))
}
optim(par=c(0.6,0.38,1,8), log_likelihood, data=data)
optim(par=c(0.4,0.1,2,6), log_likelihood, data=data)


system.time(optim(par=c(0.6,0.38,1,8), log_likelihood, data=data))
system.time(optim(par=c(0.4,0.1,2,6), log_likelihood, data=data))

# optim(par=c(0.13,0.60,1.5,6.12), log_likelihood, data=data)


#------------------------------------------------------------------------------------------------------------#
## GA
library(GA)
GA1 <- ga(type="real-valued", fitness = log_likelihood,min=c(0.1,0.2,0.5,3),max=c(0.5,0.7,2,8),data=data,popSize = 50,maxiter = 20000)
summary(GA1)
plot(GA1)

## GenSA
library(GenSA)
GenSA1 <- GenSA(fn=log_likelihood, lower=c(0.1,0.2,0.5,3), upper=c(0.5,0.7,2,8), data=data,control = list(maxit=1000) )
GenSA1$par
plot(GenSA1$trace.mat[,4],xlab="iter",ylab="Fitness function value")


system.time(ga(type="real-valued", fitness = log_likelihood,min=c(0.1,0.2,0.5,3),max=c(0.5,0.7,2,8),data=data,popSize = 50,maxiter = 20000))
system.time(GenSA(fn=log_likelihood, lower=c(0.1,0.2,0.5,3), upper=c(0.5,0.7,2,8), data=data,control = list(maxit=1000) ))

#------------------------------------------------------------------------------------------------------------#
##Results

j <- c(0:16)

alpha_hat <- 0.1221; beta_hat <- 0.5625; mu_hat <- 1.4675; lambda_hat <- 5.9389
pi_hat <- c(); t_hat <- c(); p_hat <- c()
for(i in 1:length(j)){
  pi_hat[i] <- alpha_hat*ifelse(j[i]==0,1,0)+beta_hat*mu_hat^j[i]*exp(-mu_hat)+(1-alpha_hat-beta_hat)*lambda_hat^j[i]*exp(-lambda_hat)
  
  t_hat[i] <- beta_hat*mu_hat^j[i]*exp(-mu_hat)/pi_hat[i]
  p_hat[i] <- (1-alpha_hat-beta_hat)*lambda_hat^j[i]*exp(-lambda_hat)/pi_hat[i]
}
z0_hat <- alpha_hat/pi_hat[1]

##Group1
z0_hat*data[1]*alpha_hat

##Group2
df1 <- data.frame(0:16,t_hat*data*beta_hat)
colnames(df1) <- c("num","freq")
gf1 <- goodfit(df1[, 2:1],"poisson")
plot(gf1, type = "standing", scale = "raw", main="Typical Group")

##Group3
df2 <- data.frame(0:16,p_hat*data*(1-alpha_hat-beta_hat))
colnames(df2) <- c("num","freq")
gf2 <- goodfit(df2[, 2:1],"poisson")
plot(gf2, type = "standing", scale = "raw", main="Promiscuous Group")

##Total
# z0_hat*data[1]*alpha_hat
data[1]-df1$freq[1]-df2$freq[2]
sum(df1$freq)
sum(df2$freq)
