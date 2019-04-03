#HW1
MyGmm <- function(X,threshold=10^(-10),maxiter=10000){
  #initial values
  old_pi <- 0.5; n <- length(X)
  mu <- X[sample(1:n,2)]; mu[1] <- min(mu); mu[2] <- max(mu)
  sd <- c()
  sd[1] <- sd(X); sd[2] <- sd(X)
  Y <- c(); Ey <- c()
  err <- 1; iter <- 1
  
  while(err > threshold & iter < maxiter){
    #E-step
    Ey <-old_pi*dnorm(X, mu[1], sd[1])/(old_pi*dnorm(X, mu[1], sd[1])+(1-old_pi)*dnorm(X, mu[2], sd[2]))
    Y <- ifelse(Ey>0.5,1,0)
    
    #M-step
    new_pi <- mean(Ey)
    mu[1] <- sum(Ey*X)/sum(Ey); mu[2] <- sum((1-Ey)*X)/sum(1-Ey)
    sd[1] <- sqrt(sum(Ey*(X-mu[1])^2)/sum(Ey)); sd[2] <- sqrt(sum((1-Ey)*(X-mu[2])^2)/sum(1-Ey))
    
    #Iterate until convergence
    err <- abs(old_pi-new_pi)
    old_pi <- new_pi
    iter <- iter+1
  }
  result <- list(Y=Y, new_pi=new_pi, mu=mu, sd=sd, iter=iter)
  result[[1]] <- Y; result[[2]] <- new_pi; result[[3]] <- mu; result[[4]] <- sd; result[[5]] <- iter
  return(result)
}
X <- c(rnorm(10,5,1),rnorm(10,10,2))
MyGmm(X)

#HW2
w1 <- c(8,11,16,18,6,4,20,25,9,13)
w2 <- c(10,14,16,15,20,4,18,22,NA,NA)
data <- data.frame(cbind(w1,w2))
Ew <- c(); Ew2 <- c()

#initial values
mu1 <- mean(w1); mu2 <- mean(w2, na.rm=T)
s11 <- var(w1); s22 <- var(w2, na.rm=T); s12 <- cov(w1, w2, use="complete.obs"); rho <- s12/sqrt(s11*s22)
err <- 1; iter <- 1; maxiter <- 1000; n <- length(w1)

while(err > 10^(-10) & iter < maxiter){
  #E-step
  Ew[1] <- mu2 + s12/s11*(w1[9]-mu1)
  Ew[2] <- mu2 + s12/s11*(w1[10]-mu1)
  Ew2[1] <- (Ew[1])^2 + s22*(1-rho^2)
  Ew2[2] <- (Ew[2])^2 + s22*(1-rho^2)
  
  w2[9] <- Ew[1]; w2[10] <- Ew[2]
  ww <- c(w2[1:8]^2, Ew2)
  
  #M-step
  new_mu2 <- mean(w2)
  new_s22 <- (sum(ww)- (sum(w2)^2)/n)/n
  new_s12 <- (sum(w1*w2) - sum(w1)*sum(w2)/n)/n
  new_rho <- new_s12/sqrt(s11*new_s22)
  
  #Iterate until convergence
  err <- sum(abs(mu2-new_mu2)+abs(s22-new_s22))
  mu2 <- new_mu2; s22 <- new_s22; s12 <- new_s12; rho <- new_rho
  iter <- iter+1
}
cbind(new_mu2, new_s22, new_s12, iter)

#HW3
library(dplyr)
library(randomForest)
data3 <- read.csv("C:/대학원/2018-2/1. 전공/통계계산특론1/Concrete_Data.csv")
colnames(data3) <- c("Cement","Slag","Ash","Water","Superplasticizer","Coarse","Fine","Age","strength")

n <- nrow(data3); q <- seq(0.01, 0.1, 0.01)      # q: missing
iter <- 1

origin_Y <- data3$strength
missing <- list(); for(i in 1:10) missing[[i]] <- sample(n, n*q[i])


result_lm <- data.frame(q=q, error=0, n.iter=0, MAE=0)

for(i in 1:10){
  
  data3$strength <- origin_Y
  data3$strength[missing[[i]]] <- mean(origin_Y[-missing[[i]]])
  err <- 1; iter <- 0
  
  while(err > 10^(-10) & iter < 1000){
    
    lm.step <- step(lm(strength ~ ., data3))
    old_Y <- data3$strength[missing[[i]]]
    data3$strength[missing[[i]]] <- fitted(lm.step)[missing[[i]]]
    
    err <- sum(abs(old_Y-data3$strength[missing[[i]]]))
    iter <- iter + 1
  }
  
  result_lm[i, -1] <- c(err, iter, mean(abs(origin_Y-data3$strength)))
}

result_lm

result_rf <- data.frame(q=q, error=0, n.iter=0, MAE=0)

for(i in 7:10){
  
  data3$strength <- origin_Y
  data3$strength[missing[[i]]] <- mean(origin_Y[-missing[[i]]])
  err <- 1; iter <- 0
  
  while(err > 10^(-10) & iter < 200){
    
    rf <- randomForest(strength ~ ., data3)
    old_Y <- data3$strength[missing[[i]]]
    data3$strength[missing[[i]]] <- predict(rf)[missing[[i]]]
    
    err <- sum(abs(old_Y-data3$strength[missing[[i]]]))
    iter <- iter + 1
  }
  
  result_rf[i, -1] <- c(err, iter, mean(abs(origin_Y-data3$strength)))
}

result_rf


#HW4

#initial values
nO <- 176; nA <- 182; nB <- 60; nAB <- 17; N <- 435
old_p <- nA/N; old_q <- nB/N; old_r <- 1 -(old_p + old_q)

iter <- 1; maxiter <- 1000; err1 <- 1; err2 <- 1
while(err1 > 10^(-10) & err2 > 10^(-10) & iter < maxiter){
  #E-step
  En_AA <- nA*old_p^2/(old_p^2+2*old_p*old_r)
  En_AO <- nA*2*old_p*old_r/(old_p^2+2*old_p*old_r)
  En_BB <- nB*old_q^2/(old_q^2+2*old_q*old_r)
  En_BO <- nB*2*old_q*old_r/(old_q^2+2*old_q*old_r)
  
  new_nA <- En_AA + 0.5*En_AO + 0.5*nAB
  new_nB <- En_BB + 0.5*En_BO + 0.5*nAB
  new_nO <- nO + 0.5*En_AO + 0.5*En_BO
  
  #M-step
  new_p <- new_nA/N; new_q <- new_nB/N
  new_r <- 1 - (new_p + new_q)
  
  #Iterate until convergence
  err1 <- abs(old_p - new_p); err2 <- abs(old_q - new_q)
  old_p <- new_p; old_q <- new_q; old_r <- new_r
  iter <- iter + 1
}
cbind(new_p, new_q, new_r)
