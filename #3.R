library(dplyr)

#Kmeans 함수 구현
MyKmeans <- function(x,k,iter.max=50,nstart=1){
  for(j in 1:nstart){
    n <- nrow(x)
    intial_mean <- x[sample(1:n,k),]
    membership <- c()
    sws <- c()
    error <- matrix(rep(1, k*ncol(x)), nrow=k)
    
    niter <- 1
    while(niter <= iter.max & sum(error)>0){

    
      for(i in 1:n){
        membership[i] <- which.min(as.matrix(dist(rbind(x[i,],intial_mean)))[1,-1])
      }
      old_mean <- new_mean
      group_x <- data.frame(x,membership)
      new_mean <- as.matrix(group_x %>% group_by(membership) %>% summarise_all(mean))[,-1]
    
      #if(old_mean[order(old_mean$V1),] == new_mean[order(new_mean$V1),]){break()}
    
      niter <- niter + 1
    }
    sws <- 0.5*sum(sum(new_mean))
    result_n <- as.matrix(group_x %>% group_by(membership) %>% count(ncol(group_x)))[,-2]
  }
  return(list(result_n,new_mean,membership,sws))
}

#sample x
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
colnames(x) <- c("x1", "x2")

#Using MyKmeans
c1 <- MyKmeans(x,2)
plot(x, col = c1[[3]])
points(c1[[2]], col = 1:2, pch = 8, cex = 2)

#Using kmeans in R
c2 <- kmeans(x,2)
plot(x, col=c2$cluster)
points(c2$centers, col = 1:2, pch = 8, cex = 2)

#Using GenSA in model selection
library(GenSA)
baseball <- read.table("C:/대학원/2018-2/1. 전공/통계계산특론1/baseball.txt",header=T)
baseball$freeagent <- factor(baseball$freeagent)
baseball$arbitration <- factor(baseball$arbitration)

baseball.sub <- baseball[,-1]
salary <- baseball$salary
log_salary <- log(baseball$salary)

#salary model
aic_gensa<-function(theta){
  data <- cbind(baseball.sub[,theta>0.5], salary)
  model <- lm(salary ~ .,data)
  extractAIC(model)[2]
}
baseball_gensa <- GenSA(fn=aic_gensa, lower=rep(0,27), upper=rep(1,27), control = list(maxit=100) )

#log_salary model
log_aic_gensa<-function(theta){
  data <- cbind(baseball.sub[,theta>0.5], log_salary)
  model <- lm(log_salary ~ .,data)
  extractAIC(model)[2]
}
log_baseball_gensa <- GenSA(fn=log_aic_gensa, lower=rep(0,27), upper=rep(1,27), control = list(maxit=100) )
plot(log_baseball_gensa$trace.mat[,4])

#Using GA in model selection
library(GA)
aic_ga<-function(theta){
  x<-c(FALSE,theta==1)
  -extractAIC(lm(baseball$salary~., data = baseball[,x]))[2]
}

#salary model
baseball_ga <- ga(type = "binary",fitness=aic_ga, nBits = p, popSize = 100)
#plot(baseball_ga)

#log_salary model
log_aic_ga<-function(theta){
  x<-c(FALSE,theta==1)
  -extractAIC(lm(log(baseball$salary)~., data = baseball[,x]))[2]
}

log_baseball_ga <- ga(type = "binary",fitness=log_aic_ga, nBits = p, popSize = 100)
summary(log_baseball_ga)$solution

#Graph for log_salary best aic model
plot(log_baseball_ga)
