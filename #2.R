library(leaps); library(ggplot2); library(corrplot)
baseball <- read.table("C:/대학원/2018-2/1. 전공/통계계산특론1/baseball.txt",header=T)

## 변수 탐색
str(baseball)
qplot(baseball$salary, geom="histogram") 
ggplot(data=baseball, aes(x=1, y=salary))+geom_boxplot()
baseball$freeagent <- factor(baseball$freeagent)
baseball$arbitration <- factor(baseball$arbitration)


## salary log 변환 후 탐색
baseball$ln_salary <- log(baseball$salary)
qplot(baseball$ln_salary, geom="histogram") 
ggplot(data=baseball, aes(x=1, y=ln_salary))+geom_boxplot()

M<-cor(baseball[,-c(14,15,29)])
corrplot(M, method="color")

##stepwise model
step_model <- step(lm(salary~average+obp+runs+hits+doubles+triples+homeruns+rbis+walks+sos+sbs+errors+factor(freeagent)+factor(arbitration)+runsperso+hitsperso+hrsperso+rbisperso+walksperso+obppererror+runspererror+hitspererror+hrspererror+soserrors+sbsobp+sbsruns+sbshits,data=baseball),direction = "both")
summary(lm(salary ~ runs + hits + rbis + sos + sbs + freeagent + arbitration + 
     runsperso + hitsperso + hrsperso + rbisperso + walksperso + 
     soserrors + sbsobp, data=baseball))

##stepwise ln_model
step_lnmodel <-step(lm(ln_salary~average+obp+runs+hits+doubles+triples+homeruns+rbis+walks+sos+sbs+errors+freeagent+arbitration+runsperso+hitsperso+hrsperso+rbisperso+walksperso+obppererror+runspererror+hitspererror+hrspererror+soserrors+sbsobp+sbsruns+sbshits,data=baseball),direction = "both")


##all possible regressions model 
p <- ncol(baseball)-2
regsubsets_model <- regsubsets(salary~average+obp+runs+hits+doubles+triples+homeruns+rbis+walks+sos+sbs+errors+freeagent+arbitration+runsperso+hitsperso+hrsperso+rbisperso+walksperso+obppererror+runspererror+hitspererror+hrspererror+soserrors+sbsobp+sbsruns+sbshits,data=baseball,method=c("exhaustive", "backward", "forward", "seqrep"),nvmax = p)
AIC(regsubsets(salary~average+obp+runs+hits+doubles+triples+homeruns+rbis+walks+sos+sbs+errors+freeagent+arbitration+runsperso+hitsperso+hrsperso+rbisperso+walksperso+obppererror+runspererror+hitspererror+hrspererror+soserrors+sbsobp+sbsruns+sbshits,data=baseball,method=c("exhaustive", "backward", "forward", "seqrep"),nbest=5))


summary(regsubsets_model)$cp
summary(regsubsets_model)$which

AIC(lm(salary~homeruns+rbis+walks+sos+freeagent+arbitration+walksperso+hrspererror+soserrors+sbsobp,data=baseball))

AIC(lm(salary~homeruns+rbis+walks+sos+freeagent+arbitration+walksperso+sbsobp,data=baseball))

##all possible regressions ln_model 
regsubsets_lnmodel <- regsubsets(ln_salary~average+obp+runs+hits+doubles+triples+homeruns+rbis+walks+sos+sbs+errors+freeagent+arbitration+runsperso+hitsperso+hrsperso+rbisperso+walksperso+obppererror+runspererror+hitspererror+hrspererror+soserrors+sbsobp+sbsruns+sbshits,data=baseball,method=c("exhaustive", "backward", "forward", "seqrep"),nvmax = p)

summary(regsubsets_lnmodel)$cp
summary(regsubsets_lnmodel)$which

AIC(lm(ln_salary~obp + runs + triples + rbis + sos + freeagent + arbitration + runsperso + hitsperso + soserrors + sbsobp + sbsruns,data=baseball))


##result
result <- lm(ln_salary ~ average + runs + triples + rbis + sos + freeagent + 
     arbitration + runsperso + hitsperso + soserrors + sbsobp + 
     sbsruns, data=baseball)
anova(result)
summary(result)



####extract AIC
extractAIC(lm(ln_salary ~ obp + runs + triples + rbis + sos + freeagent + arbitration + runsperso + hitsperso + soserrors + sbsobp + sbsruns,data=baseball))
