library(survival) 
library(survminer) 
library(dplyr)
library(ggplot2) 


current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
getwd()




#### 1. read data and basic preprocessing ####

data = read.csv('wilms.csv')
colnames(data)
table(data$study)

## 1-1) take only the 3rd study group
data = data[data$study == 3,]
data = data %>% select(-study) # after taking only the 3rd study group, variable "study" is not required.


## 1-2) basic EDA for each variable
str(data)
summary(data$trel)
summary(data$tsur)
table(data$relaps)
table(data$dead)
table(data$stage)
table(data$histol)
summary(data$age)
summary(data$tumdiam)


## 1-3) Define the observed time $X_i$ and the failure indicator $\delta_i$

#### censoring indicator
data$delta = 0
data[data$relaps != 0 | data$dead != 0, 'delta'] = 1

table(data$delta) ; table(data$delta) / nrow(data) # table of censored vs observed

#### observed time
data$eventtime = 0 # define the container
for(i in 1:nrow(data)){
  data$eventtime[i] = min(data$trel[i], data$tsur[i])
}


## 1-4) NA count for each variable
apply(data, 2, function(x) sum(is.na(x))) # no NA variable


## 1-5) learn about the duplicates
duplicated = table(data$eventtime)[table(data$eventtime) > 1]
data %>% filter(data$eventtime %in% as.numeric(names(duplicated))) %>% nrow # number of duplication
length(duplicated) # out of total number of duplication obtained above, 235 unique observed times
table(duplicated) # 202 : 2 duplicates, 27 : 3 duplicates, 4 : 4 duplicates, 2 : 5 duplicates


## may or may not scale the numeric variables
# data$age = scale(data$age)[,1]
# data$tumdiam = scale(data$tumdiam)[,1]







#### 2. model fitting ####

## 2-0) Before model fitting, drop useless variables for analysis

colnames(data)
data = data %>% select(-trel, -tsur, -relaps, -dead)
colnames(data)



## 2-1) As a preliminary, for numeric predictors, check functional forms of covariates (Therneau et al., 1990)

### the following commented codes do not properly work (do not know the reason why)
# par(mfrow = c(1,2))
# 
# fit0 = coxph(Surv(eventtime, delta) ~ 1, ties = "efron", data = data)
# resid0 = resid(fit0)
# 
# plot(data$age, resid0, xlab = "age", ylab = "Residual")
# lines(lowess(data$age, resid0), col = "red")
# 
# plot(data$tumdiam, resid0, xlab = "tumdiam", ylab = "Residual")
# lines(lowess(data$tumdiam, resid0), col = "red")


## about age
ggcoxfunctional(coxph(Surv(eventtime, delta) ~ age , data=data), ylim = c(-1,1),data=data)
ggcoxfunctional(coxph(Surv(eventtime, delta) ~ I(age^2) , data=data), ylim = c(-1,1),data=data)

## about tumdiam
ggcoxfunctional(coxph(Surv(eventtime, delta) ~ tumdiam ,data=data), ylim = c(-1,1),data=data)
ggcoxfunctional(coxph(Surv(eventtime, delta) ~ I(tumdiam^2) ,data=data), ylim = c(-1,1),data=data)


## opted to consider age and tumdiam^2.







## 2-2) Fit full model and see the summary

fit1 = coxph(Surv(eventtime, delta) ~ stage + histol + age + I(tumdiam^2),
             ties = "efron", data = data)

summary(fit1)


## 2-3) check the reduced model without the age variable and see the summary

fit2 = coxph(Surv(eventtime, delta) ~ stage + histol + I(tumdiam^2),
             ties = "efron", data = data)
summary(fit2)


lrt_age = summary(fit1)$logtest - summary(fit2)$logtest
pchisq(lrt_age[1], lrt_age[2], lower.tail = F) # the null model is better




## 2-4) variable to use are chosen. Then, check the proportionality hazard assumption

par(mfrow = c(1,3))
test_PH_fit2 = cox.zph(fit2)
test_PH_fit2
plot(test_PH_fit2)

ggcoxdiagnostics(fit2, type="dfbeta", linear.predictions = F, ggtheme = theme_bw()) # no influential points, but serious proportionality hazard assumption violation



## 2-5) piecewise regression coefficients by dividing time.

splitdata = survSplit(Surv(eventtime, delta) ~. , data = data, cut = c(0.6, 1.3), episode = "tgroup", id = "Patient")
head(splitdata)
nrow(data) ; nrow(splitdata)

fit3 = coxph(Surv(tstart, eventtime, delta) ~ stage:strata(tgroup) + histol:strata(tgroup) + I(tumdiam^2), ties = "efron", data = splitdata)
summary(fit3)


lrt_timedivision = summary(fit3)$logtest - summary(fit2)$logtest
pchisq(lrt_timedivision[1], lrt_timedivision[2], lower.tail = F) # time dividing is better



test_PH_fit3 = cox.zph(fit3, terms = F)
test_PH_fit3


ggcoxdiagnostics(fit3, type="dfbeta", linear.predictions = F, ggtheme = theme_bw()) # no influential points, but serious proportionality hazard assumption violation


