##imports
library(ggplot2)
library(MCMCpack)


##code
data_unfiltered <- read.csv("Jaren_Jackson_jr_data.csv")
data <- data_unfiltered[c('STOCK','Home.Away')]

##data$Home.Away <- factor(data$Home.Away)
summary(data)

## STOCK           Home.Away

## Min.   :0.000   Away:17  
## 1st Qu.:2.000   Home:16  
## Median :4.000            
## Mean   :4.152            
## 3rd Qu.:6.000            
## Max.   :9.000     

##ggplot(data) + geom_boxplot(aes(Home.Away, STOCK, fill = Home.Away)) + geom_jitter(aes(Home.Away, STOCK, shape = data$Home.Away))

tapply(data$STOCK, data$Home.Away, mean)
## Away     Home 
## 2.882353 5.500000 

tapply(data$STOCK, data$Home.Away, median)
## Away Home 
## 3.0  5.5 

tapply(data$STOCK, data$Home.Away, sd)
## Away    Home 
## 1.99632 2.42212 

t.test(STOCK ~ Home.Away, data=data, var.equal = TRUE)
## Two Sample t-test

## data:  STOCK by Home.Away
## t = -3.3965, df = 31, p-value = 0.001889
## alternative hypothesis: true difference in means between group Away and group Home is not equal to 0
## 95 percent confidence interval:
## -4.189470 -1.045824
## sample estimates:
## mean in group Away      mean in group Home 
## 2.882353                5.500000 




##  ----------------------- Gibbs Sampler

compare_2_gibbs <- function(y, ind, mu0 = 4.5, tau0 = 1/4, del0 = 0, gamma0 = 1/4, a0 = 50, b0 = 1, maxiter = 5000)
{
y1 <- y[ind == 'Home']
y2 <- y[ind == 'Away']

n1 <- length(y1) 
n2 <- length(y2)

##### starting values
mu <- (mean(y1) + mean(y2)) / 2
del <- (mean(y1) - mean(y2)) / 2

mat_store <- matrix(0, nrow = maxiter, ncol = 3)
#####

##### Gibbs sampler
an <- a0 + (n1 + n2)/2

for(s in 1 : maxiter) 
{
  
  ##update tau
  bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
  tau <- rgamma(1, an, bn)
  ##
  
  ##update mu
  taun <-  tau0 + tau * (n1 + n2)
  mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
  mu <- rnorm(1, mun, sqrt(1/taun))
  ##
  
  ##update del
  gamman <-  gamma0 + tau*(n1 + n2)
  deln <- ( del0 * gamma0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
  del<-rnorm(1, deln, sqrt(1/gamman))
  ##
  
  ## store parameter values
  mat_store[s, ] <- c(mu, del, tau)
}
colnames(mat_store) <- c("mu", "del", "tau")
return(mat_store)
}


fit <- compare_2_gibbs(data$STOCK, as.factor(data$Home.Away))

summary(fit)
##plot(as.mcmc(fit))
raftery.diag(as.mcmc(fit))
Home <- rnorm(5000, fit[, 1] + fit[, 2], sd = 1/sqrt(fit[, 3]))
Away <- rnorm(5000, fit[, 1] - fit[, 2], sd = 1/sqrt(fit[, 3]))

##ggplot(data.frame(y_sim_diff = Home - Away)) + stat_bin(aes(y_sim_diff))

## probability of home stocks mean being greater than away
mean(Home > Away) 
##[1] 0.9502
ggplot(data.frame(Home, Away)) + geom_point(aes(Home, Away), alpha = 0.3) + geom_abline(slope = 1, intercept = 0)


