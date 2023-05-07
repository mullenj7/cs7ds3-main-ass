##imports
library(ggplot2)
library(MCMCpack)


##code
data_unfiltered <- read.csv("eight_players.csv")
data_name_unfiltered <- data_unfiltered[c('STOCK','Name','Column1')]


##data<-subset(data_name_unfiltered, Name == 'Bam Adebayo') ## mean(Home > Away) - [1] 0.502
##data<-subset(data_name_unfiltered, Name == 'Jarrett Allen') ## mean(Home > Away) - [1] 0.4922
##data<-subset(data_name_unfiltered, Name == 'Jaylen Brown') ## mean(Home > Away) - [1] 0.4172
##data<-subset(data_name_unfiltered, Name == 'Robin Lopez') ## mean(Home > Away) - [1] 0.528
##data<-subset(data_name_unfiltered, Name == 'Evan Mobley') ## mean(Home > Away) - [1] 0.35
##data<-subset(data_name_unfiltered, Name == 'Jayson Tatum') ## mean(Home > Away) - [1] 0.558
##data<-subset(data_name_unfiltered, Name == 'Nikola Vucevic') ## mean(Home > Away) - [1] 0.61




data$Home.Away <- factor(data$Column1)
summary(data)


##  ----------------------- Gibbs Sampler

compare_2_gibbs <- function(y, ind, mu0 = 2.5, tau0 = 1/4, del0 = 0, gamma0 = 1/4, a0 = 50, b0 = 1, maxiter = 5000)
{
y1 <- y[ind == 'H'] ## is home
y2 <- y[ind == '@'] ## is away

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


fit <- compare_2_gibbs(data$STOCK, data$Home.Away)

summary(fit)
##plot(as.mcmc(fit))
raftery.diag(as.mcmc(fit))
Home <- rnorm(5000, fit[, 1] + fit[, 2], sd = 1/sqrt(fit[, 3]))
Away <- rnorm(5000, fit[, 1] - fit[, 2], sd = 1/sqrt(fit[, 3]))

##ggplot(data.frame(y_sim_diff = Home - Away)) + stat_bin(aes(y_sim_diff))

## probability of home stocks mean being greater than away
mean(Home > Away) 

##ggplot(data.frame(Home, Away)) + geom_point(aes(Home, Away), alpha = 0.3) + geom_abline(slope = 1, intercept = 0)


