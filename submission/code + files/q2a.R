##imports
library(ggplot2)
library(MCMCpack)


##code
data_unfiltered <- read.csv("eight_players.csv")
data <- data_unfiltered[c('STOCK','Column1', 'Name')]

data$Home.Away <- factor(data$Column1)
data$Name <- factor(data$Name)
summary(data)

##      STOCK          Column1                      Name     Home.Away 
## Min.   : 0.000   Length:552         Evan Mobley   : 75   @   :273  
## 1st Qu.: 1.000   Class :character   Nikola Vucevic: 75   NA's:279  
## Median : 2.000   Mode  :character   Robin Lopez   : 73             
## Mean   : 2.257                      Bam Adebayo   : 71             
## 3rd Qu.: 3.000                      Jayson Tatum  : 70             
## Max.   :10.000                      Jarrett Allen : 65             
                                     (Other)       :123                         
##ggplot(data.frame(size = tapply(data$STOCK, data$Name, length), 
##                  mean_score = tapply(data$STOCK, data$Name, mean)), 
##       aes(size, mean_score)) + geom_point()


compare_m_gibbs <- function(y, ind, maxiter = 5000)
{
  
### weakly informative priors
a0 <- 1/2 ; b0 <- 50 ## tau_w hyperparameters
eta0 <- 1/2 ; t0 <- 50 ## tau_b hyperparameters
mu0 <- 2.257 ; gamma0 <- 1/4

###

### starting values
m <- nlevels(ind)
ybar <- theta <- tapply(y, ind, mean)
print(ybar)
tau_w <- mean(1 / tapply(y, ind, var)) ##within group precision
mu <- mean(theta)
tau_b <-var(theta) ##between group precision
n_m <- tapply(y, ind, length)
an <- a0 + sum(n_m)/2
###

### setup MCMC
theta_mat <- matrix(0, nrow=maxiter, ncol=m)
mat_store <- matrix(0, nrow=maxiter, ncol=3)
###

### MCMC algorithm
for(s in 1:maxiter) 
{
  
  # sample new values of the thetas
  for(j in 1:m) 
  {
    taun <- n_m[j] * tau_w + tau_b
    thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
    theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
  }
  
  #sample new value of tau_w
  ss <- 0
  for(j in 1:m){
    ss <- ss + sum((y[ind == j] - theta[j])^2)
  }
  bn <- b0 + ss/2
  tau_w <- rgamma(1, an, bn)
  
  #sample a new value of mu
  gammam <- m * tau_b + gamma0
  mum <- (mean(theta) * m * tau_b + mu0 * gamma0) / gammam
  mu <- rnorm(1, mum, 1/ sqrt(gammam)) 
  
  # sample a new value of tau_b
  etam <- eta0 + m/2
  tm <- t0 + sum((theta - mu)^2) / 2
  tau_b <- rgamma(1, etam, tm)
  
  #store results
  theta_mat[s,] <- theta
  mat_store[s, ] <- c(mu, tau_w, tau_b)
}
colnames(mat_store) <- c("mu", "tau_w", "tau_b")
return(list(params = mat_store, theta = theta_mat))
}

fit2 <- compare_m_gibbs(data$STOCK, data$Name)
apply(fit2$params, 2, mean)
##        mu         tau_w      tau_b 
##        2.26675552 5.52262170 0.07920673 

## reformat samples for ggplot
theta_df <- data.frame(STOCKS = as.numeric(fit2$theta), 
Name = rep(1:ncol(fit2$theta), each = nrow(fit2$theta))) 
theta_med <- apply(theta_df, 2, mean) ## get basic posterior summary
print(theta_med)
##    Name     STOCKS
##    4.500000 2.287582 


theta_hat <- apply(fit2$theta, 2, mean)

result=data.frame(size = tapply(data$STOCK, data$Name, length),
                  theta_hat = theta_hat)
print(result)

ggplot(theta_df) + geom_boxplot(aes(x = reorder(Name, STOCKS, median), STOCKS, 
                               fill = reorder(Name, STOCKS, median)), show.legend=FALSE)






