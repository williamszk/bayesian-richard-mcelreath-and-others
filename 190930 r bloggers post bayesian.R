
#https://www.r-bloggers.com/bayesian-models-in-r-2/
#https://www.r-bloggers.com/bayesian-models-in-r-2/
#https://www.r-bloggers.com/bayesian-models-in-r-2/



library(rethinking)
library(tidyverse)
library(magrittr)
library(readxl)

# ================================================================================

rangeP <- seq(0, 1, length.out = 100)
plot(rangeP, dbinom(x = 8, prob = rangeP, size = 10), type = "l", xlab = "P(Black)", ylab = "Density")
lines(rangeP, dnorm(x = rangeP, mean = .5, sd = .1) / 15, col = "red")

lik <- dbinom(x = 8, prob = rangeP, size = 10)
prior <- dnorm(x = rangeP, mean = .5, sd = .1)
post <- lik * prior
lines(rangeP, post, col = "green")

postst <- post / sum(post)
lines(rangeP, postst, col = "blue")
legend("topleft", legend = c("Lik", "Prior", "Unstd Post", "Post"),
       text.col = 1:4, bty = "n")

# simulation =================================================================================
#DGP
set.seed(1)
sample <- rnorm(100,mean=5,sd=2^2)
#grid exploration
rmu = seq(0,10,length.out=200)
rsig = seq(1,3,length.out=200)
pnorm(0,1)
dnorm(0,mean=0,sd=1)
parmat = matrix(nrow=200,ncol=200)
rownames(parmat) <- substr(as.character(rmu),start = 1,stop = 4) 
colnames(parmat) <- substr(as.character(rsig),start = 1,stop = 4)  
mu=0
sig=1
for (ii in 1:200) {
  mu <- rmu[ii]
  for (jj in 1:200) {
    sig <- rsig[jj]
    lik <- sum(log(dnorm(sample,mean=mu,sd=sig^.5)))
    parmat[ii,jj] <- lik
  }
}

library(plot3D)
persp3D(z = parmat, theta = 120)

#product of likelihood and priors
prior_mu <- dnorm(rmu,0,5)
plot(rmu,prior_mu, type='l')
prior_sig <- dexp(rsig,1)
plot(rsig,prior_sig, type='l')

#an example of kronecker product 
# simple scalar multiplication
M <- matrix(1:6, ncol = 2) 
kronecker(4, M)
# Block diagonal matrix:
kronecker(diag(1, 3), M)
# ask for dimnames
fred <- matrix(1:12, 3, 4, dimnames = list(LETTERS[1:3], LETTERS[4:7]))
bill <- c("happy" = 100, "sad" = 1000)
class(bill)
kronecker(fred, bill, make.dimnames = TRUE)
bill <- outer(bill, c("cat" = 3, "dog" = 4))
kronecker(fred, bill, make.dimnames = TRUE)

#multiply likelihood by prior matrix
prior_matrix <- outer(prior_mu, prior_sig,FUN=function(x,y){x+y})
outer(c(1,2,3),c(1,1,1), FUN=function(x,y){x+y})
prior_matrix[1:10,1:10]
dim(prior_matrix)

post <- parmat+prior_matrix
post[1:10,1:10]
persp3D(z = post, theta = 120)

#https://www.r-bloggers.com/bayesian-models-in-r-2/
## from the site
# Define real pars mu and sigma, sample 100x
trueMu <- 5
trueSig <- 2
set.seed(100)
randomSample <- rnorm(100, trueMu, trueSig)
# Grid approximation, mu %in% [0, 10] and sigma %in% [1, 3]
grid <- expand.grid(mu = seq(0, 10, length.out = 200),
                    sigma = seq(1, 3, length.out = 200))
class(grid)
# Compute likelihood
lik <- sapply(1:nrow(grid), function(x){
  sum(dnorm(x = randomSample, mean = grid$mu[x],
            sd = grid$sigma[x], log = T))
})
class(lik)
length(lik)
# Multiply (sum logs) likelihood and priors
prod <- lik + dnorm(grid$mu, mean = 0, sd = 5, log = T) +
  dexp(grid$sigma, 1, log = T)

# Standardize the lik x prior products to sum up to 1, recover unit
prob <- exp(prod - max(prod))

# Sample from posterior dist of mu and sigma, plot
postSample <- sample(1:nrow(grid), size = 1e2, prob = prob)
#?sample()
sample_mu <- grid$mu[postSample]
sample_sig <- grid$sigma[postSample]
plot(sample_mu,sample_sig,
     xlab = "Mu", ylab = "Sigma", pch = 16, col = rgb(0,0,0,.2))
abline(v = trueMu, h = trueSig, col = "red", lty = 2)
median(sample_mu)
median(sample_sig)

#Bayesian models & MCMC
#Let’s get started with R

library(rethinking)
library(tidyverse)
library(magrittr)
library(readxl)

# Download data set from Riehl et al. 2019
cuckoo <- read_xlsx("Riehl and Strong_Social Parasitism Data_2007-2017_DRYAD.xlsx",
                    sheet=2)
class(cuckoo)
str(cuckoo)

#Zero-inflated Poisson regression of fledged egg counts
summary(cuckoo$Eggs_fledged)
hist(cuckoo$Eggs_fledged)
dim(cuckoo)
cuckoo2 <- cuckoo[!is.na(cuckoo$Eggs_fledged),]
dim(cuckoo2)

str(cuckoo2)
cuckoo3 <- cuckoo2 %>% 
  mutate(
    female_id = as.integer(factor(Female_ID_coded)),
    year_id = as.integer(factor(Year)),
    group_id = as.integer(factor(Group_ID_coded)),
    Min_age_Z = scale(Min_age),
    Group_size_Z = scale(Group_size),
    Mean_eggsize_Z = scale(Mean_eggsize))
str(cuckoo3)

ss1 <- sample(x=1:500,size=20)
cuckoo3[ss1,] %>% View

hist(cuckoo3$Min_age,breaks=25)
hist(cuckoo3$Group_size,breaks=25)

#map2stan only accepts data frames or lists
cuckoo4 <- as.data.frame(cuckoo3)

###
#https://www.r-bloggers.com/bayesian-models-in-r-2/
eggsFMod <- map2stan(alist(
  Eggs_fledged ~ dzipois(p, lambda),
  logit(p) <- ap,
  log(lambda) <- a + 
    a_fem[female_id] + a_year[year_id] + a_group[group_id] +
    Parasite*bP + Min_age_Z*bA + Group_size_Z*bGS + 
    Mean_eggsize_Z*bES +
    Parasite*Min_age_Z*bPA,
  Group_size_Z ~ dnorm(0, 3),
  Mean_eggsize_Z ~ dnorm(0, 3),
  a_fem[female_id] ~ dnorm(0, sigma1),
  a_year[year_id] ~ dnorm(0, sigma2),
  a_group[group_id] ~ dnorm(0, sigma3),
  c(sigma1, sigma2, sigma3) ~ dcauchy(0, 1),
  c(ap, a) ~ dnorm(0, 3),
  c(bP, bA, bGS, bES, bPA) ~ dnorm(0, 2)),
  data = cuckoo4,
  iter = 5e3, warmup = 1e3, chains = 4, cores = 2)

# Check posterior dists
precis(eggsFMod, prob = .95) # use depth = 2 for varying intercepts
inter1 <- precis(eggsFMod, prob = .95, depth=2) 
class(inter1)

# Sample posterior ======================
# 
post <- extract.samples(eggsFMod)
post %>% class
post[[1]] %>% dim
names(post)
post$bGS %>% dim
aa <- post$bGS
aa %>% class
hist(aa, breaks = 100)
# PI of P(no clutch at all)
dens(logistic(post$ap), show.HPDI = T, xlab = "ZIP Bernoulli(p)")

# Run simulations w/ averages of all predictors, except parasite 0 / 1
lambdaNoP <- exp(post$a + 0*post$bP + 0*post$bA +
                   0*post$bGS + 0*post$bES + 0*0*post$bPA)
simFledgeNoPar <- rpois(n = length(lambdaNoP), lambda = lambdaNoP)
aa <- post$bP
plot(aa[1:200], type='l')
lambdaNoP %>% class
hist(simFledgeNoPar)

lambdaP <- exp(post$a + 1*post$bP + 0*post$bA +
                 0*post$bGS + 0*post$bES + 1*0*post$bPA)
simFledgePar <- rpois(n = length(lambdaP), lambda = lambdaP)
hist(simFledgePar)
table(simFledgeNoPar)
table(simFledgePar)

# Simulate with varying age
rangeA <- seq(-3, 3, length.out = 100)
# No parasite
predictions <- sapply(rangeA, function(x){
  exp(post$a + 0*post$bP + x*post$bA + 0*post$bGS +
        0*post$bES + 0*x*post$bPA)
})
predictions %>% dim

hdpiPois <- apply(predictions, 2, HPDI, prob = .95)
hdpiPois %>% dim
?HPDI()
df1 <- as.data.frame(t(hdpiPois))
names(df1) <- c('lower','upper')

ggplot(data=df1) +
  geom_line(aes(x=1:100,y=upper))+
  geom_line(aes(x=1:100,y=lower))


meanPois <- colMeans(predictions)
meanPois %>% class

plot(rangeA, meanPois, type = "l", ylim = c(0, 3), yaxp = c(0, 3, 3),
     xlab = "Min Age (standardized)", ylab = expression(lambda))
shade(hdpiPois, rangeA)

# Parasite
predictionsP <- sapply(rangeA, function(x){
  exp(post$a + 1*post$bP + x*post$bA + 0*post$bGS +
        0*post$bES + x*post$bPA)
})
predictionsP %>% class
hdpiPoisP <- apply(predictionsP, 2, HPDI, prob = .95)
meanPoisP <- colMeans(predictionsP)
lines(rangeA, meanPoisP, lty = 2, col = "red")
shade(hdpiPoisP, rangeA, col = rgb(1,0,0,.25))


#Poisson regression of laid egg counts =====================================









































