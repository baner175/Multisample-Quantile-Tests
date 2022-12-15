rm(list = ls())
library(ggpubr)

rm(list = ls())

n_samp <- 200
n_rep <- 1.5e4

options('contrasts')

source('assymptotic null functions.R')

set.seed(80531)
#-----------------------------------------------------------------------------

# N(0,1) distribution:

cval <- 30
scl <- 1


samp_0.05 <- IQR_cont_norm(eps = 0.05, loc = c(1,2,3),
                           scl = 1, n_rep = n_rep,
                           n_samp = n_samp, cval = cval)$statistic

samp_0.1 <- IQR_cont_norm(eps = 0.1, loc = c(1,2,3),
                          scl = 1, n_rep = n_rep,
                          n_samp = n_samp, cval = cval)$statistic

samp_0.15 <- IQR_cont_norm(eps = 0.15, loc = c(1,2,3),
                           scl = 1, n_rep = n_rep,
                           n_samp = n_samp, cval = cval)$statistic

par(mfrow = c(2,2))

qqplot(rchisq(n_rep , df = 2), samp_0.05)
qqline(samp_0.05, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)

qqplot(rchisq(n_rep , df = 2), samp_0.1)
qqline(samp_0.1, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)

qqplot(rchisq(n_rep , df = 2), samp_0.15)
qqline(samp_0.15, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)

par(mfrow = c(1,1))

title(main = 'Samples from N(0,1)')

#-----------------------------------------------------------------------------

# t_5 distribution:


cval <- 30
scl <- 1


samp_0.05 <- IQR_cont_t(eps = 0.05, loc = c(1,2,3),
                        deg = 5, n_rep = n_rep,
                        n_samp = n_samp, cval = cval)$statistic

samp_0.1 <- IQR_cont_t(eps = 0.1, loc = c(1,2,3),
                       deg = 5, n_rep = n_rep,
                       n_samp = n_samp, cval = cval)$statistic

samp_0.15 <- IQR_cont_t(eps = 0.15, loc = c(1,2,3),
                        deg = 5, n_rep = n_rep,
                        n_samp = n_samp, cval = cval)$statistic

par(mfrow = c(2,2))

qqplot(rchisq(n_rep , df = 2), samp_0.05)
qqline(samp_0.05, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)

qqplot(rchisq(n_rep , df = 2), samp_0.1)
qqline(samp_0.1, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)

qqplot(rchisq(n_rep , df = 2), samp_0.15)
qqline(samp_0.15, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)

par(mfrow = c(1,1))

title(main = bquote(t[5]))

#-----------------------------------------------------------------------------

# Cauchy(0,1) distribution:

cval <- 30
scl <- 1


samp_0.05 <- IQR_cont_cauchy(eps = 0.05, loc = c(1,2,3),
                             scl = 1, n_rep = n_rep,
                             n_samp = n_samp, cval = cval)$statistic
samp_0.1 <- IQR_cont_cauchy(eps = 0.1, loc = c(1,2,3),
                            scl = 1, n_rep = n_rep,
                            n_samp = n_samp, cval = cval)$statistic

samp_0.15 <- IQR_cont_cauchy(eps = 0.15, loc = c(1,2,3),
                             scl = 1, n_rep = n_rep,
                             n_samp = n_samp, cval = cval)$statistic
par(mfrow = c(2,2))

qqplot(rchisq(n_rep , df = 2), samp_0.05)
qqline(samp_0.05, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)

qqplot(rchisq(n_rep , df = 2), samp_0.1)
qqline(samp_0.1, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)

qqplot(rchisq(n_rep , df = 2), samp_0.15)
qqline(samp_0.15, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)

par(mfrow = c(1,1))

title(main = 'Samples from Cauchy(0,1)')

