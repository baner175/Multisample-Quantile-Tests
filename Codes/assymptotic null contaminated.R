rm(list = ls())
library(ggpubr)

n_samp <- 200
n_rep <- 1.5e4

options('contrasts')

source('assymptotic null functions.R')


#-------------------------------------------------------------
# N(0,1) distribution:

cval <- 50
loc <- 0
scl <- 1

samp_0.05 <- fun_cont_norm(eps = 0.05,
                            loc = loc, scl = scl,
                            n_rep = n_rep, 
                            n_samp = n_samp, cval = cval)$statistic


samp_0.1 <- fun_cont_norm(eps = 0.1,
                           loc = loc, scl = scl,
                           n_rep = n_rep, 
                           n_samp = n_samp, cval = cval)$statistic

samp_0.15 <- fun_cont_norm(eps = 0.15,
                            loc = loc, scl = scl,
                            n_rep = n_rep, 
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

#-------------------------------------------------------------
# t_5 distribution:

cval <- 50
dof <- 5


samp_0.05 <- fun_cont_t(eps = 0.05,
                       deg = dof, n_rep = n_rep, 
                       n_samp = n_samp, cval = cval)$statistic


samp_0.1 <- fun_cont_t(eps = 0.1,
                       deg = dof, n_rep = n_rep, 
                       n_samp = n_samp, cval = cval)$statistic

samp_0.15 <- fun_cont_t(eps = 0.15,
                       deg = dof, n_rep = n_rep, 
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

#-------------------------------------------------------------
# Cauchy(0,1) distribution:

cval <- 50
loc <- 0
scl <- 1

samp_0.05 <- fun_cont_cauchy(eps = 0.05,
                       loc = loc, scl = scl,
                       n_rep = n_rep, 
                       n_samp = n_samp, cval = cval)$statistic


samp_0.1 <- fun_cont_cauchy(eps = 0.1,
                            loc = loc, scl = scl,
                            n_rep = n_rep, 
                            n_samp = n_samp, cval = cval)$statistic

samp_0.15 <- fun_cont_cauchy(eps = 0.15,
                            loc = loc, scl = scl,
                            n_rep = n_rep, 
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

