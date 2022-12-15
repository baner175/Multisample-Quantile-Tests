rm(list = ls())

source('assymptotic null functions.R')

n_rep <- 1.5e4
n_samp <- 200

plt_norm <- IQR_cont_norm(eps = 0, loc = c(1,2,3),
                          scl = 1, n_rep = n_rep,
                          n_samp = n_samp, cval = 0)

plt_cauchy <- IQR_cont_cauchy(eps = 0, loc = c(1,2,3),
                              scl = 1, n_rep = n_rep, 
                              n_samp = n_samp, cval = 0)

plt_t5 <- IQR_cont_t(eps = 0, n_rep = n_rep, n_samp = n_samp, deg = 5)

ggpubr::ggarrange(plt_norm$plot, plt_cauchy$plot, plt_t5$plot, 
                  labels = c('Normal', 'Cauchy', 't5'),
                  align = 'hv', hjust = c(-3,-3,-15))

par(mfrow = c(2,2))
qqplot(plt_norm$statistic, rchisq(1e4, df = 2), main = 'Normal',
       xlab = 'smaple quantiles', 
       ylab = 'theoretical quantiles')
qqline(plt_norm$statistic, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)


qqplot(plt_cauchy$statistic, rchisq(1e4, df = 2), main = 'Cauchy',
       xlab = 'smaple quantiles', 
       ylab = 'theoretical quantiles')
qqline(plt_cauchy$statistic, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)



qqplot(plt_t5$statistic, rchisq(1e4, df = 2), main = 't',
       xlab = 'smaple quantiles', 
       ylab = 'theoretical quantiles')
qqline(plt_t5$statistic, distribution = function(p) qchisq(p, df = 2),
       probs = c(0.1, 0.6), col = 2)



par(mfrow = c(1,1))


