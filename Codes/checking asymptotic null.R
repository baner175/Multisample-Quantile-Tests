rm(list = ls())
n_samp <- 200
n_rep <- 1.5e4

source('assymptotic null functions.R')

set.seed(80531)

plt_norm <- fun_cont_norm(eps = 0, loc = 0, scl = 1, n_rep = n_rep,
                          n_samp = n_samp)

plt_cauchy <- fun_cont_cauchy(eps = 0, loc = 0, scl = 1, n_rep = n_rep,
                              n_samp = n_samp)

plt_t5 <- fun_cont_t(eps = 0, deg = 5, n_rep = n_rep, 
                     n_samp = n_samp)

ggpubr::ggarrange(plt_norm$plot, plt_cauchy$plot, plt_t5$plot, 
                  labels = c('N(0,1)', 'Cauchy (0,1)', 't5'),
                  align = 'hv', hjust = c(-5,-2,-15))

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

