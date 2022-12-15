rm(list = ls())

library(doParallel)

registerDoParallel(cores = 4)

lev_cont_norm <- function(eps = 0, alpha = 0.05,
                          loc = 0, scl = 1, n_rep, n_samp,
                          cval)
{
  set.seed(1008)
  library(doParallel)
  
  decision <- foreach(i = 1:n_rep, 
                      .combine = 'rbind',
                      .packages = c('quantreg', 
                                    'RVAideMemoire')) %dopar%
    {
      y1 <- rnorm(n = n_samp, mean = loc, sd = scl)
      y2 <- rnorm(n = n_samp, mean = loc, sd = scl)
      if(eps>0){
        ind <- sample(1:n_samp, size = eps*n_samp)
        y2[ind] <- cval
      }
      y3 <- rnorm(n = n_samp, mean = loc, sd = scl)
      
      y <- c(y1,y2,y3)
      fac <- as.factor(rep(1:3, each = n_samp))
      
      df <- data.frame(y, fac)
      quant_mod <- rq(y~fac, data = df, tau = 0.5)
      sm_q <- summary(quant_mod, se = 'iid', covariance = TRUE)
      
      A <- rbind(c(1,0,0),c(1,1,0), c(1,0,1))
      B <- rbind(c(-1,1,0),
                 c(0,-1,1))
      M <- B%*%A
      
      beta_vec <- M%*%coef(quant_mod)
      
      Sig <- M%*%sm_q$cov%*%t(M)
      
      st <- t(beta_vec)%*%solve(Sig)%*%beta_vec |> as.numeric()
      dec_T <- (st > qchisq(0.95, df = 2))
      
      med_dec <- mood.medtest(y, fac)$p.value<0.05
      
      kw_dec <- kruskal.test(x = y, g = fac)$p.value < 0.05
      
      anova_test <- aov(y~fac)
      anova_dec <- summary(anova_test)[[1]][["Pr(>F)"]][1] < 0.05
      
      c(dec_T, med_dec, kw_dec, anova_dec)
      
    }
  colnames(decision) <- c('Q-test', 'Median Test', 
                          'KW Test', 'F-test')
  return(apply(decision, 2, mean) |> round(digits = 3))
}

n_samp <- 100
n_rep <- 1e4
loc <- 0
scl <- 1
cval <- 50
eps <- seq(0,0.4, by = 0.02)

type_1_err <- matrix(NA_real_, nrow = length(eps), ncol = 4)

for (i in 1:length(eps)) {
  msg <- sprintf('Iteration: %d / %d', i,
                 length(eps))
  message(msg)
  type_1_err[i,] <- lev_cont_norm(eps = eps[i],
                                  n_rep = n_rep,
                                  n_samp = n_samp,
                                  cval = cval)
  
}


colnames(type_1_err) <- c('Q-test', 'Median Test', 
                        'KW Test', 'F-test')

type_1_err <- cbind(type_1_err, eps)

write.csv(type_1_err, 'type_I_error_cont_Norm.csv')

type_1_err <- read.csv('type_I_error_cont_Norm.csv', header = T)

type_1_err <- type_1_err[,-1]

View(type_1_err)
matplot(x = type_1_err[,5],
        y = type_1_err[,1:4])

matplot(x = type_1_err[,5],
        y = type_1_err[,1:4],
        type = 'l', lwd = 2,
        lty = 1, xlab = 'contamination rate',
        ylab = 'type I error',
        col = 1:4
        )
legend('bottomright',
       col = 1:4, lty = 1, lwd = 2,
       legend = colnames(type_1_err)[1:4],
       bty = 'n')
abline(h = 0.05, col = 'red', lty = 2)
