# Functions for checking the null distribution of the test 
# statistic under null with/without contamination

fun_cont_norm <- function(eps, loc = 0, scl = 1, n_rep, 
                          n_samp, cval = 0)
{
  library(quantreg)
  library(ggplot2)
  st <- NULL
  for(i in 1:n_rep)
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
    
    st <- c(st, 
            t(beta_vec)%*%solve(Sig)%*%beta_vec
            |> as.numeric())
  }
  dat <- data.frame(st)
  plt <- ggplot(data = dat, mapping = aes(x = st)) + 
    geom_histogram(mapping = aes(y = ..density..),
                   color = 'black', 
                   fill = 'white') + geom_density(color = 'blue') + 
    stat_function(fun = 'dchisq',
                  args = list(df = 2)) + 
    xlab(label = 'T') + # using unicode for \epsilon
    ggtitle(label = paste('\u03F5', '=',100*eps,'%'))
  return(list('statistic' = st,
              'plot' = plt))
}


fun_cont_t <- function(eps, loc = 0, deg = 5, n_rep, 
                       n_samp, cval = 0)
{
  library(quantreg)
  library(ggplot2)
  st <- NULL
  for(i in 1:n_rep)
  {
    y1 <- loc + rt(n = n_samp, df = deg)
    y2 <- loc + rt(n = n_samp, df = deg)
    if(eps>0){
      ind <- sample(1:n_samp, size = eps*n_samp)
      y2[ind] <- cval
    }
    y3 <- loc + rt(n = n_samp, df = deg)
    
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
    
    st <- c(st, 
            t(beta_vec)%*%solve(Sig)%*%beta_vec
            |> as.numeric())
  }
  dat <- data.frame(st)
  plt <- ggplot(data = dat, mapping = aes(x = st)) + 
    geom_histogram(mapping = aes(y = ..density..),
                   color = 'black', 
                   fill = 'white') + geom_density(color = 'blue') + 
    stat_function(fun = 'dchisq',
                  args = list(df = 2))+ 
    xlab(label = 'T') + # using unicode for \epsilon
    ggtitle(label = paste('\u03F5', '=',100*eps,'%'))
  return(list('statistic' = st,
              'plot' = plt))
}

fun_cont_cauchy <- function(eps, loc = 0, scl = 1, n_rep, 
                            n_samp, cval = 0)
{
  library(quantreg)
  library(ggplot2)
  st <- NULL
  for(i in 1:n_rep)
  {
    y1 <- rcauchy(n = n_samp, location = loc, scale = scl)
    y2 <- rcauchy(n = n_samp, location = loc, scale = scl)
    if(eps>0){
      ind <- sample(1:n_samp, size = eps*n_samp)
      y2[ind] <- cval
    }
    y3 <- rcauchy(n = n_samp, location = loc, scale = scl)
    
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
    
    st <- c(st, 
            t(beta_vec)%*%solve(Sig)%*%beta_vec
            |> as.numeric())
  }
  dat <- data.frame(st)
  plt <- ggplot(data = dat, mapping = aes(x = st)) + 
    geom_histogram(mapping = aes(y = ..density..),
                   color = 'black', 
                   fill = 'white') + geom_density(color = 'blue') + 
    stat_function(fun = 'dchisq',
                  args = list(df = 2))+ 
    xlab(label = 'T') + # using unicode for \epsilon
    ggtitle(label = paste('\u03F5', '=',100*eps,'%'))
  return(list('statistic' = st,
              'plot' = plt))
}


IQR_cont_norm <- function(eps, loc = c(0,0,0), scl = 1, n_rep, 
                          n_samp, cval = 0)
{
  library(quantreg)
  library(ggplot2)
  N <- 3*n_samp
  T_stat <- c()
  for(i in 1:n_rep)
  {
    y1 <- loc[1] + rnorm(n_samp, mean = 0,
                         sd = scl)
    y2 <- loc[2] + rnorm(n_samp, mean = 0,
                         sd = scl)
    if(eps>0){
      ind <- sample(1:n_samp, size = eps*n_samp)
      y2[ind] <- cval
    }
    y3 <- loc[3] + rnorm(n_samp, mean = 0,
                         sd = scl)
    fac <- rep(1:3, each = n_samp) |> as.factor()
    
    df <- data.frame(y = c(y1, y2, y3), fac = fac)
    
    qmod_1 <- rq(y~fac, data = df, tau = 0.25)
    qmod_3 <- rq(y~fac, data = df, tau = 0.75)
    
    b <- c(coef(qmod_1), coef(qmod_3))
    
    mat <- rbind(
      c(3*n_samp, n_samp, n_samp),
      c(n_samp, n_samp, 0),
      c(n_samp, 0, n_samp)
    )
    
    Q <- mat/N
    
    p <- c(0.25, 0.75)
    
    f <- qnorm(p, sd = scl) |> dnorm(sd = scl)
    f <- 1/f
    
    mid <- cbind(p[1]*(1-p), (1-p[2])*p)
    
    f <- diag(f)
    
    Om <- f%*%mid%*%f
    
    Sig_b <- kronecker(Om, solve(Q))
    
    mat <- rbind(
      c(1,0,0),
      c(1,1,0),
      c(1,0,1)
    )
    A <- kronecker(rbind(c(1,0), c(0,1)), mat)
    
    B <- rbind(
      c(-1,0,0,1,0,0),
      c(0,-1,0,0,1,0),
      c(0,0,-1,0,0,1)
    )
    
    M <- B%*%A
    
    C <- rbind(
      c(-1,1,0),
      c(0,-1,1)
    )
    
    M <- C%*%M
    
    v <- M%*%b
    
    dim(v) <- c(2,1)
    
    Sig <- M%*%Sig_b%*%t(M)
    
    T_stat <- c(T_stat,
                3*n_samp * t(v) %*% solve(Sig)%*% v |> as.numeric())
    
  }
  dat <- data.frame(T_stat)
  plt <- ggplot(data = dat, mapping = aes(x = T_stat)) + 
    geom_histogram(mapping = aes(y = ..density..),
                   color = 'black', 
                   fill = 'white') + geom_density(color = 'blue') + 
    stat_function(fun = 'dchisq',
                  args = list(df = 2)) + 
    xlab(label = bquote(T[IQR])) + # using unicode for \epsilon
    ggtitle(label = paste('\u03F5', '=',100*eps,'%'))
  
  return(list('plot' = plt,
              'statistic' = T_stat))
  
}

IQR_cont_cauchy <- function(eps, loc = c(0,0,0), scl = 1, n_rep, 
                          n_samp, cval = 0)
{
  library(quantreg)
  library(ggplot2)
  T_stat <- c()
  N <- 3*n_samp
  for(i in 1:n_rep)
  {
    y1 <- loc[1] + rcauchy(n_samp, location = 0,
                         scale = scl)
    y2 <- loc[2] + rcauchy(n_samp, location = 0,
                           scale = scl)
    if(eps>0){
      ind <- sample(1:n_samp, size = eps*n_samp)
      y2[ind] <- cval
    }
    y3 <- loc[3] + rcauchy(n_samp, location = 0,
                           scale = scl)
    fac <- rep(1:3, each = n_samp) |> as.factor()
    
    df <- data.frame(y = c(y1, y2, y3), fac = fac)
    
    qmod_1 <- rq(y~fac, data = df, tau = 0.25)
    qmod_3 <- rq(y~fac, data = df, tau = 0.75)
    
    b <- c(coef(qmod_1), coef(qmod_3))
    
    mat <- rbind(
      c(3*n_samp, n_samp, n_samp),
      c(n_samp, n_samp, 0),
      c(n_samp, 0, n_samp)
    )
    
    Q <- mat/N
    
    p <- c(0.25, 0.75)
    
    f <- qcauchy(p, scale = scl) |> dcauchy(scale = scl)
    f <- 1/f
    
    mid <- cbind(p[1]*(1-p), (1-p[2])*p)
    
    f <- diag(f)
    
    Om <- f%*%mid%*%f
    
    Sig_b <- kronecker(Om, solve(Q))
    
    mat <- rbind(
      c(1,0,0),
      c(1,1,0),
      c(1,0,1)
    )
    A <- kronecker(rbind(c(1,0), c(0,1)), mat)
    
    B <- rbind(
      c(-1,0,0,1,0,0),
      c(0,-1,0,0,1,0),
      c(0,0,-1,0,0,1)
    )
    
    M <- B%*%A
    
    C <- rbind(
      c(-1,1,0),
      c(0,-1,1)
    )
    
    M <- C%*%M
    
    v <- M%*%b
    
    dim(v) <- c(2,1)
    
    Sig <- M%*%Sig_b%*%t(M)
    
    T_stat <- c(T_stat,
                3*n_samp * t(v) %*% solve(Sig)%*% v |> as.numeric())
    
  }
  dat <- data.frame(T_stat)
  plt <- ggplot(data = dat, mapping = aes(x = T_stat)) + 
    geom_histogram(mapping = aes(y = ..density..),
                   color = 'black', 
                   fill = 'white') + geom_density(color = 'blue') + 
    stat_function(fun = 'dchisq',
                  args = list(df = 2)) + 
    xlab(label = bquote(T[IQR])) + # using unicode for \epsilon
    ggtitle(label = paste('\u03F5', '=',100*eps,'%'))
  
  return(list('plot' = plt,
              'statistic' = T_stat))
}

IQR_cont_t <- function(eps, loc = c(0,0,0), deg = 5, n_rep, 
                       n_samp, cval = 0)
{
  library(quantreg)
  library(ggplot2)
  T_stat <- c()
  N <- 3*n_samp
  for(i in 1:n_rep)
  {
    y1 <- loc[1] + rt(n_samp, df = deg)
    y2 <- loc[2] + rt(n_samp, df = deg)
    if(eps>0){
      ind <- sample(1:n_samp, size = eps*n_samp)
      y2[ind] <- cval
    }
    y3 <- loc[3] + rt(n_samp, df = deg)
    fac <- rep(1:3, each = n_samp) |> as.factor()
    
    df <- data.frame(y = c(y1, y2, y3), fac = fac)
    
    qmod_1 <- rq(y~fac, data = df, tau = 0.25)
    qmod_3 <- rq(y~fac, data = df, tau = 0.75)
    
    b <- c(coef(qmod_1), coef(qmod_3))
    
    mat <- rbind(
      c(3*n_samp, n_samp, n_samp),
      c(n_samp, n_samp, 0),
      c(n_samp, 0, n_samp)
    )
    
    Q <- mat/N
    
    p <- c(0.25, 0.75)
    
    f <- qt(p, df = deg) |> dt(df = deg)
    f <- 1/f
    
    mid <- cbind(p[1]*(1-p), (1-p[2])*p)
    
    f <- diag(f)
    
    Om <- f%*%mid%*%f
    
    Sig_b <- kronecker(Om, solve(Q))
    
    mat <- rbind(
      c(1,0,0),
      c(1,1,0),
      c(1,0,1)
    )
    A <- kronecker(rbind(c(1,0), c(0,1)), mat)
    
    B <- rbind(
      c(-1,0,0,1,0,0),
      c(0,-1,0,0,1,0),
      c(0,0,-1,0,0,1)
    )
    
    M <- B%*%A
    
    C <- rbind(
      c(-1,1,0),
      c(0,-1,1)
    )
    
    M <- C%*%M
    
    v <- M%*%b
    
    dim(v) <- c(2,1)
    
    Sig <- M%*%Sig_b%*%t(M)
    
    T_stat <- c(T_stat,
                3*n_samp * t(v) %*% solve(Sig)%*% v |> as.numeric())
    
  }
  dat <- data.frame(T_stat)
  plt <- ggplot(data = dat, mapping = aes(x = T_stat)) + 
    geom_histogram(mapping = aes(y = ..density..),
                   color = 'black', 
                   fill = 'white') + geom_density(color = 'blue') + 
    stat_function(fun = 'dchisq',
                  args = list(df = 2))+  
    xlab(label = bquote(T[IQR])) + # using unicode for \epsilon
    ggtitle(label = paste('\u03F5', '=',100*eps,'%'))
  
  return(list('plot' = plt,
              'statistic' = T_stat))
}


multiq <- function(y, g, q = 0.5)
{
  require(quantreg)
  qmod <- rq(y ~ as.factor(g), tau = q)
  sm_q <- summary(qmod, se = 'iid', covariance = TRUE)
  nfac <- length(unique(g))
  A <- diag(1, nfac)
  A[,1] <- rep(1, nfac)
  v <- rep(0,nfac)
  v[1:2] <- c(-1,1)
  B <- toeplitz(v)
  B[lower.tri(B)] <- 0
  B <- B[-nfac,]
  M <- B%*%A
  
  beta_vec <- M%*%coef(qmod)
  
  Sig <- M%*%sm_q$cov%*%t(M)
  
  st <- t(beta_vec)%*%solve(Sig)%*%beta_vec |> as.numeric()
  
  p_val <- pchisq(st, df = nfac-1, lower = FALSE)
  
  return(list('statistic' = st, 'p.value' = p_val))
}
