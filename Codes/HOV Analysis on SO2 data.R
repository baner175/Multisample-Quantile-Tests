rm(list = ls())

library(car)
library(quantreg)
data = read.table(
  "http://www.stat.umn.edu/~gary/book/fcdae.data/ex6.6",header=TRUE)
#View(data)

plot(y~as.factor(source), data = data)

lmod <- lm(y~as.factor(source), data = data)

par(mfrow = c(2,2))
plot(lmod)
par(mfrow = c(1,1))
ind <- abs(rstandard(lmod))>3

tests <- c('Bartlett',
           'Fligner-Killeen',
           'Levene-Med',
           'Levene-10%',
           'IQR')


p_with <- c()


p_with <- c(p_with,
            bartlett.test(y ~ as.factor(source), data = data)$p.value)

p_with <- c(p_with, 
            fligner.test(y ~ as.factor(source), data = data)$p.value )

p_with <- c(p_with,
            leveneTest(y ~ as.factor(source), data = data,
            location = median)[3][[1]][1])

p_with <- c(p_with,
            leveneTest(y ~ as.factor(source), data = data,
           location = mean, trim = 0.1)[3][[1]][1])


# IQR test----------------------------------------------------------------------

scl <- 0.3887
qmod_1 <- rq(y~as.factor(source), data = data, tau = 0.25)
qmod_3 <- rq(y~as.factor(source), data = data, tau = 0.75)
n_vec <- table(data$source)
b <- c(coef(qmod_1), coef(qmod_3))
X <- model.matrix(lmod)
Q <- matrix(0, ncol = 4, nrow = 4)

for(i in 1:nrow(X))
{
  Q <- Q+X[i,]%*%t(X[i,])
}
Q <- Q/nrow(X)

p <- c(0.25, 0.75)

f <- qnorm(p, sd = scl) |> dnorm(sd = scl)
f <- 1/f

mid <- cbind(p[1]*(1-p), (1-p[2])*p)

f <- diag(f)

Om <- f%*%mid%*%f

Sig_b <- kronecker(Om, solve(Q))

mat <- rbind(
  c(1,0,0,0),
  c(1,1,0,0),
  c(1,0,1,0),
  c(1,0,0,1)
)
A <- kronecker(rbind(c(1,0), c(0,1)), mat)

B <- rbind(
  c(-1,0,0,0,1,0,0,0),
  c(0,-1,0,0,0,1,0,0),
  c(0,0,-1,0,0,0,1,0),
  c(0,0,0,-1,0,0,0,1)
)

M <- B%*%A

C <- rbind(
  c(-1,1,0,0),
  c(0,-1,1,0),
  c(0,0,-1,1)
)

M <- C%*%M

v <- M%*%b

dim(v) <- c(3,1)

Sig <- M%*%Sig_b%*%t(M)

T_IQR <- nrow(X) * t(v) %*% solve(Sig)%*% v |> as.numeric()
p_IQR <- pchisq(T_IQR, df = 3, lower.tail = FALSE)

p_with <- c(p_with, p_IQR)
#------------------------------------------------------------------------------


data_rm <- data[!ind,]

plot(y~as.factor(source), data = data_rm)

p_without <- c()


p_without <- c(p_without,
               bartlett.test(y ~ as.factor(source), data = data_rm)$p.value)

p_without <- c(p_without, 
               fligner.test(y ~ as.factor(source), data = data_rm)$p.value )

p_without <- c(p_without,
               leveneTest(y ~ as.factor(source), data = data_rm,
                          location = median)[3][[1]][1])

p_without <- c(p_without,
               leveneTest(y ~ as.factor(source), data = data_rm,
                          location = mean, trim = 0.1)[3][[1]][1])


# IQR test----------------------------------------------------------------------

scl <- 0.3887
qmod_1 <- rq(y~as.factor(source), data = data_rm, tau = 0.25)
qmod_3 <- rq(y~as.factor(source), data = data_rm, tau = 0.75)
n_vec <- table(data$source)
b <- c(coef(qmod_1), coef(qmod_3))
X <- model.matrix(lmod)
Q <- matrix(0, ncol = 4, nrow = 4)

for(i in 1:nrow(X))
{
  Q <- Q+X[i,]%*%t(X[i,])
}
Q <- Q/nrow(X)

p <- c(0.25, 0.75)

f <- qnorm(p, sd = scl) |> dnorm(sd = scl)
f <- 1/f

mid <- cbind(p[1]*(1-p), (1-p[2])*p)

f <- diag(f)

Om <- f%*%mid%*%f

Sig_b <- kronecker(Om, solve(Q))

mat <- rbind(
  c(1,0,0,0),
  c(1,1,0,0),
  c(1,0,1,0),
  c(1,0,0,1)
)
A <- kronecker(rbind(c(1,0), c(0,1)), mat)

B <- rbind(
  c(-1,0,0,0,1,0,0,0),
  c(0,-1,0,0,0,1,0,0),
  c(0,0,-1,0,0,0,1,0),
  c(0,0,0,-1,0,0,0,1)
)

M <- B%*%A

C <- rbind(
  c(-1,1,0,0),
  c(0,-1,1,0),
  c(0,0,-1,1)
)

M <- C%*%M

v <- M%*%b

dim(v) <- c(3,1)

Sig <- M%*%Sig_b%*%t(M)

T_IQR <- nrow(X) * t(v) %*% solve(Sig)%*% v |> as.numeric()
p_IQR <- pchisq(T_IQR, df = 3, lower.tail = FALSE)

p_without <- c(p_without, p_IQR)

p_compare <- data.frame(round(p_with, digit = 6), p_without)
rownames(p_compare) <- tests
colnames(p_compare) <- c('With outlier', 'Without outlier')
View(p_compare)
