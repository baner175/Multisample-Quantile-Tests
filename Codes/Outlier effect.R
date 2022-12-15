par(mfrow = c(1,2))

x <- runif(1000, 5, 10)
y <- 3 - 5*x + rnorm(1000, 5)

plot(y~x)

require(quantreg)

med_reg1 <- rq(y~x, tau = 0.5)

abline(med_reg1, col = 'red', lwd = 2)

quart1_reg <- rq(y~x, tau = 0.25)
abline(quart1_reg, col = 'red', lty = 2, lwd = 2)

quart3_reg <- rq(y~x, tau = 0.75)
abline(quart3_reg, col = 'red', lty = 3, lwd = 2)


olsreg1 <- lm(y~x)
abline(olsreg1, col = 'blue', lwd = 2)

legend('topright', legend = c('tau = 0.5', 'tau = 0.25', 'tau = 0.75', 'OLS'),
       col = c('red', 'red', 'red', 'blue'),
       lty = c(1,2,3,1), lwd = 2)



x <- runif(1000, 5, 10)
y <- 3 - 5*x + rnorm(1000, 5)

# ind <- which(x>9)[1:10]
# y[ind] <- -20

ind <- which(x<6)[1:10]
y[ind] <- -40

plot(y~x)

require(quantreg)

med_reg2 <- rq(y~x, tau = 0.5)

abline(med_reg1, col = 'red', lwd = 2)

quart1_reg2 <- rq(y~x, tau = 0.25)
abline(quart1_reg, col = 'red', lty = 2, lwd = 2)

quart3_reg2 <- rq(y~x, tau = 0.75)
abline(quart3_reg, col = 'red', lty = 3, lwd = 2)


olsreg2 <- lm(y~x)
abline(olsreg2, col = 'blue', lwd = 2)

legend('topright', legend = c('tau = 0.5', 'tau = 0.25', 'tau = 0.75', 'OLS'),
       col = c('red', 'red', 'red', 'blue'),
       lty = c(1,2,3,1), lwd = 2)