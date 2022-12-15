set.seed(8913)

x <- runif(1000, 5, 10)
y <- 3 - 30*x + (x-4)^1.2 *rnorm(1000, sd = 3)

plot(y~x)

require(quantreg)

med_reg1 <- rq(y~x, tau = 0.5)

t <- c(0.5, seq(0.05, 0.9, length.out = 10))

for(j in t)
{
  mod <- rq(y~x, tau = j)
  if(j == 0.5)
  {
    abline(mod, lty = 1, col = 'blue', lwd = 2)
  }else
  {
    abline(mod, lty = 2, col = 'blue', lwd = 1)
  }
  
}


olsreg1 <- lm(y~x)
abline(olsreg1, col = 'red', lwd = 2)
