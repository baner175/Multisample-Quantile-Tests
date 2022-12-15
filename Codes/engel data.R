data(engel)

plot(foodexp~income, data = engel)


require(quantreg)

t <- c(0.5, seq(0.05, 0.9, length.out = 10))

for(j in t)
{
  mod <- rq(foodexp~income, tau = j, data = engel)
  if(j == 0.5)
  {
    abline(mod, lty = 1, col = 'grey', lwd = 2)
  }else
  {
    abline(mod, lty = 2, col = 'grey', lwd = 1)
  }
  
}

olsmod <- lm(foodexp~income, data = engel)
abline(olsmod, lty = 1, lwd = 2)