rm(list = ls())

library(car)
library(quantreg)
library(RVAideMemoire)

source('assymptotic null functions.R')

data = read.table(
  "http://www.stat.umn.edu/~gary/book/fcdae.data/ex6.6",header=TRUE)
#View(data)

plot(y~as.factor(source), data = data)

lmod <- lm(y~as.factor(source), data = data)

par(mfrow = c(2,2))
plot(lmod)
par(mfrow = c(1,1))
ind <- abs(rstandard(lmod))>3

tests <- c('Kruskal-Wallis',
           "Mood's Median",
           'F-test',
           'Quant-med')

p_with <- c()

p_with <- c(p_with,
            kruskal.test(y~as.factor(source), data = data)$p.value)

p_with <- c(p_with,
            mood.medtest(y ~ as.factor(source), data = data)$p.value)
anova_test <- aov(y~as.factor(source), data = data)
p_with <- c(p_with,
            summary(anova_test)[[1]][["Pr(>F)"]][1])

p_with <- c(p_with,
            multiq(y = data$y, g = data$source)$p.value)

#-------------------------------------------------------------------------------

data_rm <- data[!ind,]


p_without <- c()

p_without <- c(p_without,
               kruskal.test(y~as.factor(source), data = data_rm)$p.value)

p_without <- c(p_without,
               mood.medtest(y ~ as.factor(source), data = data_rm)$p.value)
anova_test <- aov(y~as.factor(source), data = data_rm)
p_without <- c(p_without,
               summary(anova_test)[[1]][["Pr(>F)"]][1])

p_without <- c(p_without,
            multiq(y = data_rm$y, g = data_rm$source)$p.value)


#-------------------------------------------------------------------------------

p_compare <- data.frame(round(p_with, digit = 6), 
                        round(p_without, digit = 6))
rownames(p_compare) <- tests
colnames(p_compare) <- c('With outlier', 'Without outlier')
View(p_compare)



