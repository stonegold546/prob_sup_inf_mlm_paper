source("ps_functions_new.R")

# Hypothetical Ruscio

x <- c(6, 7, 8, 7, 9, 6, 5, 4, 7, 8, 7, 6, 9, 5, 4)
y <- c(4, 3, 5, 3, 6, 2, 2, 1, 6, 7, 4, 3, 2, 4, 3)

ps(x, y)
# [1] 0.8844444
inf.ps.boot(data.frame(y = c(x, y), x = c(rep(1, length(x)), rep(0, length(y)))))
# [1] 0.7070636 0.9642857
coef(summary(glm(ps(x, y) ~ 1, binomial, weights = length(x) + length(y))))
#             Estimate Std. Error  z value    Pr(>|z|)
# (Intercept) 2.035208  0.5710934 3.563705 0.000365657
inf.ps.logit(length(x) + length(y), ps(x, y))
# [1] 0.7142031 0.9590869

rm(x, y)
gc()

# CRT

dat <- read.csv("ToolsK_archive.csv")
# outcome must be named y
# grouping variable must be named x and coded 0-1
# cluster ID must be named ID going from 1 to number of clusters
dat
dat$y <- dat$ap_ws.G1
table(dat$x <- dat$Condition)
dat$ID <- dat$SchoolID
# remove missing data
dat <- na.omit(dat[, c("ID", "x", "y")])
head(dat)
table(dat$ID <- as.integer(factor(dat$ID)))
aggregate(ID ~ x, dat, function (x) length(unique(x)))
#   x ID
# 1 0 14
# 2 1 17
mean(table(dat$ID))
# [1] 21.32258 per sch
min(table(dat$ID))
# [1] 9 students min

plot(density(dat$y[dat$x == 1]))
lines(density(dat$y[dat$x == 0]), col = "red")

summary(lmer(y ~ (1 | ID), dat))
# Groups   Name        Variance Std.Dev.
# ID       (Intercept)  54.28    7.368  
# Residual             290.55   17.045  
# Number of obs: 661, groups:  ID, 31
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  457.671      1.505   304.1

54.28 / (54.28 + 290.55)
# [1] 0.1574109

(vargha <- ps(dat$y[dat$x == 1], dat$y[dat$x == 0]))
# [1] 0.5363497

# BCa
system.time(vargha.bca <- inf.ps.boot(dat, int = .95))
#  user  system elapsed 
# 3.603   0.087   3.691 
vargha.bca
# [1] 0.4945401 0.5821241

# recommended approach
round(ps.logit.crve <- inf.ps.logit.crve(dat, int = .95), 3)
# [1] 0.539 0.432 0.642

# placement scores mlm approach
round(ps.pl.cl <- inf.ps.pl.mlm(dat, int = .95), 3)
# [1] 0.538 0.434 0.638


# HSB

hsb <- merTools::hsb

# outcome must be named y
# grouping variable must be named x and coded 0-1
# cluster ID must be named ID going from 1 to number of clusters
dat <- na.omit(data.frame(y = hsb$mathach, x = hsb$schtype, ID = hsb$schid))
table(dat$ID <- as.integer(factor(dat$ID)))
# 160 schools
aggregate(ID ~ x, dat, function (x) length(unique(x)))
# 1 0 90
# 2 1 70
mean(table(dat$ID))
# [1] 44.90625 per sch
min(table(dat$ID))
# [1] 14 students min

plot(density(dat$y[dat$x == 1]))
lines(density(dat$y[dat$x == 0]), col = "red")

summary(lmer(y ~ (1 | ID), dat))
# Random effects:
#   Groups   Name        Variance Std.Dev.
# ID       (Intercept)  8.614   2.935   
# Residual             39.148   6.257   
# Number of obs: 7185, groups:  ID, 160
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  12.6370     0.2444   51.71

# ICC:
8.614 / (8.614 + 39.148)
# [1] 0.1803526

(vargha <- ps(dat$y[dat$x == 1], dat$y[dat$x == 0]))
# [1] 0.6158195

# recommended approach
round(ps.logit.crve <- inf.ps.logit.crve(dat, int = .95), 3)
# [1] 0.615 0.580 0.649

# placement scores mlm approach
round(ps.pl.cl <- inf.ps.pl.mlm(dat, int = .95), 3)
# [1] 0.615 0.579 0.650


# Not part of example ----

library(fitdistrplus)
# One could use these estimates in a simulation
fitdist((dat$y[dat$x == 0] + 5) / 30, dbeta)
# shape1 1.811231 0.04073575
# shape2 1.445011 0.03146809
fitdist((dat$y[dat$x == 1] + 5) / 30, dbeta)
# shape1 2.510751 0.05895978
# shape2 1.396141 0.03041828


plgtnrm <- function (q, mu, scale) pnorm(qlogis(q), mu, scale)
dlgtnrm <- function (x, mu, scale) {
  dnorm(qlogis(x), mu, scale) * (1 / (x * (1 - x)))
}

fitdist((dat$y[dat$x == 0] + 5) / 30, dlgtnrm, start = list(mu = 0, scale = 1))
# mu    0.3094311 0.02276241
# scale 1.3736886 0.01609479
fitdist((dat$y[dat$x == 1] + 5) / 30, dlgtnrm, start = list(mu = 0, scale = 1))
# mu    0.773806 0.02182698
# scale 1.299210 0.01543468


