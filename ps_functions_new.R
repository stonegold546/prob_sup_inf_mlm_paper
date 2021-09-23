library(sandwich)
library(lmtest)
library(boot)
library(lmerTest)

ps <- function(a, b) {
  # Vargha's A
  mean(outer(a, b, ">")) + .5 * mean(outer(a, b, "=="))
}

inf.ps.boot <- function(d, nboot = 2e3, int = .95) {
  boot.fx <- function(d, idx) {
    d <- d[idx, ]
    ps(d$y[d$x == 1], d$y[d$x == 0])
  }
  unname(boot::boot.ci(boot::boot(
    d, boot.fx, nboot), conf = int, type = "bca")$bca[1, 4:5])
}

inf.ps.logit <- function(n, ps, int = .95) {
  pass <- n * ps
  params <- coef(summary(glm(cbind(pass, n - pass) ~ 1, binomial)))[1, 1:2]
  unname(plogis(
    qnorm(c((1 - int) / 2, int + (1 - int) / 2), params[1], params[2])))
}

inf.ps.logit.crve <- function(d, int = .95) {
  id.0s <- unique(d$ID[d$x == 0])
  id.1s <- unique(d$ID[d$x == 1])
  n.id <- length(unique(d$ID))
  ps <- sapply(id.1s, function (id1) {
    sapply(id.0s, function (id0) {
      ps(d$y[d$ID == id1], d$y[d$ID == id0])
    })
  })
  ns <- sapply(id.1s, function (id1) {
    sapply(id.0s, function (id0) {
      sum(d$ID == id1) + sum(d$ID == id0)
    })
  })
  d.long <- data.frame(
    ps = as.vector(ps), n = as.vector(ns),
    id1.v = rep.int(1:ncol(ns), rep(nrow(ns), ncol(ns))),
    id0.v = rep.int(1:nrow(ns), ncol(ns)))
  d.long$pass <- d.long$n * d.long$ps
  params <- unname(lmtest::coeftest(
    glm(cbind(pass, n - pass) ~ 1, quasibinomial, d.long), vcov. = sandwich::vcovCL,
    type = "HC3", cluster = d.long[, c("id1.v", "id0.v")])[, 1:2])
  plogis(
    params[1] + qt(c(.5, (1 - int) / 2, int + (1 - int) / 2), n.id - 2) * params[2])
}

inf.ps.pl.mlm <- function (d, int = .95) {
  d$pl <- NA
  y0 <- d$y[d$x == 0]
  y1 <- d$y[d$x == 1]
  n.id <- length(unique(d$ID))
  # See placement example from table 2 in 10.1016/S1076-6332(97)80161-4
  # Nut keep comparison consistent across groups
  d$pl[d$x == 1] <- sapply(y1, function (x) {
    sum(x > y0) + .5 * sum(x == y0)
  }) / length(y0)
  d$pl[d$x == 0] <- sapply(y0, function (x) {
    sum(x > y1) + .5 * sum(x == y1)
  }) / length(y1)
  fit.pl <- lmerTest::lmer(pl ~ x + (1 | ID), d)
  params <- coef(summary(fit.pl))["x", 1:3]
  est <- unname((params[1] + 1) / 2)
  se <- params[2]
  lims <- 2 * asinh(qt(c((1 - int) / 2, int + (1 - int) / 2), n.id - 2) / 2 * se / (est * (1 - est)))
  c(est, plogis(qlogis(est) + lims))
}

