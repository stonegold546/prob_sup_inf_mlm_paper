# Install required packages
list.of.packages <- c("sandwich", "sandwich", "lmtest", "boot", "lmerTest", "SimDesign")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
rm(list.of.packages, new.packages)

source("ps_functions_new.R")

library(SimDesign)

calc_pop_ps <- function (mu1, mu2, sg1, sg2) {
  mu1 <- qlogis(mu1)
  mu2 <- qlogis(mu2)
  integrate(function (x) {
    dnorm(qlogis(x), mu1, sg1) * 1 / (x * (1 - x)) *
      pnorm(qlogis(x), mu2, sg2, lower.tail = FALSE)
  }, 0, 1)$value
}

curve(dnorm(qlogis(x), qlogis(.8), 1) * 1 / (x * (1 - x)), n = 1e4, col = "red")
curve(dnorm(qlogis(x), qlogis(.5), .75) * 1 / (x * (1 - x)), n = 1e4, add = TRUE)

calc_pop_ps(.5, .8, .75, 1)
pnorm((qlogis(.8) - qlogis(.5)) / sqrt(.75 ^ 2 + 1))
(Design <- rbind(
  expand.grid(
    n1 = 10, p = c(1 / 3, .5), icc = c(.05, .2), n2 = 3e1,
    mu1 = .7, mu2 = c(.7, .8), sg1 = 1, sg2 = c(1, 1.5)),
  expand.grid(
    n1 = 10, p = c(1 / 3, .5), icc = c(.05, .2), n2 = 3e1,
    mu1 = .5, mu2 = .8, sg1 = .75, sg2 = 1),
  expand.grid(
    n1 = 30, p = 1 / 3, icc = .2, n2 = 3e1,
    mu1 = .5, mu2 = .8, sg1 = .75, sg2 = 1),
  expand.grid(
    n1 = 10, p = 1 / 3, icc = .2, n2 = 5e1,
    mu1 = .5, mu2 = .8, sg1 = .75, sg2 = 1)))
Design <- Design[c(1:4, 9:nrow(Design)), ]
Design
# (condition <- list(
#   pop = pop, n1 = 100, p = 1 / 3, icc = .2, n2 = 100,
#   mu1 = qlogis(.7), mu2 = qlogis(.7), sg1 = 1.25, sg2 = 1))

#-------------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  dat <- data.frame(ID = 1:condition$n2, x = rbinom(condition$n2, 1, condition$p))
  b.var <- c(condition$sg1, condition$sg2) ^ 2 * condition$icc
  w.var <- c(condition$sg1, condition$sg2) ^ 2 * (1 - condition$icc)
  dat$g <- rnorm(
    condition$n2, qlogis(c(condition$mu1, condition$mu2))[dat$x + 1],
    sqrt(b.var[dat$x + 1]))
  dat <- dat[rep(1:condition$n2, round(rnorm(condition$n2, condition$n1))), ]
  dat$y <- plogis(rnorm(nrow(dat), dat$g, sqrt(w.var)[dat$x + 1]))
  dat
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
  pop <- calc_pop_ps(condition$mu1, condition$mu2, condition$sg1, condition$sg2)

  logit.m.res <- inf.ps.logit.crve(dat)
  logit.m <- logit.m.res[1]
  logit.m.ci <- logit.m.res[2:3]
  logit.m.ecr <- logit.m.ci[1] < pop & logit.m.ci[2] > pop
  logit.m.mse <- (logit.m - pop) ^ 2

  pl.m.res <- inf.ps.pl.mlm(dat)
  pl.m <- pl.m.res[1]
  pl.m.ci <- pl.m.res[2:3]
  pl.m.ecr <- pl.m.ci[1] < pop & pl.m.ci[2] > pop
  pl.m.mse <- (pl.m - pop) ^ 2

  c(logit.m = logit.m, logit.m.ecr = logit.m.ecr, logit.m.mse = logit.m.mse,
    pl.m = pl.m, pl.m.ecr = pl.m.ecr, pl.m.mse = pl.m.mse)
}

Summarise <- function(condition, results, fixed_objects = NULL) {
  colMeans(results)
}

#-------------------------------------------------------------------

unlink("SimDesign-results_pop-os_1", recursive = TRUE)
results <- runSimulation(
  design = Design, replications = 4e3, generate = Generate,
  analyse = Analyse, summarise = Summarise, seed = rep(123, nrow(Design)),
  save = TRUE, save_results = TRUE, progress = TRUE,
  filename = "./simdata_1", packages = c("sandwich", "lmtest", "lmerTest"))

#-------------------------------------------------------------------

list.of.packages <- c("ggplot2", "ggforce", "scales", "patchwork", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
rm(list.of.packages, new.packages)

library(ggplot2)
library(ggforce)
library(scales)
library(patchwork)
library(data.table)
theme_set(theme_bw())

calc_pop_ps <- function (mu1, mu2, sg1, sg2) {
  mu1 <- qlogis(mu1)
  mu2 <- qlogis(mu2)
  integrate(function (x) {
    dnorm(qlogis(x), mu1, sg1) * 1 / (x * (1 - x)) *
      pnorm(qlogis(x), mu2, sg2, lower.tail = FALSE)
  }, 0, 1)$value
}

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

(Design <- rbind(
  expand.grid(
    n1 = 10, p = c(1 / 3, .5), icc = c(.05, .2), n2 = 3e1,
    mu1 = .7, mu2 = c(.7, .8), sg1 = 1, sg2 = c(1, 1.5)),
  expand.grid(
    n1 = 10, p = c(1 / 3, .5), icc = c(.05, .2), n2 = 3e1,
    mu1 = .5, mu2 = .8, sg1 = .75, sg2 = 1),
  expand.grid(
    n1 = 30, p = 1 / 3, icc = .2, n2 = 3e1,
    mu1 = .5, mu2 = .8, sg1 = .75, sg2 = 1),
  expand.grid(
    n1 = 10, p = 1 / 3, icc = .2, n2 = 5e1,
    mu1 = .5, mu2 = .8, sg1 = .75, sg2 = 1)))
Design <- Design[c(1:4, 9:nrow(Design)), ]
Design
Design <- Design[1:16, ]
Design
Design$design <- 1:nrow(Design)
Design$pop <- sapply(
  1:nrow(Design), function (i) calc_pop_ps(
    Design$mu1[i], Design$mu2[i], Design$sg1[i], Design$sg2[i]))
Design

raw <- rbindlist(lapply(list.files("SimDesign-results_pop-os_1/"), function (x) {
  res <- as.data.frame(readRDS(paste0("SimDesign-results_pop-os_1/", x))$results)
  res$design <- x
  res
}))

raw.l <- reshape(
  raw, idvar = "ID", varying = list(c(1, 4), c(2, 5), c(3, 6)),
  v.names = c("estimate", "coverage", "mse"),
  direction = "long")
raw.l[, time.t := c("Fractional\nw. CRVE", "Placement\nw. MLM")[time]]
raw.l

raw.l[, design := as.integer(gsub("[A-Z]+|[a-z]+|\\.|\\-", "", design))]
raw.l <- raw.l[order(design)]
raw.l <- raw.l[design <= 16]
raw.l

sum.dat.0 <- raw.l[
  , .(mean = mean(estimate), q05 = quantile(estimate, .05),
      q25 = quantile(estimate, .25), q75 = quantile(estimate, .75),
      q95 = quantile(estimate, .95), mse = mean(mse)),
  list(design, time, time.t)]
sum.dat.0 <- merge(sum.dat.0, Design)
sum.dat.0

sum.dat.0[, Data := ((design - 1) %/% 4) + 1]
sum.dat.0 <- sum.dat.0[order(time)]
sum.dat.0

dodge <- position_dodge(width = .75)
A <- ggplot(sum.dat.0, aes(design, mean - pop, col = time.t)) +
  geom_point(aes(shape = time.t), position = dodge, size = 2.5) +
  geom_linerange(aes(ymin = q05 - pop, ymax = q95 - pop), position = dodge) +
  # geom_segment(aes(y = pop, yend = pop, x = design - .5, xend = design + .5), linetype = 2, col = 1) +
  geom_hline(yintercept = 0, linetype = 2, col = 1) +
  geom_vline(xintercept = seq(4.5, 16, 4), alpha = .5) +
  scale_x_continuous(breaks = 1:nrow(Design), trans = reverse_trans(), limits = c(nrow(Design) + .5, .5),
                     labels = paste0("Data #", sum.dat.0$Data,
                                     ", ICC:", percent(sum.dat.0$icc),
                                     ", P:", percent(sum.dat.0$p))[1:16]) +
  scale_y_continuous(labels = percent_format(1)) +
  scale_color_manual(values = cbbPalette[-1]) +
  scale_shape_manual(values = c(1, 2, 4)) +
  guides(col = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  labs(x = "Condition", y = "Distribution of (PS - Population parameter)",
       tag = "A", col = "", shape = "") +
  theme(legend.position = "bottom",
        panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank()) +
  coord_flip(xlim = c(15.75, 1.25))
A

sum.dat.1 <- raw.l[, .(mean = mean(coverage), count = .N), list(design, time, time.t)]
sum.dat.1$exp <- .95
sum.dat.1$ll <- .925
sum.dat.1$ul <- .975
sum.dat.1[, llm := mean - sqrt(mean * (1 - mean) / count) * qnorm(.95)]
sum.dat.1[, ulm := mean + sqrt(mean * (1 - mean) / count) * qnorm(.95)]
sum.dat.1[, llmjf := qbeta(.05, mean * count + .5, (1 - mean) * count + .5)]
sum.dat.1[, ulmjf := qbeta(.95, mean * count + .5, (1 - mean) * count + .5)]
sum.dat.1
sum.dat.1 <- merge(sum.dat.1, Design)
sum.dat.1
sum.dat.1 <- sum.dat.1[order(time)]
sum.dat.1

B <- ggplot(sum.dat.1, aes(design, mean, col = time.t)) +
  geom_point(aes(shape = time.t), size = 1, position = dodge) +
  geom_linerange(aes(ymin = llmjf, ymax = ulmjf), position = dodge) +
  geom_hline(aes(yintercept = exp), linetype = 2, alpha = .5) +
  geom_hline(aes(yintercept = ll), linetype = 2, col = "#CC6666") +
  geom_hline(aes(yintercept = ul), linetype = 2, col = "#CC6666") +
  geom_vline(xintercept = seq(4.5, 16, 4), alpha = .5) +
  scale_x_continuous(breaks = 1:nrow(Design), trans = reverse_trans()) +
  scale_y_continuous(labels = percent_format(1), breaks = seq(0, 1, .05),
                     sec.axis = sec_axis(trans = ~ ., labels = percent_format(),
                                         breaks = seq(.925, 1, .05))) +
  scale_color_manual(values = cbbPalette[-1]) + guides(col = "none", shape = "none") +
  scale_shape_manual(values = c(1, 2, 4)) +
  labs(y = "ECR of 95% CI", tag = "B") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  coord_flip(xlim = c(15.75, 1.25))
B
A + B + plot_layout(widths = c(2, 1))
ggsave("plots_1b_beta.png", width = 6.5, height = 6)

(Des.dt <- as.data.table(Design))
(Des.dt <- Des.dt[, .N, list(mu1, mu2, sg1, sg2, pop)])
Des.dt

pardefault <- par()
pardefault$mar
pardefault$mai

{
  png("plots_1a_data.png", width = 6.5, height = 4.5, units = "in", res = 180)
  par(mfrow = c(2, 2), mar = c(3.5, 1, 2, 2), yaxt = "n")
  sapply(1:nrow(Des.dt), function (i) {
    xs <- Des.dt[i, ]
    curve(dnorm(qlogis(x), qlogis(xs$mu2), xs$sg2) * 1 / (x * (1 - x)), n = 1e4, col = "red", lwd = 3, lty = 2,
          main = paste0("Data #", i, ", PS = ", percent(xs$pop, .1)), ylab = "density", xlab = "")
    title(xlab = paste0("medians: [", xs$mu1, ", ", xs$mu2, "], variances: [", xs$sg1 ^ 2, ", ", xs$sg2 ^ 2, "]"),
          line = 2)
    curve(dnorm(qlogis(x), qlogis(xs$mu1), xs$sg1) * 1 / (x * (1 - x)), n = 1e4, add = TRUE)
  })
  par(pardefault)
  dev.off()
}

