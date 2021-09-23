source("ps_functions_new.R")

library(SimDesign)

(Design <- expand.grid(
  pop = .5, n1 = 10, p = c(1 / 3, .5), icc = c(.05, .2), n2 = 3e1))
# (condition <- list(pop = .5, n1 = 10, p = 1 / 3, icc = .2, n2 = 3e1))

#-------------------------------------------------------------------

Generate <- function(condition, fixed_objects = NULL) {
  dat <- data.frame(ID = 1:condition$n2, x = rbinom(condition$n2, 1, condition$p),
                    g = rnorm(condition$n2, 0, sqrt(condition$icc)))
  dat <- dat[rep(1:condition$n2, round(rnorm(condition$n2, condition$n1))), ]
  dat$y <- rnorm(nrow(dat), dat$g, sqrt(1 - condition$icc))
  dat
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
  vargha <- ps(dat$y[dat$x == 1], dat$y[dat$x == 0])

  vargha.ci <- inf.ps.boot(dat)
  vargha.ecr <- vargha.ci[1] < condition$pop & vargha.ci[2] > condition$pop
  vargha.mse <- (vargha - condition$pop) ^ 2

  logit.ci <- inf.ps.logit(length(dat$y), vargha)
  logit.ecr <- logit.ci[1] < condition$pop & logit.ci[2] > condition$pop

  logit.m.res <- inf.ps.logit.crve(dat)
  logit.m <- logit.m.res[1]
  logit.m.ci <- logit.m.res[2:3]
  logit.m.ecr <- logit.m.ci[1] < condition$pop & logit.m.ci[2] > condition$pop
  logit.m.mse <- (logit.m - condition$pop) ^ 2

  pl.m.res <- inf.ps.pl.mlm(dat)
  pl.m <- pl.m.res[1]
  pl.m.ci <- pl.m.res[2:3]
  pl.m.ecr <- pl.m.ci[1] < condition$pop & pl.m.ci[2] > condition$pop
  pl.m.mse <- (pl.m - condition$pop) ^ 2

  c(vargha = vargha, logit.m = logit.m, pl.m = pl.m,
    vargha.mse = vargha.mse, logit.m.mse = logit.m.mse, pl.m.mse = pl.m.mse,
    vargha.ecr = vargha.ecr, logit.ecr = logit.ecr, logit.m.ecr = logit.m.ecr,
    pl.m.ecr = pl.m.ecr)
}

Summarise <- function(condition, results, fixed_objects = NULL) {
  colMeans(results)
}

#-------------------------------------------------------------------

unlink("SimDesign-results_pop-os", recursive = TRUE)
results <- runSimulation(
  design = Design, replications = 2e3, generate = Generate,
  analyse = Analyse, summarise = Summarise, seed = rep(123, 4),
  save = TRUE, save_results = TRUE, progress = TRUE,
  filename = "./simdata_0", packages = c("boot", "sandwich", "lmtest"))

#-------------------------------------------------------------------

library(ggplot2)
library(ggforce)
library(scales)
library(patchwork)
library(data.table)
theme_set(theme_bw())

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

(Design <- expand.grid(
  pop = .5, n1 = 10, p = c(1 / 3, .5), icc = c(.05, .2), n2 = 3e1))
Design$design <- 1:4
Design

raw <- rbindlist(lapply(list.files("SimDesign-results_pop-os/"), function (x) {
  as.data.frame(readRDS(paste0("SimDesign-results_pop-os/", x))$results)
}), idcol = "design")

raw.l <- reshape(
  raw, idvar = "ID", varying = list(c(2, 2, 3, 4), c(5, 5, 6, 7), 8:11),
  v.names = c("estimate", "mse", "coverage"),
  direction = "long")
raw.l[, time.t := c("Vargha\nw. BCa", "Fractional\nw. Wald", "Fractional\nw. CRVE",
                    "Placement\nw. MLM")[time]]
raw.l

raw.l$pop <- .5

sum.dat.0 <- raw.l[
  , .(mean = mean(estimate), q05 = quantile(estimate, .05),
      q25 = quantile(estimate, .25), q75 = quantile(estimate, .75),
      q95 = quantile(estimate, .95), mse = mean((estimate - pop) ^ 2), pop = mean(pop)),
  list(design, time.t, time)]
sum.dat.0 <- merge(sum.dat.0, Design)
sum.dat.0

dodge <- position_dodge(width = .75)
A <- ggplot(sum.dat.0[time != 2], aes(design, mean, col = time.t)) +
  geom_point(aes(shape = time.t), position = dodge, size = 2.5) +
  geom_linerange(aes(ymin = q05, ymax = q95), position = dodge) +
  geom_hline(aes(yintercept = pop), linetype = 2) +
  geom_text(aes(x = design - .125, label = number(mse * 1e3, .01), y = .56), position = dodge, size = 2.5) +
  scale_x_continuous(breaks = 1:4, trans = reverse_trans(),
                     labels = paste0("ICC = ", percent(Design$icc), "\n% Treat = ", percent(Design$p))) +
  scale_y_continuous(labels = percent_format(1)) +
  scale_color_manual(values = cbbPalette[c(2, 3, 4)]) +
  scale_shape_manual(values = c(1, 2, 4)) +
  guides(col = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  labs(x = "Condition", y = "Distribution of probability of superiority",
       tag = "A", col = "", shape = "") +
  theme(legend.position = "bottom", axis.ticks.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0),
        legend.title = element_blank(), legend.spacing.y = unit(0, "mm"),
        panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank()) +
  coord_flip()
A

sum.dat.1 <- raw.l[, .(mean = mean(coverage), count = .N), list(design, time.t, time)]
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

B <- ggplot(sum.dat.1[time != 2], aes(design, mean, col = time.t)) +
  geom_point(aes(shape = time.t), size = 1, position = dodge) +
  geom_linerange(aes(ymin = llmjf, ymax = ulmjf), position = dodge) +
  geom_hline(aes(yintercept = exp), linetype = 2, alpha = .5) +
  geom_hline(aes(yintercept = ll), linetype = 2, col = "#CC6666") +
  geom_hline(aes(yintercept = ul), linetype = 2, col = "#CC6666") +
  scale_x_continuous(trans = reverse_trans()) +
  scale_y_continuous(labels = percent_format(1), breaks = seq(0, 1, .05),
                     sec.axis = sec_axis(trans = ~ ., labels = percent_format(),
                                         breaks = seq(.925, 1, .05))) +
  scale_color_manual(values = cbbPalette[c(2, 3, 4)]) +
  scale_shape_manual(values = c(1, 2, 4)) +
  guides(col = "none", shape = "none") +
  labs(y = "ECR of 95% CI", tag = "B") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0),
        legend.title = element_blank(), legend.spacing.y = unit(0, "mm"),
        axis.text.x = element_text(angle = 30, vjust = 1),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank()) +
  coord_flip()
B
A + B + plot_layout(widths = c(1.3, 1))
ggsave("plots_0_50.png", width = 6.5, height = 3.5)

