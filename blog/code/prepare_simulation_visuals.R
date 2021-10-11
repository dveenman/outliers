suppressMessages({
  library(tidyverse)
  library(gganimate)
  library(animation)
})

source("code/animate_irwls.R")

message("Preparing simulation visuals... ", appendLF = FALSE)

img_width <- 480
set.seed(27)
nobs <- 100

base_sample <- tibble(
  x = rnorm(nobs),
  e = rnorm(nobs),
  y = x + e
)

anim_smp <- function(vertical, bs = base_sample) {
  create_sample <- function(obs, oset){
    df <- bs
    if (vertical) df$y[obs] <- df$y[obs] + oset
    else df$x[obs] <- df$x[obs] + oset
    df
  }

  dist_center <- sqrt(bs$x^2 + bs$y^2)
  center_obs <- which(dist_center == min(dist_center))
  df <- tibble(
    frame = 1:100,
    obs = center_obs,
    oset = {if(vertical) seq(0, 9.8 - bs$y[center_obs], length.out = 100)
      else seq(0, 9.8 - bs$x[center_obs], length.out = 100)}
  )

  do.call(
    rbind,
    lapply(
      1:100,
      function(x) cbind(frame = x, create_sample(df$obs[x], df$oset[x]))
    )
  )
}

sim_scatter_plot <- function(smp, coord = TRUE, reg = TRUE, diag = TRUE) {
  p <- ggplot(smp, aes(x = x, y = y)) +
    geom_point() +
    scale_x_continuous(expand=c(0,0), limits=c(-100,100)) +
    scale_y_continuous(expand=c(0,0), limits=c(-100,100)) +
    theme_classic()

  if (diag) p <- p + geom_abline(intercept = 0, slope = 1, lty = 2)
  if (reg) p <- p +
      geom_smooth(
        method = "lm", formula = "y ~ x", fullrange = TRUE, na.rm = TRUE
      )
  if (coord) p <- p + coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10))

  p
}

create_animation <- function(smp) {
  sim_scatter_plot(smp) + transition_manual(frame)
}

fig_sim_base <- sim_scatter_plot(base_sample)
ani_sim_bad_leverage <- create_animation(anim_smp(vertical = FALSE))
ani_sim_vert_outlier <- create_animation(anim_smp(vertical = TRUE))

ggsave(
  "output/fig_sim_base.png", fig_sim_base,
  width = img_width/72, height = img_width/72, dpi = 72
)
anim_save(
  "output/ani_sim_bad_leverage.gif",
  ani_sim_bad_leverage, end_pause = 20, width = img_width, height = img_width
)
anim_save(
  "output/ani_sim_vert_outlier.gif",
  ani_sim_vert_outlier, end_pause = 20, width = img_width, height = img_width
)

# --- RR Simulation ------------------------------------------------------------


# Brute force randomize until we have a sample where robust regression
# has a huge effect on the coefficient
largest_delta <- 0
rcoef <- 0
while (largest_delta < 0.4 | rcoef < 0.90) {
  nobs <- 100
    os <- tibble(
    x = rnorm(nobs),
    e = rnorm(nobs),
    y = x + e
  )
  outlier_x <- sample(nobs, 10)
  outlier_y <- sample(nobs, 10)
  os$x[outlier_x] <- os$x[outlier_x] + 3*rnorm(10)
  os$y[outlier_y] <- os$y[outlier_y] + 2*os$e[outlier_y]

  rcoef <- suppressWarnings(MASS::rlm(y ~ x, data = os)$coefficients[2])
  lscoef <- lm(y ~ x, data = os)$coefficients[2]
  delta_coef <- rcoef - lscoef
  if (delta_coef > largest_delta) {
    largest_delta <- delta_coef
    best_sample <- os
  }
}

abs_max <- max(abs(c(best_sample$x, best_sample$y)))
xrange <- extendrange(r = c(-abs_max, abs_max))
yrange <- xrange

fig_sim_irwls_base <- sim_scatter_plot(best_sample, coord = FALSE) +
  coord_fixed(xlim = c(xrange), ylim = c(yrange))

ggsave(
  "output/fig_sim_irwls_base.png", fig_sim_irwls_base,
  width = img_width/72, height = img_width/72, dpi = 72
)

df <- animate_rrm(
  y ~ x, best_sample, "ani_sim_irwls.gif", k = 1.345, return_df = TRUE,
  full_range = TRUE, fixed_coord = TRUE
)

# Use the below to verify that our animation code yields the same coefficients
# (but not the same SEs) as MASS:rlm()


if (FALSE) {
  fs <- df %>% filter(round == max(round))
  summary(lm(y ~ x, data = fs))
  summary(lm(y ~ x, weights = weights, data = fs))
  mod_rlm <- MASS::rlm(y ~ x, data = fs, acc = 1e-3)
  summary(mod_rlm)
}

est_nl_models <- function(n = 100) {
  nl <- tibble(
    x = rnorm(n),
    y = x^3 + x^2 + x + rnorm(n)
  )
  ols <- lm(y ~x, data = nl)
  rr <- MASS::rlm(y ~x, data = nl)
  list(
    nl,
    tibble(
      ols_const = unname(ols$coefficients[1]),
      ols_x = unname(ols$coefficients[2]),
      rr_const = unname(rr$coefficients[1]),
      rr_x = unname(rr$coefficients[2])
    )
  )
}

mod_list <- lapply(1:100, function(x) est_nl_models())
dta <- do.call(rbind, lapply(mod_list, function(x) x[[1]]))
coefs <- do.call(rbind, lapply(mod_list, function(x) x[[2]]))
fig_sim_nl_100runs <- ggplot(coefs) +
  geom_point(data = dta, aes(x = x, y = y), alpha = 0.1) +
  geom_abline(
    aes(intercept = ols_const, slope = ols_x), alpha = 0.1, color = "blue"
  ) +
  geom_abline(
    aes(intercept = rr_const, slope = rr_x), alpha = 0.1, color = "red"
  ) +
  theme_classic()

ggsave(
  "output/fig_sim_nl_100runs.png", fig_sim_nl_100runs,
  width = img_width/72, height = img_width/72, dpi = 72
)


delta_coefs <- 0
while(delta_coefs < 2) {
 lst <- est_nl_models()
 nl <- lst[[1]]
 delta_coefs <- lst[[2]]$ols_x - lst[[2]]$rr_x
}

xrange <- extendrange(nl$x)
yrange <- extendrange(nl$y)

fig_sim_nl_base <- sim_scatter_plot(nl, coord = FALSE, diag = FALSE) +
  coord_cartesian(xlim = c(xrange), ylim = c(yrange))

ggsave(
  "output/fig_sim_nl_base.png", fig_sim_nl_base,
  width = img_width/72, height = img_width/72, dpi = 72
)

if (FALSE) {
  fig_sim_nl_smooth <- sim_scatter_plot(nl, reg = FALSE, coord = FALSE, diag = FALSE) +
    geom_smooth(method = 'loess', formula = 'y ~ x') +
    coord_cartesian(xlim = c(xrange), ylim = c(yrange))

  ggsave(
    "output/fig_sim_nl_smooth.png", fig_sim_nl_smooth,
    width = img_width/72, height = img_width/72, dpi = 72
  )
}

animate_rrm(
  y ~ x, nl, "ani_sim_nl_irwls.gif", fixed_coord = FALSE, full_range = TRUE
)

save(
  list = ls(pattern = "(ani|fig|tab)_"),
  file = "output/simulation_visuals.rda"
)

message("done")
