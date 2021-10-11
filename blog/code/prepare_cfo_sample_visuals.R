suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(huxtable)
  library(ExPanDaR)
})

source("code/animate_irwls.R")
img_width <- 480

message("Preparing cfo/tacc visuals... ", appendLF = FALSE)

simple_scatter_plot <- function(df, alpha = 1) {
  ggplot(df, aes(x = cfo, y = tacc)) +
    geom_point(alpha = alpha) +
    scale_x_continuous(expand=c(0,0), limits=c(-100,100)) +
    scale_y_continuous(expand=c(0,0), limits=c(-100,100)) +
    theme_classic()  +
    coord_cartesian(xlim = extendrange(df$cfo), ylim = extendrange(df$tacc))
}


os <- readRDS("data/generated/sample.rds")

fig_cfo_scatter <- simple_scatter_plot(os)
ggsave(
  "output/fig_cfo_scatter.png", fig_cfo_scatter,
  width = img_width/72, height = img_width/72, dpi = 72
)


animate_rrm(tacc ~ cfo, os, "ani_cfo_linear_irwls.gif", zoom = 0.1)

os$rr_weight <- MASS::rlm(tacc ~cfo, data = os)$w


tab_cfo_rrweights_by_ind <- os %>%
  group_by(`Fama/French Industry` = ff12_ind) %>%
  summarise(
    `% down-weighted` = sum(rr_weight != 1)/n()
  ) %>%
  arrange(`% down-weighted`) %>%
  as_hux() %>%
  theme_basic() %>%
  set_number_format(-1, 2, fmt_percent(1))

quick_html(
  tab_cfo_rrweights_by_ind, file = "output/tab_cfo_rrweights_by_ind.html",
  open = FALSE
)

lollipop_rrweights <- function(var, varstr) {
  os %>%
    mutate(decile = factor(ntile({{var}}, 10))) %>%
    group_by(decile) %>%
    summarise(pct_dw = sum(rr_weight != 1)/n()) %>%
    ggplot(aes(x = decile, y = pct_dw)) +
    geom_point(color = "blue", size = 4) +
    geom_segment(
      aes(x = decile, xend = decile, y=0, yend = pct_dw)
    ) +
    labs(
      title = sprintf("Share of down-weighted observations by %s decile", varstr),
      x = NULL, y = NULL
    ) +
    scale_y_continuous(
      expand = expansion(add = c(0,0.05)),
      labels = scales::percent_format(accuracy = 1)
    ) +
    theme_classic()
}

fig_cfo_rrweights_by_size <- lollipop_rrweights(ln_ta, "size")
fig_cfo_rrweights_by_sales_growth <- lollipop_rrweights(sales_gr, "sales growth")

ggsave(
  "output/fig_cfo_rrweights_by_size.png", fig_cfo_rrweights_by_size,
  width = img_width/72, height = img_width/72, dpi = 72
)

ggsave(
  "output/fig_cfo_rrweights_by_sales_growth.png",
  fig_cfo_rrweights_by_sales_growth,
  width = img_width/72, height = img_width/72, dpi = 72
)


if (FALSE) {
  animate_rrm(tacc ~ (cfo<0)*cfo, os, "ani_cfo_iacted_irwls.gif", zoom = 0.1)

  os <- os  %>% treat_outliers()

  animate_rrm(tacc ~ cfo, os, "ani_cfo_linear_irwls_win.gif")
  df <- animate_rrm(
    tacc ~ (cfo<0)*cfo, os, "ani_cfo_iacted_irwls_win.gif", return_df = TRUE
  )

}

# Use the below to verify that our animation code yields the same coefficients
# (but not the same SEs) as MASS:rlm()

if (FALSE) {
  fs <- df %>% filter(round == max(round))
  summary(lm(tacc ~ (cfo<0)*cfo, data = fs))
  summary(lm(tacc ~ (cfo<0)*cfo, weights = weights, data = fs))
  mod_rlm <- MASS::rlm(tacc ~ (cfo<0)*cfo, data = fs, acc = 1e-3)
  summary(mod_rlm)
}

smp <- readRDS("data/generated/sample.rds")

smp_win <- smp %>%
  treat_outliers()

pred <- function(mod) {
  as_tibble(cbind(
    cfo = smp_win$cfo,
    predict(mod, newdata = smp_win[, "cfo"], interval = "confidence")
  ))
}
rrmod <- MASS::rlm(tacc ~ cfo, data = smp)
rr_pred <- pred(rrmod)

fig_cfo_linear <- simple_scatter_plot(smp_win, alpha = 0.1) +
  geom_smooth(method = 'lm', formula = 'y ~ x', color = "blue") +
  geom_line(
    data = rr_pred, aes(x = cfo, y = fit),
    color = "red", size = 1.5, show.legend = FALSE
  ) +
  geom_ribbon(
    data = rr_pred, aes(x = cfo, y = fit, ymin = lwr, ymax = upr),
    alpha = 0.3
  )

ggsave(
  "output/fig_cfo_linear.png", fig_cfo_linear,
  width = img_width/72, height = img_width/72, dpi = 72
)

olsmod <- lm(tacc ~ (cfo < 0)*cfo, data = smp_win)
ols_pred <- pred(olsmod)

rrmod <- MASS::rlm(tacc ~ (cfo < 0)*cfo, data = smp)
rr_pred <- pred(rrmod)

fig_cfo_iacted <- simple_scatter_plot(smp_win, alpha = 0.1) +
  geom_line(
    data = ols_pred, aes(x = cfo, y = fit),
    color = "blue", size = 1.5, show.legend = FALSE
  ) +
  geom_ribbon(
    data = ols_pred, aes(x = cfo, y = fit, ymin = lwr, ymax = upr),
    alpha = 0.3
  ) +
  geom_line(
    data = rr_pred, aes(x = cfo, y = fit),
    color = "red", size = 1.5, show.legend = FALSE
  ) +
  geom_ribbon(
    data = rr_pred, aes(x = cfo, y = fit, ymin = lwr, ymax = upr),
    alpha = 0.3
  )

ggsave(
  "output/fig_cfo_iacted.png", fig_cfo_iacted,
  width = img_width/72, height = img_width/72, dpi = 72
)

save(
  list = ls(pattern = "(ani|fig|tab)_"),
  file = "output/cfo_sample_visuals.rda"
)

message("done")
