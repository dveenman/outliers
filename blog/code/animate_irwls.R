# ------------------------------------------------------------------------------
# Code to create an animated gif visualizing the IRWLS process to generate
# a robust M-estimator - generates the same coefficients as MASS::rlm()
# but not the same SEs as they are not adjusted for a robust scale estimate
# ------------------------------------------------------------------------------

suppressMessages({
  library(tidyverse)
  library(animation)
})

huber_weights <- function(resids, k, rhohat) {
  pmin(1, k/abs(resids/rhohat))
}

animate_rrm <- function(
  formula, data, fname, k = 1.345, acc = 1e-3, 
  loss_func = huber_weights, zoom = NULL, return_df = FALSE,
  fixed_coord = FALSE, full_range = FALSE
) {
  data$weights <- 1
  out_df <- NULL
  round <- 1
  delta_resids <- Inf
  while(delta_resids > acc) {
    out_df <- rbind(out_df, cbind(round = round, data))
    mod <- lm(formula, data = data, weights =  weights)
    new_resids <- resid(mod)
    rhohat <- median(abs(new_resids))/0.6745
    new_weights <- loss_func(new_resids, k, rhohat)
    if (round > 1) {
      delta_resids <- sqrt(sum((new_resids - old_resids)^2)/sum(old_resids^2))
    }
    old_resids <- new_resids
    data$weights <- new_weights
    round <- round + 1
  }

  rv <- suppressMessages(saveGIF({
    for (i in 1:max(out_df$round)) {
      step_plot <- function(df, reg) {
        xstr <- all.vars(formula)[2]
        ystr <- all.vars(formula)[1]
        mod <- lm(formula, data = df, weights =  weights)
        if (fixed_coord) {
          abs_max <- max(abs(c(out_df$x, out_df$y)))
          xrange <- extendrange(r = c(-abs_max, abs_max))
          yrange <- xrange
        } else {
          if (!is.null(zoom)) {
            xrange <- extendrange(out_df[df$weights > zoom, xstr])
            yrange <- extendrange(out_df[df$weights > zoom, ystr])
          } else {
            xrange <- extendrange(out_df[, xstr])
            yrange <- extendrange(out_df[, ystr])
          }
        }
        if (full_range) {
          df_x <- tibble(x = seq(xrange[1], xrange[2], length.out = nrow(df)))
          df <- cbind(
            df,
            x_pred = df_x$x,
            as_tibble(
              predict(mod, newdata = df_x, interval = "confidence")
            )
          )
        } else {
          df <- cbind(
            df, x_pred = df[, xstr], predict(mod, interval = "confidence")
          )
        }
        p <- ggplot(
          df, aes_string(x = xstr, y = ystr)
        ) + 
          geom_point(aes(
            color = new_weights < 1,
            alpha = new_weights
          )) + 
          scale_x_continuous(expand = c(0,0), limits = c(-1000, 1000)) +
          scale_y_continuous(expand = c(0,0), limits = c(-1000, 1000)) +
          scale_color_manual(
            name = "Observation down-weighted?",
            values = c("black", "red"), 
            breaks = c(FALSE, TRUE),
            labels = c("No", "Yes"),
            drop = FALSE
          ) +
          scale_alpha(guide = "none") +
          theme_classic()  + 
          theme(legend.position = "bottom")
        
        if(fixed_coord) {
          p <- p + 
            coord_fixed(xlim = c(xrange), ylim = c(yrange)) +
            geom_abline(intercept = 0, slope = 1, lty = 2) 
        } else {
          p <- p + coord_cartesian(xlim = c(xrange), ylim = c(yrange))
        }
        if(!reg) p
        else {
          p + geom_line(
            data = df, aes(x = x_pred, y = fit), color = "blue", 
            size = 1.5, show.legend = FALSE
          ) + geom_ribbon(
            data = df, aes(x = x_pred, ymin = lwr, ymax = upr), alpha = 0.3
          )
        }
      }
      
      df <- out_df %>% filter(round == i)
      df$new_weights <- df$weights
      if (i == 1) {
        print(
          step_plot(df, reg = FALSE) + labs(title = "Iteration 1: Data")
        )
        print(
          step_plot(df, reg = TRUE) + labs(title = "Iteration 1: OLS estimation")
        )
      } 
      if (i < max(out_df$round)) {
        if (i > 1) print(
          step_plot(df, reg = TRUE) + labs(
            title = sprintf("Iteration %d: Reweighted least square estimation", i)
          )
        ) 
        df$new_weights <- out_df$weights[out_df$round == i + 1]
        print(
          step_plot(df, reg = TRUE) + labs(
            title = sprintf("Iteration %d: Adjust weights", i)
          )
        ) 
      } else for(c in 1:2) print(
        step_plot(df, reg = TRUE) + labs(
          title = sprintf("Iteration %d: Final robust estimation", i)
        )
      ) 
    }
  }, movie.name = fname, interval = 2))
  
  # For whatever reason, saveGIF won't store in sub folders...
  if(!file.rename(fname, file.path("output", fname))) 
    stop(sprintf("Could not move '%s' animation file", fname))
  
  if (return_df) out_df
}
