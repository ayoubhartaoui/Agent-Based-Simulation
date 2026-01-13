# ----- Libraries-----

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggnewscale)
library(cowplot)
library(patchwork)
library(tidyverse)
library(ggExtra)
library(truncnorm)


# ----- Paths to the CSV (BE CAREFUL OF THE SIZE)-----
base_path <- "Your file path"

paths <- list(
  low_var = list(
    timeseries = file.path(base_path, "Minimal_law_low_variability_timeseries_summary_ALL.csv"),
    boxplot    = file.path(base_path, "Minimal_law_low_variability_boxplot_lastwindow_ALL.csv"),
    samples    = file.path(base_path, "Minimal_law_low_variability_lastwindow_bias_lambda_samples_ALL.csv")
  ),
  uncorrelated = list(
    timeseries = file.path(base_path, "Minimal_law_uncorrelated_timeseries_summary_ALL.csv"),
    boxplot    = file.path(base_path, "Minimal_law_uncorrelated_boxplot_lastwindow_ALL.csv"),
    samples    = file.path(base_path, "Minimal_law_uncorrelated_lastwindow_bias_lambda_samples_ALL.csv")
  ),
  trade_off = list(
    timeseries = file.path(base_path, "Minimal_law_trade_off_timeseries_summary_ALL.csv"),
    boxplot    = file.path(base_path, "Minimal_law_trade_off_boxplot_lastwindow_ALL.csv"),
    samples    = file.path(base_path, "Minimal_law_trade_off_lastwindow_bias_lambda_samples_ALL.csv")
  ),
  rich_poor = list(
    timeseries = file.path(base_path, "Minimal_law_rich_poor_timeseries_summary_ALL.csv"),
    boxplot    = file.path(base_path, "Minimal_law_rich_poor_boxplot_lastwindow_ALL.csv"),
    samples    = file.path(base_path, "Minimal_law_rich_poor_lastwindow_bias_lambda_samples_ALL.csv")
  )
)


read_one <- function(path, scenario_name) {
  df <- as.data.frame(fread(path))
  df$scenario <- scenario_name
  df
}

read_bp <- function(path, scenario_name) {
  df <- as.data.frame(fread(path))
  df$scenario <- scenario_name
  df
}

read_samples <- function(path, scenario_name) {
  df <- as.data.frame(fread(path))
  df$scenario <- scenario_name
  df
}


# ---- Dataframe conversion----

ts_low <- fread(paths$low_var$timeseries)
bp_to  <- fread(paths$trade_off$boxplot)
sp_rp  <- fread(paths$rich_poor$samples)

df_ts_all <- bind_rows(
  read_one(paths$low_var$timeseries, "low_var"),
  read_one(paths$uncorrelated$timeseries, "uncorrelated"),
  read_one(paths$trade_off$timeseries, "trade_off"),
  read_one(paths$rich_poor$timeseries, "rich_poor")
)

df_bp_all <- bind_rows(
  read_one(paths$low_var$boxplot, "low_var"),
  read_one(paths$uncorrelated$boxplot, "uncorrelated"),
  read_one(paths$trade_off$boxplot, "trade_off"),
  read_one(paths$rich_poor$boxplot, "rich_poor")
)

df_samples_all <- bind_rows(
  read_one(paths$low_var$samples, "low_var"),
  read_one(paths$uncorrelated$samples, "uncorrelated"),
  read_one(paths$trade_off$samples, "trade_off"),
  read_one(paths$rich_poor$samples, "rich_poor")
)

#---- Labels names -----

cost_levels <- c(1.25, 4.0, 7.0)
cost_colors <- c("1.25" = "#d95f02", "4" = "#1b9e77", "7" = "#7570b3")
cost_labels <- c("1.25" = "Low cost (1.25)", "4" = "Intermediate cost (4)", "7" = "High cost (7)")

scenario_labels <- c(
  low_var      = "Low variability",
  uncorrelated = "Uncorrelated variability",
  trade_off    = "Trade-off variability",
  rich_poor    = "Rich-poor variability"
)

#----- Time series of mean belligerence 𝑏 and mean dissatisfaction threshold 𝜆 ----
plot_vs <- function(df_all, scenario_A, scenario_B) {
  
  df_pair <- df_all %>%
    dplyr::filter(scenario %in% c(scenario_A, scenario_B)) %>%
    dplyr::mutate(
      scenario = factor(
        scenario,
        levels = c(scenario_A, scenario_B),
        labels = c(scenario_labels[scenario_A], scenario_labels[scenario_B])
      )
    ) %>%
    dplyr::select(t, Cw, scenario, mean_bias, mean_lambda) %>%
    tidyr::pivot_longer(
      cols = c(mean_lambda, mean_bias),
      names_to = "trait",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      trait = factor(
        trait,
        levels = c("mean_lambda", "mean_bias"),
        labels = c("Mean \u03BB", "Mean b")
      )
    )
  
  panel_labels <- expand.grid(
    trait = levels(df_pair$trait),
    scenario = levels(df_pair$scenario),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      lab = c("(d)", "(b)", "(c)", "(a)"),
      x = 0,
      y = 1
    )
  
  ggplot(df_pair, aes(x = t, y = value, color = factor(Cw))) +
    geom_line(linewidth = 0.8) +
    facet_grid(trait ~ scenario, switch = "y") +
    geom_text(
      data = panel_labels,
      aes(x = x, y = y, label = lab),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.2,
      size = 4
    ) +
    scale_color_manual(
      values = cost_colors,
      labels = cost_labels,
      name = "Conflict costs:"
    ) +
    scale_x_continuous(
      labels = function(x) x / 1000,
      breaks = seq(0, 200000, by = 20000)
    ) +
    coord_cartesian(xlim = c(0, 100000), ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    labs(
      x = "Time step (×1000)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      strip.text.y.right = element_blank(),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
      strip.text = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15)
    )
}

p_low_vs_unc <- plot_vs(df_ts_all, "low_var", "uncorrelated")
p_to_vs_rp   <- plot_vs(df_ts_all, "trade_off", "rich_poor")

print(p_low_vs_unc)
print(p_to_vs_rp)

#---- Boxplot of action proportions among-non-avoiders------
base_colors <- c(war="firebrick", raid="darkorange3", trade="darkgreen")

prep_nonavoid_allcosts <- function(df, scenario_name) {
  df %>%
    filter(scenario == scenario_name, Cw %in% cost_levels) %>%
    mutate(
      Cw = factor(Cw, levels = cost_levels),
      nonavoid = war_mean + raid_mean + trade_mean,
      war_na   = ifelse(nonavoid > 0, war_mean   / nonavoid, NA_real_),
      raid_na  = ifelse(nonavoid > 0, raid_mean  / nonavoid, NA_real_),
      trade_na = ifelse(nonavoid > 0, trade_mean / nonavoid, NA_real_)
    ) %>%
    select(run_id, Cw, war_na, raid_na, trade_na) %>%
    pivot_longer(cols = c(war_na, raid_na, trade_na),
                 names_to = "action", values_to = "value") %>%
    mutate(
      action = recode(action, war_na="war", raid_na="raid", trade_na="trade"),
      action = factor(action, levels = names(base_colors))
    )
}

plot_box_vs_lowvar <- function(df_all, scenario_right) {
  long_low  <- prep_nonavoid_allcosts(df_all, "low_var")
  long_var  <- prep_nonavoid_allcosts(df_all, scenario_right)
  
  action_levels <- levels(long_low$action)
  nudge <- 0.18
  
  long_var <- long_var %>%
    mutate(x_base = as.numeric(factor(action, levels = action_levels)),
           x_pos  = x_base + nudge)
  
  long_low <- long_low %>%
    mutate(x_base = as.numeric(factor(action, levels = action_levels)),
           x_pos  = x_base - nudge)
  
  ggplot() +
    geom_boxplot(
      data = long_var,
      aes(x = x_pos, y = value, fill = action, colour = action),
      width = 0.30,
      outlier.shape = NA,
      linewidth = 0.7
    ) +
    geom_jitter(
      data = long_var,
      aes(x = x_pos, y = value, colour = action),
      width = 0.06,
      alpha = 0.45,
      size = 1.3,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = base_colors, guide = "none") +
    scale_color_manual(
      name = scenario_labels[scenario_right],
      values = base_colors,
      breaks = names(base_colors),
      labels = c(war="War", raid="Raid", trade="Trade")
    ) +
    ggnewscale::new_scale_color() +
    geom_boxplot(
      data = long_low,
      aes(x = x_pos, y = value, fill = action, colour = action),
      width = 0.30,
      outlier.shape = NA,
      linewidth = 0.7,
      linetype = "dashed",
      alpha = 0.35
    ) +
    geom_jitter(
      data = long_low,
      aes(x = x_pos, y = value, colour = action),
      width = 0.06,
      alpha = 0.20,
      size = 1.3,
      show.legend = FALSE
    ) +
    scale_color_manual(
      name = scenario_labels["low_var"],
      values = scales::alpha(base_colors, 0.6),
      breaks = names(base_colors),
      labels = c(war="War", raid="Raid", trade="Trade")
    ) +
    facet_wrap(~ Cw, nrow = 1, labeller = labeller(Cw = cost_labels)) +
    scale_x_continuous(breaks = seq_along(action_levels), labels = action_levels) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
    labs(
      x = "",
      y = "Proportion (%)",
      title = paste0("Action proportions among non-avoiders: Low variability vs ", scenario_labels[scenario_right])
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14),
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 11),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
      strip.text = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)
    )
}

p_unc <- plot_box_vs_lowvar(df_bp_all, "uncorrelated")
p_rp  <- plot_box_vs_lowvar(df_bp_all, "rich_poor")
p_to  <- plot_box_vs_lowvar(df_bp_all, "trade_off")

print(p_unc)
print(p_rp)
print(p_to)


# --- Proportion of engagement ----

df_box_low <- df_bp_all %>% filter(scenario == "low_var")
df_box_unc <- df_bp_all %>% filter(scenario == "uncorrelated")
df_box_to  <- df_bp_all %>% filter(scenario == "trade_off")
df_box_rp  <- df_bp_all %>% filter(scenario == "rich_poor")


df_plot <- df_box_to   # or df_box_low, df_box_unc, df_box_rp

df_sum <- df_plot %>%
  filter(Cw %in% cost_levels) %>%
  mutate(
    nonavoid = war_mean + raid_mean + trade_mean,
    avoid    = inaction_mean
  ) %>%
  group_by(Cw) %>%
  summarise(
    nonavoid = mean(nonavoid, na.rm = TRUE),
    avoid    = mean(avoid,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(nonavoid, avoid),
    names_to = "type",
    values_to = "prop"
  ) %>%
  mutate(
    type = factor(type, levels = c("nonavoid", "avoid"),
                  labels = c("Engage", "Avoid")),
    cost = factor(Cw, levels = cost_levels,
                  labels = cost_labels[as.character(cost_levels)])
  )

ggplot(df_sum, aes(x = cost, y = prop, fill = type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Engage" = "steelblue3",
                               "Avoid"  = "grey70")) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "Proportion (%)",
    fill = NULL,
    title = "d) Frequency of engagement versus avoidance in trade-off variability"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 11)
  )
## ----- the exact purcentage ----

df_box_rp %>%
  filter(Cw %in% cost_levels) %>%
  mutate(nonavoid = war_mean + raid_mean + trade_mean) %>%
  group_by(
    cost = factor(Cw, levels = cost_levels,
                  labels = cost_labels[as.character(cost_levels)])
  ) %>%
  summarise(
    mean_engage = mean(nonavoid, na.rm = TRUE),
    mean_avoid  = mean(inaction_mean, na.rm = TRUE),
    .groups = "drop"
  )


#----- Within population anylsis -----

scenario_name <- "trade_off" #choose the environement accordingly
cw_value <- 1.25 # you can choose the Conflict costs accordingly

ts   <- as.data.frame(fread(paths_ts[[scenario_name]]))
snap <- as.data.frame(fread(paths_snap[[scenario_name]]))

ts_cw <- ts %>% filter(Cw == cw_value)
snap_cw <- snap %>% filter(Cw == cw_value)

t_min_eq <- max(snap_cw$gen_in_window, na.rm = TRUE) - 2000 + 1
t_max_eq <- max(snap_cw$gen_in_window, na.rm = TRUE)

snap_eq <- snap_cw %>% filter(gen_in_window >= t_min_eq, gen_in_window <= t_max_eq)

trait_summary <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  list(mean = m, sd = s, sd_low = max(0, m - s), sd_high = min(1, m + s))
}

s_b <- trait_summary(snap_eq$bias)
s_l <- trait_summary(snap_eq$lambda)

## ---- Historgram of belligerence 𝑏 and dissatisfaction threshold 𝜆----

p_hist_bias <- ggplot(snap_eq, aes(x = bias)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.025, boundary = 0, closed = "left",
                 fill = "blue", color = "black", alpha = 0.35) +
  geom_vline(xintercept = s_b$mean, linewidth = 1.0) +
  geom_vline(xintercept = c(s_b$sd_low, s_b$sd_high),
             linetype = "dotted", linewidth = 0.9) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Belligerence (b)", y = "Frequency", title = "b) Histogram") +
  theme_bw()

p_hist_lambda <- ggplot(snap_eq, aes(x = lambda)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.025, boundary = 0, closed = "left",
                 fill = "blue", color = "black", alpha = 0.35) +
  geom_vline(xintercept = s_l$mean, linewidth = 1.0) +
  geom_vline(xintercept = c(s_l$sd_low, s_l$sd_high),
             linetype = "dotted", linewidth = 0.9) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Dissatisfaction threshold (\u03BB)", y = "", title = "c) Histogram ") +
  theme_bw()

##---- mean trait with its variance ----

set.seed(1)

trait_to_plot <- "bias"
mean_col <- if (trait_to_plot == "bias") "mean_bias" else "mean_lambda"
sd_col   <- if (trait_to_plot == "bias") "sd_bias_within" else "sd_lambda_within"

line_df <- ts_cw %>%
  transmute(
    t = as.numeric(t),
    mean = .data[[mean_col]],
    var  = (.data[[sd_col]])^2
  )

p_trait_time <- ggplot(line_df, aes(x = t/1000)) +
  geom_line(aes(y = mean), linewidth = 1.1, color = "black") +
  geom_line(aes(y = var),  linewidth = 1.0, color = "red") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Time (generations \u00D7 1000)",
    y = if (trait_to_plot == "bias") "Belligerence (b)" else "Dissatisfaction threshold (\u03BB)",
    title = "a) Dynamics of belligerence distribution"
  ) +
  theme_bw()

##---- Heatmap ----

bins_2d <- 100 # choose this number accordingly to the number of agent you sampled in the julia

p_heatmap <- ggplot(snap_eq, aes(x = bias, y = lambda)) +
  stat_bin2d(bins = bins_2d) +
  scale_fill_continuous(trans = "reverse") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Belligerence (b)",
    y = "Dissatisfaction threshold (\u03BB)",
    title = "d) Heatmap distribution of traits b and \u03BB",
    fill = "Density"
  ) +
  theme_bw()

# (p_trait_time / (p_hist_bias | p_hist_lambda)) | p_heatmap #THIS IS THE BASIC ARCHITECTURE
scenario_name <- "low_var";      print((p_trait_time / (p_hist_bias | p_hist_lambda)) | p_heatmap)
scenario_name <- "uncorrelated"; print((p_trait_time / (p_hist_bias | p_hist_lambda)) | p_heatmap)
scenario_name <- "trade_off";    print((p_trait_time / (p_hist_bias | p_hist_lambda)) | p_heatmap)
scenario_name <- "rich_poor";    print((p_trait_time / (p_hist_bias | p_hist_lambda)) | p_heatmap)





#----- Distribution of resources (r1,r2) -----

set.seed(1)

n_heat <- 120000
n_pts  <- 3000
bins2d <- 60

# env_key is stable; title is just display text
env_titles <- c(
  lowvar   = "a) Low variability",
  uncorr   = "b) Uncorrelated",
  richpoor = "c) Rich--poor",
  tradeoff = "d) Trade-off"
)

draw_env <- function(env_key, n) {
  
  if (env_key == "lowvar") {
    r1 <- rtruncnorm(n, 0, 1, mean = 0.5, sd = 0.1)
    r2 <- rtruncnorm(n, 0, 1, mean = 0.5, sd = 0.1)
  } else if (env_key == "uncorr") {
    r1 <- runif(n)
    r2 <- runif(n)
  } else if (env_key == "richpoor") {
    rich <- rbinom(n, 1, 0.5) == 1
    r1 <- ifelse(rich, rbeta(n, 15.0, 1.5), rbeta(n, 1.5, 15.0))
    r2 <- ifelse(rich, rbeta(n, 15.0, 1.5), rbeta(n, 1.5, 15.0))
  } else if (env_key == "tradeoff") {
    swap <- rbinom(n, 1, 0.5) == 1
    high <- rbeta(n, 5, 1.5)
    low  <- rbeta(n, 1.5, 5)
    r1 <- ifelse(swap, high, low)
    r2 <- ifelse(swap, low, high)
  } else {
    stop("Unknown env_key: ", env_key)
  }
  
  tibble(r1 = r1, r2 = r2)
}

make_panel <- function(env_key) {
  df <- draw_env(env_key, n_heat)
  
  idx <- sample.int(nrow(df), size = min(n_pts, nrow(df)), replace = FALSE)
  pts <- df[idx, , drop = FALSE]
  
  p <- ggplot(df, aes(r1, r2)) +
    geom_bin2d(bins = bins2d) +
    scale_fill_continuous(trans = "reverse") +
    geom_point(data = pts, alpha = 0) +
    coord_fixed() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = env_titles[[env_key]], x = NULL, y = NULL, fill = "Count") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "none"
    )
  
  ggMarginal(
    p,
    type = "histogram",
    size = 5,
    fill = "steelblue",
    colour = "steelblue"
  )
}

p1 <- make_panel("lowvar")
p2 <- make_panel("uncorr")
p3 <- make_panel("richpoor")
p4 <- make_panel("tradeoff")

grid_2x2 <- plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv")

# Legend (use any env_key that exists!)
df_leg <- draw_env("lowvar", 40000)
p_leg <- ggplot(df_leg, aes(r1, r2)) +
  geom_bin2d(bins = bins2d) +
  scale_fill_continuous(trans = "reverse") +
  coord_fixed() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(fill = "Count") +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )

legend <- get_legend(p_leg + theme(legend.position = "right"))

with_legend <- plot_grid(grid_2x2, legend, ncol = 2, rel_widths = c(1, 0.12))

final <- ggdraw(with_legend) +
  draw_label("Resource 1", x = 0.46, y = 0.02, size = 12) +
  draw_label("Resource 2", x = 0.02, y = 0.52, angle = 90, size = 12)

final

