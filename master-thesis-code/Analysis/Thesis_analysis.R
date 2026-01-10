base_pkgs <- c(
  "package:stats", "package:graphics", "package:grDevices",
  "package:utils", "package:datasets", "package:methods",
  "package:base"
)

for (p in search()) {
  if (startsWith(p, "package:") && !(p %in% base_pkgs)) {
    detach(p, unload = TRUE, character.only = TRUE)
  }
}


# ----- 1. Libraries --------
#window()
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggnewscale)
library(cowplot)


# ---- Paths ----
## Principal model paths 
paths <- list(
  low_var      = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_no_var/Minimal_Law_low_variation.csv",
  uncorrelated = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_uniform/Minimal_Law_uniform.csv",
  trade_off    = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_law_trade_off/Minimial_Law_trade_off.csv",
  rich_poor    = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_rich_poor/Minimal_law_rich_poor.csv"
)

## Non-avoider model paths 
paths <- list(
  low_var      = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_no_var/Minimal_no_var_normal.csv",
  uncorrelated = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_uniform/ML_uniform.csv",
  trade_off    = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_law_trade_off/ML_trade_off.csv",
  rich_poor    = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_rich_poor/Minimal_law_rich_poor_nonavoid.csv"
)

## cobb-douglas
paths <- list(
  low_var      = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/CB_no_varr/cb_novar_normal_a_test.csv",
  uncorrelated = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/CB_uniform/cb_uniform_a.csv",
  trade_off    = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/CB_trade_off/CobbDouglas_trade_off_bias_lambda.csv",
  rich_poor    = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/CB_rich_poor/CobbDouglas_rich_poor_bias_lambda.csv"
)
  
read_one <- function(path, scenario_name) {
  df <- as.data.frame(fread(path))
  df$scenario <- scenario_name
  df
}

df_all <- bind_rows(
  read_one(paths$low_var, "low_var"),
  read_one(paths$uncorrelated, "uncorrelated"),
  read_one(paths$trade_off, "trade_off"),
  read_one(paths$rich_poor, "rich_poor")
)

# ---- Shared legend settings ----
cost_colors <- c("1.25" = "#d95f02", "4" = "#1b9e77", "7" = "#7570b3")
cost_labels <- c("1.25" = "Low cost (1.25)", "4" = "Intermediate cost (4)", "7" = "High cost (7)")

scenario_labels <- c(
  low_var = "Low variability",
  uncorrelated = "Uncorrelated variability",
  trade_off = "Trade-off variability",
  rich_poor = "Rich-poor variability"
)

# ---- Function: Belligerence and Dissatisfaction treshold Figures in 2x2 grid----
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
      # order is: (row1,col1) (row1,col2) (row2,col1) (row2,col2)
      lab = c("(d)", "(b)", "(c)", "(a)"),
      x = 0,
      y = 1
    )
  
  ggplot(df_pair, aes(x = t, y = value, color = factor(Cw / 2))) +
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



p_low_vs_unc <- plot_vs(df_all, "low_var", "uncorrelated")
p_to_vs_rp   <- plot_vs(df_all, "trade_off", "rich_poor" )

print(p_low_vs_unc)
print(p_to_vs_rp)





# save
#ggsave("traits_low_vs_uncorrelated.png", p_low_unc, width = 12,height = 7,dpi = 600)

#ggsave("traits_low_vs_tradeoff.png", p_low_to, width = 12,height = 7,dpi = 600)

#ggsave("traits_low_vs_richpoor.png",p_low_rp,width = 12,height = 7,dpi = 600)



# ----3. Boxplot Figure of actions  ----

paths_box <- list(
  low_var   = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_no_var/boxplot_ML_no_varriation.csv",
  uncorr    = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_uniform/boxplot_ML_uniform.csv",
  tradeoff  = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_law_trade_off/boxplot_ML_trade_off.csv",
  richpoor  = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_rich_poor/boxplot_rich_poor.csv"
)

## cobb-douglas paths
paths_box <- list(
  low_var   = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/CB_no_varr/boxplot_low_variation_cobb.csv",
  uncorr    = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/CB_uniform/boxplot_uncorrelated_cobb.csv",
  tradeoff  = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/CB_trade_off/boxplot_trade_off_cobb.csv",
  richpoor  = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/CB_rich_poor/boxplot_rich_poor_cobb.csv"
)


read_box <- function(path) as.data.frame(fread(path))

df_low <- read_box(paths_box$low_var)
df_unc <- read_box(paths_box$uncorr)
df_to  <- read_box(paths_box$tradeoff)
df_rp  <- read_box(paths_box$richpoor)

base_colors <- c(war="firebrick", raid="darkorange3", trade="darkgreen")

cost_levels <- c(2.5, 8, 14)
cost_labels <- c(
  "2.5" = "(a) Low cost",
  "8"   = "(b) Intermediate cost",
  "14"  = "(c) High cost"
)

prep_nonavoid_allcosts <- function(df, model_label) {
  df %>%
    filter(Cw %in% cost_levels) %>%
    mutate(
      nonavoid = war + raid + trade,
      war_na   = ifelse(nonavoid > 0, war   / nonavoid, NA_real_),
      raid_na  = ifelse(nonavoid > 0, raid  / nonavoid, NA_real_),
      trade_na = ifelse(nonavoid > 0, trade / nonavoid, NA_real_)
    ) %>%
    select(run_id, Cw, war_na, raid_na, trade_na) %>%
    pivot_longer(cols = c(war_na, raid_na, trade_na),
                 names_to = "action", values_to = "value") %>%
    mutate(
      model = model_label,
      action = recode(action, war_na="war", raid_na="raid", trade_na="trade"),
      action = factor(action, levels = names(base_colors)),
      cost = factor(Cw, levels = cost_levels, labels = cost_labels[as.character(cost_levels)])
    )
}

plot_box_allcosts <- function(df_low, df_var, var_name) {
  
  long_var <- prep_nonavoid_allcosts(df_var, var_name)
  long_low <- prep_nonavoid_allcosts(df_low, "Low variability")
  
  action_levels <- levels(long_var$action)
  nudge <- 0.18
  
  long_var <- long_var %>%
    mutate(
      x_base = as.numeric(factor(action, levels = action_levels)),
      x_pos  = x_base + nudge   # ← var on the RIGHT
    )
  
  long_low <- long_low %>%
    mutate(
      x_base = as.numeric(factor(action, levels = action_levels)),
      x_pos  = x_base - nudge   # ← low var on the LEFT
    )

  
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
      size = 1.5,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = base_colors, guide = "none") +
    scale_color_manual(
      name = var_name,
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
      size = 1.5,
      show.legend = FALSE
    ) +
    scale_color_manual(
      name = "Low variability",
      values = scales::alpha(base_colors, 0.6),
      breaks = names(base_colors),
      labels = c(war="War", raid="Raid", trade="Trade")
    ) +
    facet_wrap(~ cost, nrow = 1) +
    scale_x_continuous(
      breaks = seq_along(action_levels),
      labels = action_levels
    ) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
    labs(
      x = "",
      y = "Proportion (%)",
      title = paste0("Action proportions among non-avoiders: Low variability vs ", var_name)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),
      legend.position = "right",
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 12),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
      strip.text = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15)
    )
}



p_unc_all <- plot_box_allcosts(df_low, df_unc, "Uncorrelated variability")
p_to_all  <- plot_box_allcosts(df_low, df_to,  "Trade-off variability")
p_rp_all  <- plot_box_allcosts(df_low, df_rp,  "Rich-poor variability")

print(p_unc_all)
print(p_to_all)
print(p_rp_all)






# ---- engagment proportions-----
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
cost_labels <- c(
  "2.5" = "Low cost",
  "8"   = "Intermediate cost",
  "14"  = "High cost"
)
df_plot <- df_to   # or df_unc, df_to, df_rp

df_sum <- df_plot %>%
  filter(Cw %in% cost_levels) %>%
  mutate(
    nonavoid = war + raid + trade,
    avoid    = 1 - nonavoid
  ) %>%
  group_by(Cw) %>%
  summarise(
    nonavoid = mean(nonavoid, na.rm = TRUE),
    avoid    = mean(avoid,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(nonavoid, avoid),
               names_to = "type", values_to = "prop") %>%
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
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "Proportion (%)",
    fill = NULL,
    title ="d) Frequency of engagement versus avoidance in trade-off variability"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

## exact numbers
df_to %>%
  filter(Cw %in% cost_levels) %>%
  mutate(nonavoid = war + raid + trade) %>%
  group_by(cost = factor(Cw, levels = cost_levels,
                         labels = cost_labels[as.character(cost_levels)])) %>%
  summarise(mean_engage = mean(nonavoid, na.rm = TRUE),
            .groups = "drop")


## model for avoidance as an action.

## Non-avoider model paths 
paths_box <- list(
  low_var   = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_no_var/boxplot_low_variation_nonavoid.csv",
  uncorr    = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_uniform/boxplot_uniform_unavoid.csv",
  tradeoff  = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_law_trade_off/boxplot_trade_nonavoid.csv",
  richpoor  = "C:/Users/sabof/OneDrive - Université de Lausanne/Bureau/UNIL_CEE/Project/Graph and models/Minimal_Law_rich_poor/boxplot_rich_poor_nonavoid.csv"
)

read_box <- function(path) as.data.frame(fread(path))

df_low <- read_box(paths_box$low_var)
df_unc <- read_box(paths_box$uncorr)
df_to  <- read_box(paths_box$tradeoff)
df_rp  <- read_box(paths_box$richpoor)

base_colors <- c(
  war      = "firebrick",
  raid     = "darkorange3",
  peaceful = "darkgreen"
)

cost_levels <- c(2.5, 8, 14)
cost_labels <- c(
  "2.5" = "(a) Low cost",
  "8"   = "(b) Intermediate cost",
  "14"  = "(c) High cost"
)

prep_actions_allcosts <- function(df, model_label) {
  df %>%
    filter(Cw %in% cost_levels) %>%
    mutate(
      peaceful = trade + inaction
    ) %>%
    select(run_id, Cw, war, raid, peaceful) %>%
    pivot_longer(
      cols = c(war, raid, peaceful),
      names_to = "action",
      values_to = "value"
    ) %>%
    mutate(
      model = model_label,
      action = factor(action, levels = names(base_colors)),
      cost = factor(Cw, levels = cost_levels, labels = cost_labels[as.character(cost_levels)])
    )
}

plot_box_allcosts <- function(df_low, df_var, var_name) {
  
  long_var <- prep_actions_allcosts(df_var, var_name)
  long_low <- prep_actions_allcosts(df_low, "Low variability")
  
  action_levels <- levels(long_var$action)
  nudge <- 0.18
  
  # Var on the RIGHT, low-var on the LEFT (swap +/- if you ever want the opposite)
  long_var <- long_var %>%
    mutate(
      x_base = as.numeric(factor(action, levels = action_levels)),
      x_pos  = x_base + nudge
    )
  
  long_low <- long_low %>%
    mutate(
      x_base = as.numeric(factor(action, levels = action_levels)),
      x_pos  = x_base - nudge
    )
  
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
      size = 1.5,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = base_colors, guide = "none") +
    scale_color_manual(
      name = var_name,
      values = base_colors,
      breaks = names(base_colors),
      labels = c(
        war      = "War",
        raid     = "Raid",
        peaceful = "Peaceful"
      )
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
      size = 1.5,
      show.legend = FALSE
    ) +
    scale_color_manual(
      name = "Low variability",
      values = scales::alpha(base_colors, 0.6),
      breaks = names(base_colors),
      labels = c(
        war      = "War",
        raid     = "Raid",
        peaceful = "Peaceful"
      )
    ) +
    facet_wrap(~ cost, nrow = 1) +
    scale_x_continuous(
      breaks = seq_along(action_levels),
      labels = c(war = "war", raid = "raid", peaceful = "peaceful")
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      x = "",
      y = "Proportion (%)",
      title = paste0("Action proportions: Low variability vs ", var_name)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),
      legend.position = "right",
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 12),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
      strip.text = element_text(size = 15),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 15)
    )
}


# Example calls:
# p_unc <- plot_box_allcosts(df_low, df_unc, "Uncorrelated variability")
# p_to  <- plot_box_allcosts(df_low, df_to,  "Trade-off variability")
# p_rp  <- plot_box_allcosts(df_low, df_rp,  "Rich-poor variability")
# print(p_unc); print(p_to); print(p_rp)


#------------------ graph method----------
library(dplyr)
library(ggplot2)
library(truncnorm)
library(ggExtra)
library(cowplot)

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
