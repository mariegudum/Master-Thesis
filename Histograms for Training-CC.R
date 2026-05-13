# =============================
# Training-conditional coverage
# Histogram across simulation runs
# =============================

library(dplyr)
library(ggplot2)

alpha <- 0.1

# Keep only CP methods
df_tc <- results_all_4 %>%
  filter(method %in% c("CP_LM", "CP_RF_tuned"))

# =============================
# Beta density helper
# =============================
beta_df <- expand.grid(
  x = seq(0.75, 1.0, length.out = 1000),
  n = unique(df_tc$n),
  method = unique(df_tc$method)
) %>%
  mutate(
    a = (1 - alpha) * (n/2 + 1),   # calibration size = n/2
    b = alpha * (n/2 + 1),
    density = dbeta(x, a, b)
  )

# =============================
# Histogram + Beta overlay
# =============================
ggplot(df_tc, aes(x = coverage)) +
  
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 20,
    fill = "skyblue",
    color = "white",
    alpha = 0.7
  ) +
  
  geom_line(
    data = beta_df,
    aes(x = x, y = density),
    color = "red",
    linewidth = 1
  ) +
  
  facet_grid(method ~ n, scales = "free_y") +
  
  geom_vline(
    xintercept = 0.9,
    linetype = "dashed"
  ) +
  
  labs(
    title = "Training-conditional coverage across simulation runs",
    subtitle = "Histogram of empirical coverage with Beta approximation overlay",
    x = "Empirical coverage",
    y = "Density"
  ) +
  
  theme_minimal(base_size = 13)