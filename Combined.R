results_extended_new_1$experiment = "Linear Gaussian"
results_all_2$experiment = "Linear Asymmetric"
results_all_3$experiment = "Non-linear Gaussian"
results_all_4$experiment = "Heteroskedastic Asymmetric"

df_all = rbind(results_extended_new_1,
               results_all_2,
               results_all_3,
               results_all_4)

df_all = df_all %>%
  mutate(
    method = recode(method,
                    "CP_LM" = "CP (LM)",
                    "CP_RF_tuned" = "CP (RF)",
                    "Gaussian_LM" = "Gaussian (LM)",
                    "QR_LM" = "QR (LM)",
                    "QR_RF_tuned" = "QR (RF)",
                    "RQ_LM" = "RQ (LM)",
                    "RQ_RF_tuned" = "RQ (RF)"
    )
  )

df_plot = df_all %>%
  group_by(experiment, method, n) %>%
  summarise(
    coverage = mean(coverage),
    length = mean(length),
    .groups = "drop"
  )

ggplot(df_plot,
       aes(x = coverage, y = length,
           color = method, shape = factor(n))) +
  
  geom_point(size = 3) +
  geom_path(aes(group = method), linewidth = 0.5) +
  
  geom_vline(xintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, scales = "free_y") +
  
  labs(
    x = "Empirical coverage",
    y = "Relative interval length",
    color = "Method",
    shape = "n"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )



ggplot(df_plot,
       aes(x = coverage, y = length,
           color = method)) +
  
  geom_point(size = 2) +
  geom_path(aes(group = method), linewidth = 0.8) +
  
  geom_vline(xintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, scales = "free_y") +
  
  theme_minimal()


ggplot(df_plot,
       aes(x = coverage, y = length,
           color = n, group = interaction(method))) +
  
  geom_point(size = 2) +
  geom_path(aes(color = method), linewidth = 0.8) +
  
  scale_color_viridis_c() +
  
  facet_wrap(~ experiment) +
  theme_minimal()


ggplot(df_plot,
       aes(x = coverage, y = length,
           color = method, group = method)) +
  
  geom_point(size = 2) +
  geom_path(linewidth = 0.7) +
  
  geom_vline(xintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment) +
  
  theme_minimal()

# alpha grading
ggplot(df_plot,
       aes(x = coverage, y = length,
           color = method, alpha = n, group = method)) +
  
  geom_point(size = 2) +
  geom_path(linewidth = 0.7) +
  
  scale_alpha(range = c(0.4, 1)) +
  
  facet_wrap(~ experiment) +
  theme_minimal()


# grading

library(dplyr)

# vælg basefarver (én pr metode)
base_colors = c(
  "CP_LM" = "red",
  "CP_RF_tuned" = "pink",
  "Gaussian_LM" = "green",
  "QR_LM" = "blue",
  "QR_RF_tuned" = "purple",
  "RQ_LM" = "yellow",
  "RQ_RF_tuned" = "orange"
)

library(colorspace)
df_plot = df_plot %>%
  group_by(method) %>%
  mutate(
    n_scaled = (n - min(n)) / (max(n) - min(n)),
    color = colorspace::lighten(base_colors[method], amount = 0.5 * (1 - n_scaled))
  ) %>%
  ungroup()

ggplot(df_plot,
       aes(x = coverage, y = length,
           group = method)) +
  
  geom_path(aes(color = color), linewidth = 0.3) +
  geom_point(aes(color = color), size = 1) +
  
  geom_vline(xintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, scales = "free_y") +
  
  scale_color_identity() +
  
  theme_minimal() +
  theme(legend.position = "none")

ggplot(df_plot,
       aes(x = coverage, y = length,
           color = method, alpha = n, group = method)) +
  
  geom_point(size = 1) +
  geom_path(linewidth = 0.4) +
  
  scale_alpha(range = c(0.4, 1)) +
  
  facet_wrap(~ experiment) +
  theme_minimal()





# kun laveste og højeste n med length ratio

df_plot = df_plot %>%
  mutate(
    method = recode(method,
                    "QR_RF_tuned_per_n" = "QR (RF)")
  )


df_plot <- df_plot %>%
  filter(method != "QR_RF") %>%
  filter(method != "CP_RF") %>%
  filter(method != "RQ_RF") 

df_plot = df_plot %>%
  group_by(experiment, n) %>%
  mutate(
    length_oracle = length[method == "Oracle"],
    length_ratio = length / length_oracle
  ) %>%
  ungroup()

df_plot <- df_plot %>%
  filter(method != "Oracle")
head(df_plot)

df_plot$method = factor(df_plot$method,
                       levels = c("Gaussian (LM)", "CP (LM)", "CP (RF)",
                                  "QR (LM)", "QR (RF)",
                                  "RQ (LM)", "RQ (RF)")
)

df_extreme = df_plot %>%
  group_by(experiment, method) %>%
  filter(n == min(n) | n == max(n)) %>%
  ungroup()

df_extreme$experiment = factor(df_extreme$experiment,
                            levels = c(
                              "Linear Gaussian",
                              "Linear Asymmetric",
                              "Non-linear Gaussian",
                              "Heteroskedastic Asymmetric"
                            )
)

ggplot(df_extreme,
       aes(x = coverage, y = length_ratio,
           color = method, group = method)) +
  
  geom_point(size = 1.5, alpha = 0.7) +
  geom_path(linewidth = 0.5) +
  
  geom_vline(xintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, scales = "free_y") +
  
  theme_minimal(base_size = 13)

# byttet om på akser

ggplot(df_plot,
       aes(x = length_ratio, y = coverage,
           color = method, group = method)) +
  
  geom_point(size = 1.5, alpha = 0.7) +
  geom_path(aes(group = method),
            arrow = arrow(length = unit(0.2, "cm"))) +
  
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, scales = "free_y") +
  
  labs(
    x = "Relative interval length (normalized by Oracle)",
    y = "Empirical coverage",
    color = "Method",
    title = "Coverage vs. Relative Interval Length Across Sample Sizes",
    subtitle = "Arrows run from lowest sample size to highest sample size"
  ) +
  
  theme_minimal(base_size = 13)


# legeplads

df_plot$experiment = factor(df_plot$experiment,
                               levels = c(
                                 "Linear Gaussian",
                                 "Linear Asymmetric",
                                 "Non-linear Gaussian",
                                 "Heteroskedastic Asymmetric"
                               )
)

ggplot(df_plot,
       aes(x = length_ratio, y = coverage,
           color = method, group = method)) +
  geom_point(aes(size = n)) +
  scale_size_continuous(range = c(1, 2)) +
  geom_path(aes(group = method),
            arrow = arrow(length = unit(0.2, "cm"))) +
  
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, scales = "free_y") +
  
  theme_minimal(base_size = 13)


ggplot(df_plot,
       aes(x = length_ratio, y = coverage,
           color = method, group = method)) +
  
  geom_path(linewidth = 0.8) +   # forbinder i rækkefølge
  
  geom_point(size = 2.5) +
  
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, ncol = 2) +
  
  theme_minimal(base_size = 13) +
  
  geom_text(aes(label = n_id), vjust = -0.8)



df_plot <- df_plot %>%
  group_by(experiment, method) %>%
  arrange(n) %>%
  mutate(n_id = row_number()) %>%
  ungroup()

ggplot(df_plot,
       aes(x = length_ratio, y = coverage,
           color = method)) +
  
  geom_text(aes(label = n_id), size = 5) +
  
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, ncol = 2) +
  
  labs(
    x = "Relative interval length (normalized by Oracle)",
    y = "Empirical coverage",
    color = "Method"
  ) +
  
  theme_minimal(base_size = 13) +
  labs(
    title = "Coverage vs. Relative Interval Length Across Sample Sizes",
    subtitle = "Numbers indicate sample sizes (1 = smallest, 4 = largest)"
  )

# linje gennem

df_arrow <- df_plot %>%
  group_by(experiment, method) %>%
  arrange(n_id) %>%
  summarise(
    x_start = first(length_ratio),
    y_start = first(coverage),
    x_end   = last(length_ratio),
    y_end   = last(coverage),
    .groups = "drop"
  )

ggplot(df_plot,
       aes(x = length_ratio, y = coverage,
           color = method)) +
  
  geom_text(aes(label = n_id), size = 5) +
  
  geom_segment(
    data = df_arrow,
    aes(x = x_start, y = y_start,
        xend = x_end, yend = y_end,
        color = method),
    arrow = arrow(length = unit(0.15, "cm")),
    linewidth = 0.5,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, ncol = 2) +
  
  labs(
    x = "Relative interval length (normalized by Oracle)",
    y = "Empirical coverage",
    color = "Method",
    title = "Coverage vs. Relative Interval Length Across Sample Sizes",
    subtitle = "Numbers indicate sample sizes (1 = smallest, 4 = largest)"
  ) +
  
  theme_minimal(base_size = 13)



# lm linje


ggplot(df_plot,
       aes(x = length_ratio, y = coverage,
           color = method)) +
  
  geom_text(aes(label = n_id), size = 5) +
  
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, ncol = 2) +
  
  labs(
    x = "Relative interval length (normalized by Oracle)",
    y = "Empirical coverage",
    color = "Method",
    title = "Coverage vs. Relative Interval Length Across Sample Sizes",
    subtitle = "Numbers indicate sample sizes (1 = smallest, 4 = largest)"
  ) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", ) +
  
  theme_minimal(base_size = 13)

# final

ggplot(df_plot,
       aes(x = length_ratio, y = coverage,
           color = method)) +
  
  geom_text(aes(label = n_id), size = 5) +
  
  geom_segment(
    data = df_arrow,
    aes(x = x_start, y = y_start,
        xend = x_end, yend = y_end,
        color = method),
    linewidth = 0.5,
    linetype = "dashed",
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  
  facet_wrap(~ experiment, ncol = 2) +
  
  labs(
    x = "Relative interval length (normalized by Oracle)",
    y = "Empirical coverage",
    color = "Method",
    title = "Coverage vs. Relative Interval Length Across Sample Sizes",
    subtitle = "Numbers indicate sample sizes (1 = smallest, 4 = largest)"
  ) +
  
  theme_minimal(base_size = 13)


