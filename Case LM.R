# =========================================================
# Classical Gaussian prediction intervals
# Linear regression benchmark
# =========================================================

library(dplyr)
library(ggplot2)

# =========================================================
# Stratified split
# =========================================================

cph_make_split <- function(df) {
  
  df %>%
    
    group_by(cph_salgsaar) %>%
    
    mutate(
      
      cph_n = n(),
      
      cph_idx = sample.int(cph_n[1]),
      
      split = case_when(
        
        cph_idx <= round(0.6 * cph_n)
        ~ "train",
        
        cph_idx <= round(0.8 * cph_n)
        ~ "cal",
        
        TRUE
        ~ "test"
      )
    ) %>%
    
    ungroup()
}

set.seed(123)

cph_split_r <- cph_make_split(cph_df)

cph_train_r <- cph_split_r %>%
  filter(split == "train") %>%
  select(-split, -cph_n, -cph_idx)

cph_test_r <- cph_split_r %>%
  filter(split == "test") %>%
  select(-split, -cph_n, -cph_idx)
# =========================================================
# Fit linear model
# =========================================================

lm_model <- lm(
  log_KOEBESUM ~ .,
  data = cph_train_r
)

# =========================================================
# Classical prediction intervals
# =========================================================

lm_pred <- predict(
  lm_model,
  newdata = cph_test_r,
  
  interval = "prediction",
  level = 0.90
)

# =========================================================
# Construct dataframe
# =========================================================

lm_analysis_df <- data.frame(
  
  faktisk = exp(cph_test_r$log_KOEBESUM),
  
  pred = exp(lm_pred[, "fit"]),
  
  lower = exp(lm_pred[, "lwr"]),
  
  upper = exp(lm_pred[, "upr"]),
  
  covered =
    cph_test_r$log_KOEBESUM >= lm_pred[, "lwr"] &
    cph_test_r$log_KOEBESUM <= lm_pred[, "upr"],
  
  width =
    exp(lm_pred[, "upr"]) -
    exp(lm_pred[, "lwr"]),
  
  rel_width =
    (
      exp(lm_pred[, "upr"]) -
        exp(lm_pred[, "lwr"])
    ) / exp(cph_test_r$log_KOEBESUM),
  
  cph_salgsaar = cph_test_r$cph_salgsaar,
  
  KOMMUNE_NAVN = cph_test_r$KOMMUNE_NAVN
)

# =========================================================
# Overall metrics
# =========================================================

cat("\n============================\n")
cat("Linear model benchmark\n")
cat("============================\n\n")

cat("Coverage:\n")
print(mean(lm_analysis_df$covered))

cat("\nMean width:\n")
print(mean(lm_analysis_df$width))

cat("\nMean relative width:\n")
print(mean(lm_analysis_df$rel_width))

# =========================================================
# Price bins
# =========================================================

breaks_price <- quantile(
  lm_analysis_df$faktisk / 1e6,
  probs = seq(0,1,length.out = 11),
  na.rm = TRUE
)

lm_analysis_df <- lm_analysis_df %>%
  mutate(
    
    price_bin = cut(
      faktisk / 1e6,
      breaks = breaks_price,
      include.lowest = TRUE
    )
  )

# =========================================================
# Coverage by year
# =========================================================

coverage_year_lm <- lm_analysis_df %>%
  group_by(cph_salgsaar) %>%
  summarise(
    coverage = mean(covered),
    .groups = "drop"
  )

p_cov_year_lm <- ggplot(
  coverage_year_lm,
  
  aes(
    x = cph_salgsaar,
    y = coverage
  )
) +
  
  geom_line() +
  geom_point(size = 2) +
  
  geom_hline(
    yintercept = 0.90,
    color = "red",
    linetype = "dashed"
  ) +
  
  ylim(0.6,1) +
  
  theme_minimal() +
  
  labs(
    title = "Gaussian prediction intervals: coverage by year",
    x = "Sale year",
    y = "Coverage"
  )

p_cov_year_lm

# =========================================================
# Coverage by price
# =========================================================

coverage_price_lm <- lm_analysis_df %>%
  group_by(price_bin) %>%
  summarise(
    coverage = mean(covered),
    .groups = "drop"
  )

p_cov_price_lm <- ggplot(
  coverage_price_lm,
  
  aes(
    x = price_bin,
    y = coverage,
    group = 1
  )
) +
  
  geom_line() +
  geom_point(size = 2) +
  
  geom_hline(
    yintercept = 0.90,
    color = "red",
    linetype = "dashed"
  ) +
  
  ylim(0.6,1) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    title = "Gaussian prediction intervals: coverage by price",
    x = "Price interval (million DKK)",
    y = "Coverage"
  )

p_cov_price_lm

# =========================================================
# Coverage by municipality
# =========================================================

coverage_area_lm <- lm_analysis_df %>%
  group_by(KOMMUNE_NAVN) %>%
  summarise(
    coverage = mean(covered),
    n = n(),
    .groups = "drop"
  ) %>%
  
  filter(n > 50)

p_cov_area_lm <- ggplot(
  coverage_area_lm,
  
  aes(
    x = reorder(KOMMUNE_NAVN, coverage),
    y = coverage
  )
) +
  
  geom_point(size = 2) +
  
  geom_hline(
    yintercept = 0.90,
    color = "red",
    linetype = "dashed"
  ) +
  
  coord_flip() +
  
  ylim(0.6,1) +
  
  theme_minimal() +
  
  labs(
    title = "Gaussian prediction intervals: coverage by municipality",
    x = NULL,
    y = "Coverage"
  )

p_cov_area_lm