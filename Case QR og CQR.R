# =========================================================
# Quantile Regression Forest (QR) + CQR with ranger
# =========================================================

library(ranger)
library(dplyr)
library(ggplot2)

# =========================================================
# Parameters
# =========================================================

cph_alpha <- 0.10
cph_n_rep <- 100

# =========================================================
# Data preparation
# =========================================================

cph_df <- CphHousingPrices %>%
  na.omit() %>%
  mutate(
    log_KOEBESUM  = log(KONTANT_KOEBESUM),
    cph_salgsaar  = as.integer(format(SALG_KVARTAL, "%Y")),
    
    ENH_BAD_KODE  = as.factor(ENH_BAD_KODE),
    LANDSDEL_NAVN = as.factor(LANDSDEL_NAVN),
    KOMMUNE_NAVN  = as.factor(KOMMUNE_NAVN)
  ) %>%
  select(-KONTANT_KOEBESUM, -SALG_KVARTAL)

# =========================================================
# Stratified split
# =========================================================

cph_make_split <- function(df) {
  
  df %>%
    group_by(cph_salgsaar) %>%
    mutate(
      cph_n   = n(),
      cph_idx = sample.int(cph_n[1]),
      
      split = case_when(
        cph_idx <= round(0.6 * cph_n) ~ "train",
        cph_idx <= round(0.8 * cph_n) ~ "cal",
        TRUE                          ~ "test"
      )
    ) %>%
    ungroup()
}

# =========================================================
# Storage
# =========================================================

qr_results  <- list()
cqr_results <- list()

# =========================================================
# Main loop
# =========================================================

for(i in 1:cph_n_rep) {
  
  cat("Replication:", i, "\n")
  
  # =======================================================
  # Split data
  # =======================================================
  
  cph_split <- cph_make_split(cph_df)
  
  cph_train <- cph_split %>%
    filter(split == "train") %>%
    select(-split, -cph_idx, -cph_n)
  
  cph_cal <- cph_split %>%
    filter(split == "cal") %>%
    select(-split, -cph_idx, -cph_n)
  
  cph_test <- cph_split %>%
    filter(split == "test") %>%
    select(-split, -cph_idx, -cph_n)
  
  # =======================================================
  # Train quantile forest
  # =======================================================
  
  qr_model <- ranger(
    formula = log_KOEBESUM ~ .,
    data = cph_train,
    
    num.trees = 500,
    quantreg = TRUE,
    
    num.threads = parallel::detectCores() - 1,
    
    seed = sample(1:10000, 1)
  )
  
  # =======================================================
  # QR predictions - calibration
  # =======================================================
  
  cal_pred <- predict(
    qr_model,
    data = cph_cal,
    type = "quantiles",
    quantiles = c(cph_alpha/2, 1 - cph_alpha/2)
  )$predictions
  
  cal_lower <- cal_pred[,1]
  cal_upper <- cal_pred[,2]
  
  # =======================================================
  # CQR conformity scores
  # =======================================================
  
  cqr_scores <- pmax(
    cal_lower - cph_cal$log_KOEBESUM,
    cph_cal$log_KOEBESUM - cal_upper
  )
  
  n_cal <- length(cqr_scores)
  
  qhat <- quantile(
    cqr_scores,
    probs = ceiling((n_cal + 1) * (1 - cph_alpha)) / n_cal
  )
  
  # =======================================================
  # Test predictions
  # =======================================================
  
  test_pred <- predict(
    qr_model,
    data = cph_test,
    type = "quantiles",
    quantiles = c(cph_alpha/2, 1 - cph_alpha/2)
  )$predictions
  
  qr_lower <- test_pred[,1]
  qr_upper <- test_pred[,2]
  
  # =======================================================
  # QR intervals
  # =======================================================
  
  qr_df <- data.frame(
    
    faktisk = exp(cph_test$log_KOEBESUM),
    
    lower = exp(qr_lower),
    upper = exp(qr_upper),
    
    covered =
      cph_test$log_KOEBESUM >= qr_lower &
      cph_test$log_KOEBESUM <= qr_upper,
    
    width = exp(qr_upper) - exp(qr_lower),
    
    method = "QR",
    
    cph_salgsaar = cph_test$cph_salgsaar,
    KOMMUNE_NAVN = cph_test$KOMMUNE_NAVN,
    LANDSDEL_NAVN = cph_test$LANDSDEL_NAVN,
    ENH_BEBO_ARL = cph_test$ENH_BEBO_ARL,
    
    replication = i
  )
  
  # =======================================================
  # CQR intervals
  # =======================================================
  
  cqr_lower <- qr_lower - qhat
  cqr_upper <- qr_upper + qhat
  
  cqr_df <- data.frame(
    
    faktisk = exp(cph_test$log_KOEBESUM),
    
    lower = exp(cqr_lower),
    upper = exp(cqr_upper),
    
    covered =
      cph_test$log_KOEBESUM >= cqr_lower &
      cph_test$log_KOEBESUM <= cqr_upper,
    
    width = exp(cqr_upper) - exp(cqr_lower),
    
    method = "CQR",
    
    cph_salgsaar = cph_test$cph_salgsaar,
    KOMMUNE_NAVN = cph_test$KOMMUNE_NAVN,
    LANDSDEL_NAVN = cph_test$LANDSDEL_NAVN,
    ENH_BEBO_ARL = cph_test$ENH_BEBO_ARL,
    
    replication = i
  )
  
  # =======================================================
  # Store
  # =======================================================
  
  qr_results[[i]]  <- qr_df
  cqr_results[[i]] <- cqr_df
}

# =========================================================
# Merge all replications
# =========================================================

qr_analysis_df <- bind_rows(qr_results)

cqr_analysis_df <- bind_rows(cqr_results)

# =========================================================
# Summary
# =========================================================

qr_analysis_df %>%
  summarise(
    coverage = mean(covered),
    width = mean(width)
  )

cqr_analysis_df %>%
  summarise(
    coverage = mean(covered),
    width = mean(width)
  )

# =========================================================
# Example:
# Coverage by year
# =========================================================

coverage_year <- cqr_analysis_df %>%
  group_by(cph_salgsaar) %>%
  summarise(
    coverage = mean(covered),
    width = mean(width),
    .groups = "drop"
  )

ggplot(coverage_year,
       aes(cph_salgsaar, coverage)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.90,
             linetype = "dashed",
             color = "red") +
  theme_minimal()

# =========================================================
# Example:
# Coverage by price bins
# =========================================================

cqr_analysis_df <- cqr_analysis_df %>%
  mutate(
    price_bin = cut(
      faktisk,
      breaks = quantile(
        faktisk,
        probs = seq(0,1,length.out = 11)
      ),
      include.lowest = TRUE
    )
  )

coverage_price <- cqr_analysis_df %>%
  group_by(price_bin) %>%
  summarise(
    coverage = mean(covered),
    width = mean(width),
    .groups = "drop"
  )

ggplot(coverage_price,
       aes(price_bin,
           coverage,
           group = 1)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.90,
             linetype = "dashed",
             color = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )