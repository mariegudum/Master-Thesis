# =========================================================
# MONDRIAN / LOCALIZED CP
# =========================================================
# Variants:
# 1. Global CP
# 2. Municipality-bin CP
# 3. Year-bin CP
# 4. Prediction-bin CP
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
    log_KOEBESUM = log(KONTANT_KOEBESUM),
    
    cph_salgsaar = as.integer(format(SALG_KVARTAL, "%Y")),
    
    ENH_BAD_KODE = as.factor(ENH_BAD_KODE),
    LANDSDEL_NAVN = as.factor(LANDSDEL_NAVN),
    KOMMUNE_NAVN = as.factor(KOMMUNE_NAVN)
  ) %>%
  select(-KONTANT_KOEBESUM, -SALG_KVARTAL)

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
        cph_idx <= round(0.6 * cph_n) ~ "train",
        cph_idx <= round(0.8 * cph_n) ~ "cal",
        TRUE ~ "test"
      )
    ) %>%
    ungroup()
}

# =========================================================
# Storage
# =========================================================

global_results      <- list()
kommune_results     <- list()
year_results        <- list()
prediction_results  <- list()

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
  # Fit RF
  # =======================================================
  
  cph_rf <- ranger(
    formula = log_KOEBESUM ~ .,
    data = cph_train,
    
    num.trees = 500,
    
    num.threads = parallel::detectCores() - 1,
    
    seed = sample(1:10000, 1)
  )
  
  # =======================================================
  # Predictions
  # =======================================================
  
  cal_pred <- predict(cph_rf, data = cph_cal)$predictions
  
  test_pred <- predict(cph_rf, data = cph_test)$predictions
  
  # =======================================================
  # Residuals
  # =======================================================
  
  cph_cal$residual <- abs(
    cph_cal$log_KOEBESUM - cal_pred
  )
  
  # =======================================================
  # GLOBAL CP
  # =======================================================
  
  global_qhat <- quantile(
    cph_cal$residual,
    probs = ceiling(
      (nrow(cph_cal)+1)*(1-cph_alpha)
    ) / nrow(cph_cal)
  )
  
  global_df <- data.frame(
    
    faktisk = exp(cph_test$log_KOEBESUM),
    
    lower = exp(test_pred - global_qhat),
    upper = exp(test_pred + global_qhat),
    
    covered =
      cph_test$log_KOEBESUM >= (test_pred - global_qhat) &
      cph_test$log_KOEBESUM <= (test_pred + global_qhat),
    
    width =
      exp(test_pred + global_qhat) -
      exp(test_pred - global_qhat),
    
    method = "Global CP",
    
    cph_salgsaar = cph_test$cph_salgsaar,
    KOMMUNE_NAVN = cph_test$KOMMUNE_NAVN,
    
    replication = i
  )
  
  global_results[[i]] <- global_df
  
  # =======================================================
  # 1. MONDRIAN CP - MUNICIPALITY BIN
  # =======================================================
  
  kommune_prices <- cph_cal %>%
    group_by(KOMMUNE_NAVN) %>%
    summarise(
      mean_price = mean(exp(log_KOEBESUM)),
      .groups = "drop"
    )
  
  kommune_prices <- kommune_prices %>%
    mutate(
      kommune_bin = cut(
        mean_price,
        
        breaks = quantile(
          mean_price,
          probs = seq(0,1,length.out = 4),
          na.rm = TRUE
        ),
        
        labels = c(
          "Low",
          "Medium",
          "High"
        ),
        
        include.lowest = TRUE
      )
    )
  
  cph_cal <- cph_cal %>%
    left_join(
      kommune_prices %>%
        select(KOMMUNE_NAVN, kommune_bin),
      by = "KOMMUNE_NAVN"
    )
  
  cph_test <- cph_test %>%
    left_join(
      kommune_prices %>%
        select(KOMMUNE_NAVN, kommune_bin),
      by = "KOMMUNE_NAVN"
    )
  
  kommune_qhat <- cph_cal %>%
    group_by(kommune_bin) %>%
    summarise(
      qhat = quantile(
        residual,
        probs = ceiling(
          (n()+1)*(1-cph_alpha)
        ) / n()
      ),
      .groups = "drop"
    )
  
  cph_test_kommune <- cph_test %>%
    left_join(
      kommune_qhat,
      by = "kommune_bin"
    )
  
  kommune_df <- data.frame(
    
    faktisk = exp(cph_test_kommune$log_KOEBESUM),
    
    lower = exp(test_pred - cph_test_kommune$qhat),
    upper = exp(test_pred + cph_test_kommune$qhat),
    
    covered =
      cph_test_kommune$log_KOEBESUM >=
      (test_pred - cph_test_kommune$qhat) &
      
      cph_test_kommune$log_KOEBESUM <=
      (test_pred + cph_test_kommune$qhat),
    
    width =
      exp(test_pred + cph_test_kommune$qhat) -
      exp(test_pred - cph_test_kommune$qhat),
    
    method = "Mondrian CP - Municipality",
    
    kommune_bin = cph_test_kommune$kommune_bin,
    
    replication = i
  )
  
  kommune_results[[i]] <- kommune_df
  
  # =======================================================
  # 2. MONDRIAN CP - YEAR
  # =======================================================
  
  year_qhat <- cph_cal %>%
    group_by(cph_salgsaar) %>%
    summarise(
      qhat = quantile(
        residual,
        probs = ceiling(
          (n()+1)*(1-cph_alpha)
        ) / n()
      ),
      .groups = "drop"
    )
  
  cph_test_year <- cph_test %>%
    left_join(
      year_qhat,
      by = "cph_salgsaar"
    )
  
  year_df <- data.frame(
    
    faktisk = exp(cph_test_year$log_KOEBESUM),
    
    lower = exp(test_pred - cph_test_year$qhat),
    upper = exp(test_pred + cph_test_year$qhat),
    
    covered =
      cph_test_year$log_KOEBESUM >=
      (test_pred - cph_test_year$qhat) &
      
      cph_test_year$log_KOEBESUM <=
      (test_pred + cph_test_year$qhat),
    
    width =
      exp(test_pred + cph_test_year$qhat) -
      exp(test_pred - cph_test_year$qhat),
    
    method = "Mondrian CP - Year",
    
    cph_salgsaar = cph_test_year$cph_salgsaar,
    
    replication = i
  )
  
  year_results[[i]] <- year_df
  
  # =======================================================
  # 3. MONDRIAN CP - PREDICTION BINS
  # =======================================================
  
  cph_cal$predicted_value <- cal_pred
  
  pred_breaks <- quantile(
    cal_pred,
    probs = seq(0,1,length.out = 4),
    na.rm = TRUE
  )
  
  cph_cal$pred_bin <- cut(
    cal_pred,
    breaks = pred_breaks,
    labels = c("Low","Medium","High"),
    include.lowest = TRUE
  )
  
  cph_test$pred_bin <- cut(
    test_pred,
    breaks = pred_breaks,
    labels = c("Low","Medium","High"),
    include.lowest = TRUE
  )
  
  pred_qhat <- cph_cal %>%
    group_by(pred_bin) %>%
    summarise(
      qhat = quantile(
        residual,
        probs = ceiling(
          (n()+1)*(1-cph_alpha)
        ) / n()
      ),
      .groups = "drop"
    )
  
  cph_test_pred <- cph_test %>%
    left_join(
      pred_qhat,
      by = "pred_bin"
    )
  
  pred_df <- data.frame(
    
    faktisk = exp(cph_test_pred$log_KOEBESUM),
    
    lower = exp(test_pred - cph_test_pred$qhat),
    upper = exp(test_pred + cph_test_pred$qhat),
    
    covered =
      cph_test_pred$log_KOEBESUM >=
      (test_pred - cph_test_pred$qhat) &
      
      cph_test_pred$log_KOEBESUM <=
      (test_pred + cph_test_pred$qhat),
    
    width =
      exp(test_pred + cph_test_pred$qhat) -
      exp(test_pred - cph_test_pred$qhat),
    
    method = "Mondrian CP - Prediction",
    
    pred_bin = cph_test_pred$pred_bin,
    
    replication = i
  )
  
  prediction_results[[i]] <- pred_df
}

# =========================================================
# Merge results
# =========================================================

global_df_all <- bind_rows(global_results)

kommune_df_all <- bind_rows(kommune_results)

year_df_all <- bind_rows(year_results)

prediction_df_all <- bind_rows(prediction_results)

all_cp_results <- bind_rows(
  global_df_all,
  kommune_df_all,
  year_df_all,
  prediction_df_all
)

# =========================================================
# Add metrics
# =========================================================

all_cp_results <- all_cp_results %>%
  mutate(
    rel_width = width / faktisk
  )

# =========================================================
# Summary table
# =========================================================

summary_table <- all_cp_results %>%
  group_by(method) %>%
  summarise(
    coverage = mean(covered),
    mean_width = mean(width),
    mean_rel_width = mean(rel_width),
    .groups = "drop"
  )

print(summary_table)

# =========================================================
# Coverage-width plot
# =========================================================

ggplot(
  summary_table,
  
  aes(
    x = mean_width,
    y = coverage,
    color = method
  )
) +
  
  geom_point(size = 4) +
  
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Coverage vs width",
    x = "Mean interval width",
    y = "Coverage"
  )