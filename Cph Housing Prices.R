library(ggplot2)
library(dplyr)
library(ranger)
library(missRanger)

View(CphHousingPrices)
CphHousingPrices <- Region_H_2000_og_frem_boligdata_Korrekt[, 1:21]
summary(CphHousingPrices)
hist(CphHousingPrices$BYG_ALDER_AAR)

# Missing values
colSums(is.na(cph_df))

# Boligprisudvikling over tid (per kvartal)

CphHousingPrices %>%
  group_by(SALG_KVARTAL) %>%
  summarise(mean_pris = mean(KONTANT_KOEBESUM, na.rm = TRUE)) %>%
  ggplot(aes(x = SALG_KVARTAL, y = mean_pris)) +
  geom_line() +
  geom_smooth(method = "loess", se = TRUE, color = "steelblue") +
  labs(title = "Meanboligpris over tid",
       x = "Kvartal",
       y = "Mean kontantpris (DKK)") +
  theme_minimal()

# Log-transformeret fordeling
CphHousingPrices$log_pris <- log(CphHousingPrices$KONTANT_KOEBESUM)

ggplot(CphHousingPrices, aes(x = log_pris)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  labs(title = "Log-transformeret boligpris",
       x = "log(Kontantpris)",
       y = "Frekvens") +
  theme_minimal()




# Mean pris per postnummer
CphHousingPrices %>%
  group_by(KOMMUNE_NAVN) %>%
  summarise(mean_pris = mean(KONTANT_KOEBESUM, na.rm = TRUE),
            antal = n()) %>%
  filter(antal > 30) %>%
  arrange(desc(mean_pris)) %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(KOMMUNE_NAVN, mean_pris), y = mean_pris)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 dyreste postnumre",
       x = "Postnummer",
       y = "Mean kontantpris (DKK)") +
  theme_minimal()

library(corrplot)

# Vælg numeriske variable
num_vars <- CphHousingPrices %>%
  select(where(is.numeric)) %>%
  select(-log_pris)  # undgå redundans

# Korrelationsmatrix
cor_matrix <- cor(num_vars, use = "complete.obs")

corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         tl.cex = 0.7,
         addCoef.col = "black",
         number.cex = 0.6,
         title = "Korrelationsmatrix")


# Scatterplot: pris vs. størrelse
ggplot(CphHousingPrices, aes(x = ENH_BEBO_ARL, y = KONTANT_KOEBESUM)) +
  geom_point(alpha = 0.1, color = "steelblue") +
  geom_smooth(method = "lm", color = "red") +
  scale_y_log10() +
  scale_x_log10() +
  labs(title = "Boligpris vs. størrelse (log-log)",
       x = "Boligareal (m²)",
       y = "Kontantpris (DKK)") +
  theme_minimal()







# ── Pakker ────────────────────────────────────────────────────────────────────
library(ranger)
library(dplyr)
library(ggplot2)

# ── Data forberedelse ─────────────────────────────────────────────────────────
cph_df <- CphHousingPrices %>%
  select(KONTANT_KOEBESUM, ENH_BEBO_ARL, ENH_VAERELSE_ANT,
         BYG_ALDER_AAR, kyst_AFSTAND, ENH_ANVEND_KODE,
         SALG_KVARTAL) %>%
  na.omit() %>%
  mutate(SALG_KVARTAL = as.Date(SALG_KVARTAL))

# ── Tidsbaseret split ─────────────────────────────────────────────────────────
cph_train <- cph_df %>% filter(SALG_KVARTAL <= as.Date("2019-12-31"))
cph_cal   <- cph_df %>% filter(SALG_KVARTAL >  as.Date("2019-12-31") & 
                                 SALG_KVARTAL <= as.Date("2022-12-31"))
cph_test  <- cph_df %>% filter(SALG_KVARTAL >  as.Date("2022-12-31"))

# Fjern SALG_KVARTAL før modellering
cph_train_model <- cph_train %>% select(-SALG_KVARTAL)
cph_cal_model   <- cph_cal   %>% select(-SALG_KVARTAL)
cph_test_model  <- cph_test  %>% select(-SALG_KVARTAL)

# ── Træn Random Forest ────────────────────────────────────────────────────────
cph_rf_model <- ranger(
  formula   = KONTANT_KOEBESUM ~ .,
  data      = cph_train_model,
  num.trees = 500,
  seed      = 123
)

# ── Conformal Prediction ──────────────────────────────────────────────────────
# Residualer på kalibreringsdata
cph_cal_preds <- predict(cph_rf_model, data = cph_cal_model)$predictions
cph_cal_resid <- abs(cph_cal_model$KONTANT_KOEBESUM - cph_cal_preds)

# Konfidensniveau
cph_alpha <- 0.10  # 90% dækningsgrad

# Konformal kvantil
cph_q_hat <- quantile(cph_cal_resid, probs = 1 - cph_alpha)

# Prediktionsintervaller på testdata
cph_test_preds <- predict(cph_rf_model, data = cph_test_model)$predictions

cph_intervals <- data.frame(
  faktisk = cph_test_model$KONTANT_KOEBESUM,
  pred    = cph_test_preds,
  lower   = cph_test_preds - cph_q_hat,
  upper   = cph_test_preds + cph_q_hat
)

# ── Evaluering ────────────────────────────────────────────────────────────────
cph_coverage <- mean(cph_intervals$faktisk >= cph_intervals$lower &
                       cph_intervals$faktisk <= cph_intervals$upper)
cat("Dækningsgrad:", round(cph_coverage * 100, 2), "%\n")

cph_avg_width <- mean(cph_intervals$upper - cph_intervals$lower)
cat("Gns. intervalbredde:", round(cph_avg_width, 0), "DKK\n")

# ── Visualisering ─────────────────────────────────────────────────────────────
ggplot(cph_intervals[1:200, ], aes(x = 1:200)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
  geom_point(aes(y = faktisk), color = "black", size = 0.8) +
  geom_line(aes(y = pred), color = "red") +
  labs(title = "Konformal prediktion - boligpriser (første 200 testobservationer)",
       x = "Observation",
       y = "Kontantpris (DKK)") +
  theme_minimal()

# Tjek om residualerne vokser over tid (tegn på ikke-stationaritet)
cph_test_intervals <- cph_intervals %>%
  mutate(dato = cph_test$SALG_KVARTAL,
         resid = abs(faktisk - pred))

ggplot(cph_test_intervals, aes(x = dato, y = resid)) +
  geom_point(alpha = 0.2, color = "steelblue") +
  geom_smooth(method = "loess", color = "red") +
  labs(title = "Residualer over tid (testperiode)",
       x = "Kvartal", y = "Absolut residual") +
  theme_minimal()






# ── Data forberedelse med log-transformation ──────────────────────────────────
cph_df <- CphHousingPrices %>%
  select(KONTANT_KOEBESUM, ENH_BEBO_ARL, ENH_VAERELSE_ANT,
         BYG_ALDER_AAR, kyst_AFSTAND, ENH_ANVEND_KODE,
         SALG_KVARTAL) %>%
  na.omit() %>%
  mutate(
    SALG_KVARTAL      = as.Date(SALG_KVARTAL),
    log_KOEBESUM      = log(KONTANT_KOEBESUM)  # log-transformer prisen
  )

# ── Tidsbaseret split ─────────────────────────────────────────────────────────
cph_train <- cph_df %>% filter(SALG_KVARTAL <= as.Date("2019-12-31"))
cph_cal   <- cph_df %>% filter(SALG_KVARTAL >  as.Date("2019-12-31") & 
                                 SALG_KVARTAL <= as.Date("2022-12-31"))
cph_test  <- cph_df %>% filter(SALG_KVARTAL >  as.Date("2022-12-31"))

cph_train_model <- cph_train %>% select(-SALG_KVARTAL, -KONTANT_KOEBESUM)
cph_cal_model   <- cph_cal   %>% select(-SALG_KVARTAL, -KONTANT_KOEBESUM)
cph_test_model  <- cph_test  %>% select(-SALG_KVARTAL, -KONTANT_KOEBESUM)

# ── Træn Random Forest på log-pris ────────────────────────────────────────────
cph_rf_model <- ranger(
  formula   = log_KOEBESUM ~ .,
  data      = cph_train_model,
  num.trees = 500,
  seed      = 123
)

# ── Conformal Prediction ──────────────────────────────────────────────────────
cph_cal_preds <- predict(cph_rf_model, data = cph_cal_model)$predictions
cph_cal_resid <- abs(cph_cal_model$log_KOEBESUM - cph_cal_preds)

cph_alpha <- 0.10
cph_q_hat <- quantile(cph_cal_resid, probs = 1 - cph_alpha)

cph_test_preds <- predict(cph_rf_model, data = cph_test_model)$predictions

# Transformer tilbage til DKK
cph_intervals <- data.frame(
  faktisk = exp(cph_test_model$log_KOEBESUM),
  pred    = exp(cph_test_preds),
  lower   = exp(cph_test_preds - cph_q_hat),
  upper   = exp(cph_test_preds + cph_q_hat)
)

# ── Evaluering ────────────────────────────────────────────────────────────────
cph_coverage <- mean(cph_intervals$faktisk >= cph_intervals$lower &
                       cph_intervals$faktisk <= cph_intervals$upper)
cat("Dækningsgrad:", round(cph_coverage * 100, 2), "%\n")

cph_avg_width <- mean(cph_intervals$upper - cph_intervals$lower)
cat("Gns. intervalbredde:", round(cph_avg_width, 0), "DKK\n")





# time based split 100 times



# ── Gentaget conformal prediction ─────────────────────────────────────────────
cph_n_rep <- 100
cph_results <- replicate(cph_n_rep, {
  
  # Tilfældig seed hver gang
  set.seed(sample(1:10000, 1))
  
  # Tidsbaseret split med lidt variation i grænserne
  cph_train_r <- cph_df %>% filter(SALG_KVARTAL <= as.Date("2019-12-31"))
  cph_cal_r   <- cph_df %>% filter(SALG_KVARTAL >  as.Date("2019-12-31") &
                                     SALG_KVARTAL <= as.Date("2022-12-31"))
  cph_test_r  <- cph_df %>% filter(SALG_KVARTAL >  as.Date("2022-12-31"))
  
  cph_train_r <- cph_train_r %>% select(-SALG_KVARTAL, -KONTANT_KOEBESUM)
  cph_cal_r   <- cph_cal_r   %>% select(-SALG_KVARTAL, -KONTANT_KOEBESUM)
  cph_test_r  <- cph_test_r  %>% select(-SALG_KVARTAL, -KONTANT_KOEBESUM)
  
  # Træn model
  cph_rf_r <- ranger(
    formula   = log_KOEBESUM ~ .,
    data      = cph_train_r,
    num.trees = 500,
    seed      = sample(1:10000, 1)
  )
  
  # Konformal kalibrering
  cph_cal_preds_r <- predict(cph_rf_r, data = cph_cal_r)$predictions
  cph_cal_resid_r <- abs(cph_cal_r$log_KOEBESUM - cph_cal_preds_r)
  cph_q_hat_r <- quantile(cph_cal_resid_r, 
                          probs = ceiling((length(cph_cal_resid_r) + 1) * 
                                            (1 - cph_alpha)) / length(cph_cal_resid_r))
  
  # Test
  cph_test_preds_r <- predict(cph_rf_r, data = cph_test_r)$predictions
  
  cph_intervals_r <- data.frame(
    faktisk = exp(cph_test_r$log_KOEBESUM),
    lower   = exp(cph_test_preds_r - cph_q_hat_r),
    upper   = exp(cph_test_preds_r + cph_q_hat_r)
  )
  
  # Returnér dækningsgrad og intervalbredde
  c(
    coverage = mean(cph_intervals_r$faktisk >= cph_intervals_r$lower &
                      cph_intervals_r$faktisk <= cph_intervals_r$upper),
    avg_width = mean(cph_intervals_r$upper - cph_intervals_r$lower)
  )
})

# ── Resultater ────────────────────────────────────────────────────────────────
cph_coverage_vec <- cph_results["coverage", ]
cph_width_vec    <- cph_results["avg_width", ]

cat("Gns. dækningsgrad:", round(mean(cph_coverage_vec) * 100, 2), "%\n")
cat("Std. dækningsgrad:", round(sd(cph_coverage_vec) * 100, 2), "%\n")
cat("Gns. intervalbredde:", round(mean(cph_width_vec), 0), "DKK\n")


# Visualisering af dækningsgrader på tværs af gentagelser
hist(cph_coverage_vec, 
     main = "Dækningsgrad over 100 gentagelser",
     xlab = "Dækningsgrad",
     col  = "steelblue")
abline(v = 0.90, col = "red", lwd = 2, lty = 2)  # 90% målstreg






# random split 100 times


cph_n_rep <- 100
cph_results_random <- replicate(cph_n_rep, {
  
  set.seed(sample(1:10000, 1))
  
  cph_n <- nrow(cph_df)
  cph_train_idx <- sample(1:cph_n, size = 0.6 * cph_n)
  cph_cal_idx   <- sample(setdiff(1:cph_n, cph_train_idx), size = 0.2 * cph_n)
  cph_test_idx  <- setdiff(1:cph_n, c(cph_train_idx, cph_cal_idx))
  
  cph_train_r <- cph_df[cph_train_idx, ] %>% select(-SALG_KVARTAL, -KONTANT_KOEBESUM)
  cph_cal_r   <- cph_df[cph_cal_idx, ]   %>% select(-SALG_KVARTAL, -KONTANT_KOEBESUM)
  cph_test_r  <- cph_df[cph_test_idx, ]  %>% select(-SALG_KVARTAL, -KONTANT_KOEBESUM)
  
  cph_rf_r <- ranger(
    formula   = log_KOEBESUM ~ .,
    data      = cph_train_r,
    num.trees = 500,
    seed      = sample(1:10000, 1)
  )
  
  cph_cal_preds_r <- predict(cph_rf_r, data = cph_cal_r)$predictions
  cph_cal_resid_r <- abs(cph_cal_r$log_KOEBESUM - cph_cal_preds_r)
  cph_q_hat_r <- quantile(cph_cal_resid_r, 
                          probs = ceiling((length(cph_cal_resid_r) + 1) * 
                                            (1 - cph_alpha)) / length(cph_cal_resid_r))
  
  cph_test_preds_r <- predict(cph_rf_r, data = cph_test_r)$predictions
  
  cph_intervals_r <- data.frame(
    faktisk = exp(cph_test_r$log_KOEBESUM),
    lower   = exp(cph_test_preds_r - cph_q_hat_r),
    upper   = exp(cph_test_preds_r + cph_q_hat_r)
  )
  
  c(
    coverage  = mean(cph_intervals_r$faktisk >= cph_intervals_r$lower &
                       cph_intervals_r$faktisk <= cph_intervals_r$upper),
    avg_width = mean(cph_intervals_r$upper - cph_intervals_r$lower)
  )
})

cph_coverage_vec_random <- cph_results_random["coverage", ]
cph_width_vec_random    <- cph_results_random["avg_width", ]

cat("Gns. dækningsgrad:", round(mean(cph_coverage_vec_random) * 100, 2), "%\n")
cat("Std. dækningsgrad:", round(sd(cph_coverage_vec_random) * 100, 2), "%\n")
cat("Gns. intervalbredde:", round(mean(cph_width_vec_random), 0), "DKK\n")

hist(cph_coverage_vec_random,
     main = "Dækningsgrad over 100 gentagelser (random split)",
     xlab = "Dækningsgrad",
     col  = "steelblue", breaks = 15)
abline(v = 0.90, col = "red", lwd = 2, lty = 2)




# =========================
# 1. Data preprocessing
# =========================
cph_data <- CphHousingPrices

cph_data <- cph_data %>%
  mutate(
    cph_log_price = log(KONTANT_KOEBESUM),
    BYG_ALDER_AAR = pmax(0, BYG_ALDER_AAR),
    KOMMUNE_NAVN = as.factor(KOMMUNE_NAVN),
    LANDSDEL_NAVN = as.factor(LANDSDEL_NAVN)
  ) %>%
  select(-KONTANT_KOEBESUM, -SALG_KVARTAL)

# imputering
cph_data <- missRanger(cph_data, verbose = FALSE)


# =========================
# 2. Split: train / cal / test
# =========================
set.seed(1)

cph_n <- nrow(cph_data)
cph_idx <- sample(1:cph_n)

cph_train_idx <- cph_idx[1:floor(0.6*cph_n)]
cph_cal_idx   <- cph_idx[(floor(0.6*cph_n)+1):floor(0.8*cph_n)]
cph_test_idx  <- cph_idx[(floor(0.8*cph_n)+1):cph_n]

cph_train <- cph_data[cph_train_idx, ]
cph_cal   <- cph_data[cph_cal_idx, ]
cph_test  <- cph_data[cph_test_idx, ]

# =========================
# 3. Random Forest (ranger)
# =========================
cph_rf_model <- ranger(
  cph_log_price ~ .,
  data = cph_train,
  num.trees = 500,
  mtry = floor(sqrt(ncol(cph_train) - 1)),
  min.node.size = 5,
  seed = 1
)

# =========================
# 4. Split Conformal Prediction
# =========================

# predictions (calibration)
cph_pred_cal <- predict(cph_rf_model, data = cph_cal)$predictions

# nonconformity scores
cph_scores <- abs(cph_cal$cph_log_price - cph_pred_cal)
head(cph_scores)
# kvantil (finite-sample korrekt)
cph_alpha <- 0.1
cph_n_cal <- base::length(cph_scores)

cph_q_hat <- quantile(
  cph_scores,
  probs = ceiling((cph_n_cal + 1) * (1 - cph_alpha)) / cph_n_cal
)

# prediction (test)
cph_pred_test <- predict(cph_rf_model, data = cph_test)$predictions

cph_lower <- cph_pred_test - cph_q_hat
cph_upper <- cph_pred_test + cph_q_hat

# =========================
# 5. Evaluering (standard CP)
# =========================
cph_coverage <- mean(
  cph_test$cph_log_price >= cph_lower &
    cph_test$cph_log_price <= cph_upper
)

cph_mean_width <- mean(cph_upper - cph_lower)

# =========================
# 6. Adaptive / Normalized CP
# =========================

# residual model
cph_pred_train <- predict(cph_rf_model, data = cph_train)$predictions
cph_res_train <- abs(cph_train$cph_log_price - cph_pred_train)

cph_train_sigma <- cph_train
cph_train_sigma$cph_residual <- cph_res_train

cph_rf_sigma <- ranger(
  cph_residual ~ .,
  data = cph_train_sigma,
  num.trees = 500,
  seed = 1
)

# normaliserede scores
cph_sigma_cal <- predict(cph_rf_sigma, data = cph_cal)$predictions
cph_scores_norm <- abs(cph_cal$cph_log_price - cph_pred_cal) / cph_sigma_cal

cph_q_hat_norm <- quantile(
  cph_scores_norm,
  probs = ceiling((cph_n_cal + 1) * (1 - cph_alpha)) / cph_n_cal
)

# adaptive intervals
cph_sigma_test <- predict(cph_rf_sigma, data = cph_test)$predictions

cph_lower_adapt <- cph_pred_test - cph_q_hat_norm * cph_sigma_test
cph_upper_adapt <- cph_pred_test + cph_q_hat_norm * cph_sigma_test

# =========================
# 7. Evaluering (adaptive CP)
# =========================
cph_coverage_adapt <- mean(
  cph_test$cph_log_price >= cph_lower_adapt &
    cph_test$cph_log_price <= cph_upper_adapt
)

cph_width_adapt <- mean(cph_upper_adapt - cph_lower_adapt)

# =========================
# 8. Resultater
# =========================
cph_results <- data.frame(
  method = c("Split CP", "Adaptive CP"),
  coverage = c(cph_coverage, cph_coverage_adapt),
  width = c(cph_mean_width, cph_width_adapt)
)

print(cph_results)




# FLERE RUNS




# =========================
# 0. Pakker
# =========================
library(ranger)
library(dplyr)
library(ggplot2)

# =========================
# 1. Preprocessing
# =========================
cph_prep_data <- function(data) {
  data <- data %>%
    mutate(
      cph_log_price = log(KONTANT_KOEBESUM),
      BYG_ALDER_AAR = pmax(0, BYG_ALDER_AAR),
      KOMMUNE_NAVN = as.factor(KOMMUNE_NAVN),
      LANDSDEL_NAVN = as.factor(LANDSDEL_NAVN)
    ) %>%
    select(-KONTANT_KOEBESUM, -SALG_KVARTAL)
  
  # simpel imputering
  cph_num_vars <- sapply(data, is.numeric)
  data[cph_num_vars] <- lapply(data[cph_num_vars], function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })
  
  cph_cat_vars <- sapply(data, is.character)
  data[cph_cat_vars] <- lapply(data[cph_cat_vars], function(x) {
    x[is.na(x)] <- "missing"
    as.factor(x)
  })
  
  return(data)
}

# =========================
# 2. Ét run
# =========================
cph_run_once <- function(data, alpha = 0.1) {
  
  cph_n <- nrow(data)
  cph_idx <- sample(1:cph_n)
  
  cph_train <- data[cph_idx[1:floor(0.6*cph_n)], ]
  cph_cal   <- data[cph_idx[(floor(0.6*cph_n)+1):floor(0.8*cph_n)], ]
  cph_test  <- data[cph_idx[(floor(0.8*cph_n)+1):cph_n], ]
  
  # model
  cph_rf <- ranger(cph_log_price ~ ., data = cph_train, num.trees = 300)
  
  # predictions
  cph_pred_train <- predict(cph_rf, cph_train)$predictions
  cph_pred_cal   <- predict(cph_rf, cph_cal)$predictions
  cph_pred_test  <- predict(cph_rf, cph_test)$predictions
  
  # =====================
  # Standard CP
  # =====================
  cph_scores <- abs(cph_cal$cph_log_price - cph_pred_cal)
  cph_n_cal <- length(cph_scores)
  
  cph_q_hat <- as.numeric(quantile(
    cph_scores,
    probs = ceiling((cph_n_cal + 1) * (1 - alpha)) / cph_n_cal
  ))
  
  cph_lower <- cph_pred_test - cph_q_hat
  cph_upper <- cph_pred_test + cph_q_hat
  
  cph_coverage <- mean(cph_test$cph_log_price >= cph_lower &
                         cph_test$cph_log_price <= cph_upper)
  
  cph_width <- mean(cph_upper - cph_lower)
  
  # =====================
  # Adaptive CP
  # =====================
  cph_res_train <- abs(cph_train$cph_log_price - cph_pred_train)
  cph_train$cph_residual <- cph_res_train
  
  cph_rf_sigma <- ranger(cph_residual ~ ., data = cph_train, num.trees = 300)
  
  cph_sigma_cal  <- predict(cph_rf_sigma, cph_cal)$predictions
  cph_sigma_test <- predict(cph_rf_sigma, cph_test)$predictions
  
  # stabilitet
  cph_sigma_cal[cph_sigma_cal <= 1e-8] <- 1e-8
  cph_sigma_test[cph_sigma_test <= 1e-8] <- 1e-8
  
  cph_scores_norm <- abs(cph_cal$cph_log_price - cph_pred_cal) / cph_sigma_cal
  
  cph_q_hat_norm <- as.numeric(quantile(
    cph_scores_norm,
    probs = ceiling((cph_n_cal + 1) * (1 - alpha)) / cph_n_cal
  ))
  
  cph_lower_adapt <- cph_pred_test - cph_q_hat_norm * cph_sigma_test
  cph_upper_adapt <- cph_pred_test + cph_q_hat_norm * cph_sigma_test
  
  cph_coverage_adapt <- mean(cph_test$cph_log_price >= cph_lower_adapt &
                               cph_test$cph_log_price <= cph_upper_adapt)
  
  cph_width_adapt <- mean(cph_upper_adapt - cph_lower_adapt)
  
  return(data.frame(
    cph_cp_coverage = cph_coverage,
    cph_cp_width = cph_width,
    cph_adapt_coverage = cph_coverage_adapt,
    cph_adapt_width = cph_width_adapt
  ))
}

# =========================
# 3. Gentagelser
# =========================
set.seed(1)

cph_data <- cph_prep_data(CphHousingPrices)

cph_B <- 30

cph_results_runs <- bind_rows(
  lapply(1:cph_B, function(i) cph_run_once(cph_data))
)

# =========================
# 4. Summary
# =========================
cph_summary <- data.frame(
  method = c("Split CP", "Adaptive CP"),
  coverage_mean = c(mean(cph_results_runs$cph_cp_coverage),
                    mean(cph_results_runs$cph_adapt_coverage)),
  coverage_sd   = c(sd(cph_results_runs$cph_cp_coverage),
                    sd(cph_results_runs$cph_adapt_coverage)),
  width_mean    = c(mean(cph_results_runs$cph_cp_width),
                    mean(cph_results_runs$cph_adapt_width))
)

print(cph_summary)

# =========================
# 5. Plots
# =========================

# Coverage
cph_df_cov <- data.frame(
  value = c(cph_results_runs$cph_cp_coverage,
            cph_results_runs$cph_adapt_coverage),
  method = rep(c("Split CP", "Adaptive CP"), each = cph_B)
)

ggplot(cph_df_cov, aes(x = value, fill = method)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0.9, linetype = "dashed") +
  labs(title = "CP Coverage (Copenhagen Housing)",
       x = "Coverage")

# Width
cph_df_width <- data.frame(
  value = c(cph_results_runs$cph_cp_width,
            cph_results_runs$cph_adapt_width),
  method = rep(c("Split CP", "Adaptive CP"), each = cph_B)
)

ggplot(cph_df_width, aes(x = method, y = value, fill = method)) +
  geom_boxplot(alpha = 0.6) +
  labs(title = "Interval Width (Copenhagen Housing)",
       y = "Width")

# Coverage over runs
cph_df_runs <- data.frame(
  run = 1:cph_B,
  cp = cph_results_runs$cph_cp_coverage,
  adaptive = cph_results_runs$cph_adapt_coverage
)

ggplot(cph_df_runs, aes(x = run)) +
  geom_line(aes(y = cp, color = "Split CP")) +
  geom_line(aes(y = adaptive, color = "Adaptive CP")) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  labs(title = "Coverage Stability (Copenhagen Housing)",
       y = "Coverage")


