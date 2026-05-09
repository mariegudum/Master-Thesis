# ── Pakker ────────────────────────────────────────────────────────────────────
library(ranger)
library(dplyr)
library(ggplot2)

# ── Parametre ─────────────────────────────────────────────────────────────────
cph_alpha <- 0.10
cph_n_rep <- 100

# ── Data forberedelse ─────────────────────────────────────────────────────────
cph_df <- CphHousingPrices %>%
  na.omit() %>%
  mutate(
    log_KOEBESUM   = log(KONTANT_KOEBESUM),
    cph_salgsaar   = as.integer(format(SALG_KVARTAL, "%Y")),
    ENH_BAD_KODE   = as.factor(ENH_BAD_KODE),
    LANDSDEL_NAVN  = as.factor(LANDSDEL_NAVN),
    KOMMUNE_NAVN   = as.factor(KOMMUNE_NAVN)
  ) %>%
  select(-KONTANT_KOEBESUM, -SALG_KVARTAL)

# ── Split funktion ────────────────────────────────────────────────────────────
cph_make_split <- function(df) {
  df %>%
    group_by(cph_salgsaar) %>%
    mutate(
      cph_n   = n(),
      cph_idx = sample.int(cph_n[1]),
      split   = case_when(
        cph_idx <= round(0.6 * cph_n) ~ "train",
        cph_idx <= round(0.8 * cph_n) ~ "cal",
        TRUE                          ~ "test"
      )
    ) %>%
    ungroup()
}

# ── Gentaget conformal prediction ─────────────────────────────────────────────
cph_coverage_vec_strat <- numeric(cph_n_rep)
cph_width_vec_strat    <- numeric(cph_n_rep)

for (i in 1:cph_n_rep) {
  
  cph_split <- cph_make_split(cph_df)
  
  cph_train_r <- cph_split %>% filter(split == "train") %>%
    select(-split, -cph_idx, -cph_n)
  cph_cal_r   <- cph_split %>% filter(split == "cal")   %>%
    select(-split, -cph_idx, -cph_n)
  cph_test_r  <- cph_split %>% filter(split == "test")  %>%
    select(-split, -cph_idx, -cph_n)
  
  cph_rf_r <- ranger(
    formula     = log_KOEBESUM ~ .,
    data        = cph_train_r,
    num.trees   = 500,
    num.threads = parallel::detectCores() - 1,
    seed        = sample(1:10000, 1)
  )
  
  cph_cal_preds_r <- predict(cph_rf_r, data = cph_cal_r)$predictions
  cph_cal_resid_r <- abs(cph_cal_r$log_KOEBESUM - cph_cal_preds_r)
  cph_q_hat_r     <- quantile(cph_cal_resid_r,
                              probs = ceiling((base::length(cph_cal_resid_r) + 1) *
                                                (1 - cph_alpha)) / base::length(cph_cal_resid_r))
  
  cph_test_preds_r <- predict(cph_rf_r, data = cph_test_r)$predictions
  
  cph_intervals_r <- data.frame(
    faktisk = exp(cph_test_r$log_KOEBESUM),
    lower   = exp(cph_test_preds_r - cph_q_hat_r),
    upper   = exp(cph_test_preds_r + cph_q_hat_r)
  )
  
  cph_coverage_vec_strat[i] <- mean(cph_intervals_r$faktisk >= cph_intervals_r$lower &
                                      cph_intervals_r$faktisk <= cph_intervals_r$upper)
  cph_width_vec_strat[i]    <- mean(cph_intervals_r$upper - cph_intervals_r$lower)
  
  cat("Iteration", i, "færdig - dækningsgrad:",
      round(cph_coverage_vec_strat[i] * 100, 2), "%\n")
}

# ── Resultater ────────────────────────────────────────────────────────────────
cat("Gns. dækningsgrad:", round(mean(cph_coverage_vec_strat) * 100, 2), "%\n")
cat("Std. dækningsgrad:", round(sd(cph_coverage_vec_strat) * 100, 2), "%\n")
cat("Gns. intervalbredde:", round(mean(cph_width_vec_strat), 0), "DKK\n")

# ── Visualisering ─────────────────────────────────────────────────────────────
hist(cph_coverage_vec_strat,
     main = "Coverage over 100 replications (stratifed split)", breaks = 27,
     xlab = "Coverage",
     col  = "steelblue")
abline(v = 0.90, col = "red", lwd = 2, lty = 2)
abline(v = mean(cph_coverage_vec_strat), lwd = 2)
abline(v = mean(cph_coverage_vec_strat)-sd(cph_coverage_vec_strat), lwd = 2, lty = 2)
abline(v = mean(cph_coverage_vec_strat)+sd(cph_coverage_vec_strat), lwd = 2, lty = 2)


hist(cph_width_vec_strat,
     main = "Interval width over 100 replications (stratified split)",
     breaks = 27,
     xlab = "Interval width (DKK)",
     col  = "steelblue")

abline(v = mean(cph_width_vec_strat), lwd = 2)
abline(v = mean(cph_width_vec_strat) - sd(cph_width_vec_strat),
       lwd = 2, lty = 2)
abline(v = mean(cph_width_vec_strat) + sd(cph_width_vec_strat),
       lwd = 2, lty = 2)


plot(cph_width_vec_strat,
     cph_coverage_vec_strat,
     xlab = "Width",
     ylab = "Coverage",
     main = "Coverage vs Width")

abline(h = 0.9, col = "red", lty = 2)

plot(cph_intervals_r$faktisk, 
     cph_intervals_r$upper - cph_intervals_r$lower)


cph_intervals_r$bin <- cut(cph_intervals_r$faktisk, 10)

aggregate((faktisk >= lower & faktisk <= upper) ~ bin,
          data = cph_intervals_r,
          mean)

cph_intervals_r$year <- cph_test_r$cph_salgsaar

aggregate((faktisk >= lower & faktisk <= upper) ~ year,
          data = cph_intervals_r,
          mean)

hist(cph_width_vec_strat)


rel_width <- (cph_intervals_r$upper - cph_intervals_r$lower) /
  cph_intervals_r$faktisk

plot(cph_intervals_r$faktisk, rel_width)



plot(log(cph_intervals_r$faktisk),
     log(cph_intervals_r$upper - cph_intervals_r$lower))




cph_intervals_r$residual <- abs(cph_test_r$log_KOEBESUM - cph_test_preds_r)

plot(cph_intervals_r$residual,
     cph_intervals_r$upper - cph_intervals_r$lower,
     xlab = "Model error",
     ylab = "Interval width")


# =========================
# 1. Lav coverage indikator + bins
# =========================
cph_intervals_r <- cph_intervals_r %>%
  mutate(
    covered = (faktisk >= lower & faktisk <= upper),
    price_bin = cut(faktisk, breaks = 10)
  )

# =========================
# 2. Beregn coverage per bin
# =========================
cph_cov_bin <- cph_intervals_r %>%
  group_by(price_bin) %>%
  summarise(
    coverage = mean(covered),
    count = n(),
    .groups = "drop"
  )

print(cph_cov_bin)

# =========================
# 3. Plot (ggplot - anbefalet)
# =========================
ggplot(cph_cov_bin, aes(x = price_bin, y = coverage, group = 1)) +
  geom_point(size = 2) +
  geom_line() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  labs(
    title = "Conditional Coverage by Price Bin",
    x = "Price bin",
    y = "Coverage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(cph_cov_bin, aes(x = price_bin, y = count)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Observations per Price Bin",
    x = "Price bin",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# lav numeriske breaks
breaks <- quantile(cph_intervals_r$faktisk, probs = seq(0, 1, length.out = 11))

# lav labels i mio DKK
labels <- paste0(
  round(breaks[-base::length(breaks)] / 1e6, 1), "–",
  round(breaks[-1] / 1e6, 1), " mio"
)

# brug cut med labels
cph_intervals_r$price_bin <- cut(
  cph_intervals_r$faktisk,
  breaks = breaks,
  labels = labels,
  include.lowest = TRUE
)

cph_cov_bin <- cph_intervals_r %>%
  mutate(covered = (faktisk >= lower & faktisk <= upper)) %>%
  group_by(price_bin) %>%
  summarise(coverage = mean(covered), .groups = "drop")

ggplot(cph_cov_bin, aes(x = price_bin, y = coverage, group = 1)) +
  geom_point(size = 2) +
  geom_line() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  labs(
    title = "Conditional Coverage by Price Bin",
    x = "Price (million DKK)",
    y = "Coverage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# =========================
# 1. Intervalbredder
# =========================

# log-skala (konstant!)
cph_width_log <- (cph_test_preds_r + cph_q_hat_r) - 
  (cph_test_preds_r - cph_q_hat_r)

# original skala (ikke konstant)
cph_width_orig <- cph_intervals_r$upper - cph_intervals_r$lower

# relativ bredde
cph_rel_width <- cph_width_orig / cph_intervals_r$faktisk


# =========================
# 2. Tjek variation
# =========================

cat("Std log-width:", sd(cph_width_log), "\n")
cat("Std original width:", sd(cph_width_orig), "\n")
cat("Std relative width:", sd(cph_rel_width), "\n")


# =========================
# 3. Plots
# =========================

par(mfrow = c(1,3))

# (A) Log-scale width
plot(cph_intervals_r$faktisk, cph_width_log,
     main = "Width (log-scale)",
     xlab = "Price",
     ylab = "Width",
     col = "steelblue")

# (B) Original width
plot(cph_intervals_r$faktisk, cph_width_orig,
     main = "Width (original scale)",
     xlab = "Price",
     ylab = "Width",
     col = "darkgreen")

# (C) Relative width
plot(cph_intervals_r$faktisk, cph_rel_width,
     main = "Relative width",
     xlab = "Price",
     ylab = "Width / Price",
     col = "firebrick")




cph_analysis_df <- cph_intervals_r %>%
  mutate(
    covered = (faktisk >= lower & faktisk <= upper),
    width   = upper - lower,
    residual = abs(cph_test_r$log_KOEBESUM - cph_test_preds_r),
    
    # Quantile bins (pæne og stabile)
    price_bin = cut(
      faktisk,
      breaks = quantile(faktisk, probs = seq(0,1,length.out = 11)),
      include.lowest = TRUE
    )
  )


cph_analysis_df <- cph_analysis_df %>%
  mutate(
    size_bin = cut(
      ENH_BEBO_ARL,
      breaks = quantile(ENH_BEBO_ARL, probs = seq(0,1,length.out = 11)),
      include.lowest = TRUE
    )
  )

cph_cov_size <- cph_analysis_df %>%
  group_by(size_bin) %>%
  summarise(
    coverage = mean(covered),
    .groups = "drop"
  )


ggplot(cph_cov_size, aes(x = size_bin, y = coverage, group = 1)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Coverage vs house size")





library(dplyr)
library(ggplot2)

# =========================
# 1. Saml data
# =========================
cph_analysis_df <- cph_intervals_r %>%
  mutate(
    covered = (faktisk >= lower & faktisk <= upper),
    width   = upper - lower
  ) %>%
  bind_cols(cph_test_r)

# =========================
# 2. Vælg numeriske covariater
# =========================
num_vars <- names(cph_analysis_df)[sapply(cph_analysis_df, is.numeric)]

num_vars <- setdiff(num_vars, c(
  "faktisk", "lower", "upper",
  "width", "covered", "log_KOEBESUM"
))

num_vars <- names(cph_analysis_df)[sapply(cph_analysis_df, is.numeric)]

num_vars <- setdiff(num_vars, c(
  "faktisk", "lower", "upper",
  "width", "covered",
  "residual",     
  "year"          
))

# =========================
# 3. Robust bin-funktion
# =========================
make_bins <- function(x, n_bins = 5) {
  
  # fjern NA
  x <- x[!is.na(x)]
  
  # hvis for få unikke værdier → skip
  if(length(unique(x)) < 10) return(NULL)
  
  breaks <- unique(quantile(x, probs = seq(0,1,length.out = n_bins+1), na.rm = TRUE))
  
  # hvis stadig for få bins → skip
  if(length(breaks) < 3) return(NULL)
  
  return(breaks)
}

# =========================
# 4. Coverage per variabel
# =========================
compute_cov_bins <- function(var_name) {
  
  x <- cph_analysis_df[[var_name]]
  breaks <- make_bins(x)
  
  if(is.null(breaks)) return(NULL)
  
  bins <- cut(x, breaks = breaks, include.lowest = TRUE)
  
  df <- data.frame(
    bin = bins,
    covered = cph_analysis_df$covered
  )
  
  out <- df %>%
    group_by(bin) %>%
    summarise(
      coverage = mean(covered),
      n = n(),
      .groups = "drop"
    )
  
  out$variable <- var_name
  return(out)
}

# =========================
# 5. Width per variabel
# =========================
compute_width_bins <- function(var_name) {
  
  x <- cph_analysis_df[[var_name]]
  breaks <- make_bins(x)
  
  if(is.null(breaks)) return(NULL)
  
  bins <- cut(x, breaks = breaks, include.lowest = TRUE)
  
  df <- data.frame(
    bin = bins,
    width = cph_analysis_df$width
  )
  
  out <- df %>%
    group_by(bin) %>%
    summarise(
      mean_width = mean(width),
      .groups = "drop"
    )
  
  out$variable <- var_name
  return(out)
}

# =========================
# 6. Loop over variabler
# =========================
cov_list <- lapply(num_vars, compute_cov_bins)
cov_list <- cov_list[!sapply(cov_list, is.null)]
cov_all  <- bind_rows(cov_list)

width_list <- lapply(num_vars, compute_width_bins)
width_list <- width_list[!sapply(width_list, is.null)]
width_all  <- bind_rows(width_list)

# =========================
# 7. Coverage plot
# =========================
p_cov <- ggplot(cov_all, aes(x = bin, y = coverage, group = 1)) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
  facet_wrap(~variable, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8)
  ) +
  labs(
    title = "Conditional coverage across covariates",
    x = "Bins",
    y = "Coverage"
  )

# =========================
# 8. Width plot
# =========================
p_width <- ggplot(width_all, aes(x = bin, y = mean_width, group = 1)) +
  geom_point(size = 1.5, color = "darkgreen") +
  geom_line(color = "darkgreen") +
  facet_wrap(~variable, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8)
  ) +
  labs(
    title = "Conditional interval width across covariates",
    x = "Bins",
    y = "Width"
  )

# =========================
# 9. Print plots
# =========================
print(p_cov)
print(p_width)

# Ideer fra Anton

ggplot(bin_results, aes(x = bin_mid, y = mean_width, color = coverage)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    y = "Interval width",
    color = "Coverage"
  )

cph_df %>%
  group_by(cph_salgsaar) %>%
  summarise(
    mean_width = mean(width),
    coverage = mean(covered)
  )

ggplot(df, aes(x = cph_salgsaar, y = mean_width, color = coverage)) +
  geom_line() +
  geom_point()

df_analysis <- cph_intervals %>%
  mutate(
    residual = abs(y_true - y_pred),
    width = upper - lower
  )

cor(df_analysis$width, df_analysis$residual)

ggplot(df_analysis, aes(x = residual, y = width)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm")

# residuals

plot(cph_cal_preds_r, cph_cal_resid_r,
     pch = 16, cex = 0.5,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, col = "red")
