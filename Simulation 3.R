# =============================
# Setup
# =============================
set.seed(1)

library(ggplot2)
library(dplyr)
library(ranger)

# Parameters
n_list <- c(20, 80, 320, 1280)
M      <- 200
n_test <- 100000
alpha  <- 0.1
sigma  <- 1

f <- function(x) sin(2 * pi * x)

# Storage
results_all <- data.frame()

# =============================
# Simulation
# =============================
for (n in n_list) {
  cat("Running n =", n, "\n")
  n_train <- n / 2
  n_cal   <- n / 2
  
  for (m in 1:M) {
    gen_data <- function(n) {
      x <- runif(n, -1, 1)
      data.frame(x = x, y = f(x) + rnorm(n, 0, sigma))
    }
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # Fit ranger
    rf_fit     <- ranger(y ~ x, data = train)
    pred_train <- predict(rf_fit, train)$predictions
    pred_cal   <- predict(rf_fit, cal)$predictions
    pred_test  <- predict(rf_fit, test)$predictions
    
    # =========================
    # (1) Oracle
    # =========================
    z            <- qnorm(1 - alpha / 2)
    oracle_lower <- f(test$x) - z * sigma
    oracle_upper <- f(test$x) + z * sigma
    
    # =========================
    # (2) Gaussian
    # =========================
    sigma_hat   <- sqrt(mean((train$y - pred_train)^2))
    gauss_lower <- pred_test - qnorm(1 - alpha / 2) * sigma_hat
    gauss_upper <- pred_test + qnorm(1 - alpha / 2) * sigma_hat
    
    # =========================
    # (3) Residual quantile (RQ)
    # =========================
    res_cal  <- cal$y - pred_cal
    rq_lower <- pred_test + quantile(res_cal, alpha / 2)
    rq_upper <- pred_test + quantile(res_cal, 1 - alpha / 2)
    
    # =========================
    # (4) Conformal (CP)
    # =========================
    k        <- ceiling((n_cal + 1) * (1 - alpha))
    q_cp     <- sort(abs(cal$y - pred_cal))[k]
    cp_lower <- pred_test - q_cp
    cp_upper <- pred_test + q_cp
    
    # =========================
    # (5) Quantile regression (QR)
    # =========================
    qr_fit   <- ranger(y ~ x, data = cal, quantreg = TRUE)
    qr_preds <- predict(qr_fit, test, type = "quantiles",
                        quantiles = c(alpha / 2, 1 - alpha / 2))$predictions
    qr_lower <- qr_preds[, 1]
    qr_upper <- qr_preds[, 2]
    
    # Metrics
    cov_fn <- function(l, u) mean(test$y >= l & test$y <= u)
    len_fn <- function(l, u) mean(u - l)
    
    results_all <- rbind(results_all, data.frame(
      n      = n,
      m      = m,
      method = c("Oracle", "Gaussian", "RQ", "CP", "QR"),
      coverage = c(cov_fn(oracle_lower, oracle_upper),
                   cov_fn(gauss_lower,  gauss_upper),
                   cov_fn(rq_lower,     rq_upper),
                   cov_fn(cp_lower,     cp_upper),
                   cov_fn(qr_lower,     qr_upper)),
      length = c(len_fn(oracle_lower, oracle_upper),
                 len_fn(gauss_lower,  gauss_upper),
                 len_fn(rq_lower,     rq_upper),
                 len_fn(cp_lower,     cp_upper),
                 len_fn(qr_lower,     qr_upper))
    ))
  }
}

# =============================
# Post-processing
# =============================
results_all <- results_all |>
  mutate(n = as.numeric(n))

oracle_length <- results_all |>
  filter(method == "Oracle") |>
  group_by(n) |>
  summarise(L_oracle = mean(length))

results_all <- left_join(results_all, oracle_length, by = "n")

# Summary for længdeplot
summary_diff <- results_all |>
  filter(method != "Oracle") |>
  mutate(length_diff = sqrt(n) * (length - L_oracle)) |>
  group_by(n, method) |>
  summarise(
    mean_len = mean(length_diff),
    se_len   = sd(length_diff) / sqrt(M),
    .groups  = "drop"
  )

# Summary for coverageplot
mean_cov <- results_all |>
  filter(method != "Oracle") |>
  group_by(n, method) |>
  summarise(mean_cov = mean(coverage), .groups = "drop")

# =============================
# Plot: Coverage
# =============================
ggplot(results_all |> filter(method != "Oracle"),
       aes(x = n, y = coverage, color = method)) +
  geom_point(alpha = 0.3, size = 0.8,
             position = position_jitter(width = 0.03)) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "linerange", linewidth = 0.3, alpha = 0.5) +
  geom_hline(yintercept = 1 - alpha, color = "black") +
  geom_text(data = mean_cov,
            aes(x = n, y = 1.02,
                label = scales::percent(mean_cov, accuracy = 1)),
            size = 3) +
  coord_cartesian(ylim = c(NA, 1.05)) +
  scale_x_log10(breaks = n_list, labels = n_list) +
  facet_wrap(~ method) +
  labs(y = "Coverage", x = "Sample size n") +
  theme_bw() +
  theme(legend.position = "none")

# =============================
# Plot: Length (scaled)
# =============================
ggplot(summary_diff, aes(x = n, y = mean_len, color = method, group = method)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = mean_len - 2 * se_len,
                    ymax = mean_len + 2 * se_len), width = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_log10(breaks = n_list, labels = n_list) +
  labs(y = expression(sqrt(n) * (L[n] - L[oracle])), x = "Sample size n") +
  theme_bw()
