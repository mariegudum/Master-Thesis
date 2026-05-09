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
    rf_fit     <- ranger(y ~ x, data = train, num.trees = 500)
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
    
    # Metrics
    cov_fn <- function(l, u) mean(test$y >= l & test$y <= u)
    len_fn <- function(l, u) mean(u - l)
    
    results_all <- rbind(results_all, data.frame(
      n      = n,
      m      = m,
      method = c("Oracle", "Gaussian", "RQ", "CP"),
      coverage = c(cov_fn(oracle_lower, oracle_upper),
                   cov_fn(gauss_lower,  gauss_upper),
                   cov_fn(rq_lower,     rq_upper),
                   cov_fn(cp_lower,     cp_upper)),
      length = c(len_fn(oracle_lower, oracle_upper),
                 len_fn(gauss_lower,  gauss_upper),
                 len_fn(rq_lower,     rq_upper),
                 len_fn(cp_lower,     cp_upper))
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

summary_diff <- results_all |>
  filter(method != "Oracle") |>
  mutate(length_diff = sqrt(n) * (length - L_oracle)) |>
  group_by(n, method) |>
  summarise(
    mean_len = mean(length_diff),
    se_len   = sd(length_diff) / sqrt(M),
    .groups  = "drop"
  )

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





# =============================
# Plot: Bånd på data
# =============================
set.seed(42)
n_show  <- 1000
x_show  <- runif(n_show, -1, 1)
y_show  <- f(x_show) + rnorm(n_show, 0, sigma)

# Træn på separat data
n_band  <- 80
n_train_band <- n_band / 2
n_cal_band   <- n_band / 2

train_band <- data.frame(x = runif(n_train_band, -1, 1)) |>
  mutate(y = f(x) + rnorm(n_train_band, 0, sigma))
cal_band   <- data.frame(x = runif(n_cal_band, -1, 1)) |>
  mutate(y = f(x) + rnorm(n_cal_band, 0, sigma))

rf_band     <- ranger(y ~ x, data = train_band, num.trees = 500)
pred_train_band <- predict(rf_band, train_band)$predictions
pred_cal_band   <- predict(rf_band, cal_band)$predictions

# Grid af x-værdier
x_grid     <- seq(-1, 1, length.out = 300)
grid_df    <- data.frame(x = x_grid)
pred_grid  <- predict(rf_band, grid_df)$predictions

# Oracle
z <- qnorm(1 - alpha / 2)

# Gaussian
sigma_hat <- sqrt(mean((train_band$y - pred_train_band)^2))

# CP
k     <- ceiling((n_cal_band + 1) * (1 - alpha))
q_cp  <- sort(abs(cal_band$y - pred_cal_band))[k]

# RQ
res_cal_band <- cal_band$y - pred_cal_band
q_rq_low     <- quantile(res_cal_band, alpha / 2)
q_rq_high    <- quantile(res_cal_band, 1 - alpha / 2)

# Saml bånd
band_df <- data.frame(
  x           = x_grid,
  pred        = pred_grid,
  true_f      = f(x_grid),
  oracle_low  = f(x_grid) - z * sigma,
  oracle_high = f(x_grid) + z * sigma,
  gauss_low   = pred_grid - qnorm(1 - alpha / 2) * sigma_hat,
  gauss_high  = pred_grid + qnorm(1 - alpha / 2) * sigma_hat,
  cp_low      = pred_grid - q_cp,
  cp_high     = pred_grid + q_cp,
  rq_low      = pred_grid + q_rq_low,
  rq_high     = pred_grid + q_rq_high
)

# Plot funktion
plot_band <- function(low_col, high_col, title, fill_col) {
  ggplot(band_df, aes(x = x)) +
    geom_ribbon(aes(ymin = .data[[low_col]], ymax = .data[[high_col]]),
                fill = fill_col, alpha = 0.3) +
    geom_line(aes(y = true_f), color = "black", linetype = "dashed") +
    geom_line(aes(y = pred), color = "grey40") +
    geom_point(data = data.frame(x = x_show, y = y_show),
               aes(x = x, y = y), size = 0.8, alpha = 0.4) +
    labs(title = title, x = expression(x), y = expression(y)) +
    theme_bw()
}

library(patchwork)
plot_band("oracle_low", "oracle_high", "Oracle",   "turquoise") +
  plot_band("gauss_low",  "gauss_high",  "Gaussian", "green4") +
  plot_band("cp_low",     "cp_high",     "CP",       "red") +
  plot_band("rq_low",     "rq_high",     "RQ",       "purple")
