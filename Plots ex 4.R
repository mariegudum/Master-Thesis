# ============================================================
# Heteroskedastic setting:
# Local interval behavior and conditional coverage
# ============================================================

library(ranger)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================
# Setup
# =============================
n <- 1280
n_train <- n/2
n_cal   <- n/2
n_test  <- 1000
alpha   <- 0.1

beta <- 2
a <- 1

f <- function(x) beta*x + a
sigma_fun <- function(x) 1 + x

# =============================
# Data generation
# =============================
gen_data <- function(n) {
  
  x <- runif(n, -0.99, 1)
  
  sigma_x <- sigma_fun(x)
  
  y <- f(x) + sigma_x * rnorm(n)
  
  data.frame(x = x, y = y)
}

train <- gen_data(n_train)
cal   <- gen_data(n_cal)
test  <- gen_data(n_test)

# ============================================================
# MODELS
# ============================================================

# =============================
# Linear model
# =============================
lm_fit <- lm(y ~ x, data = train)

lm_pred_train <- predict(lm_fit, train)
lm_pred_test  <- predict(lm_fit, test)

# =============================
# Random Forest
# =============================
rf_fit <- ranger(
  y ~ x,
  data = train,
  quantreg = TRUE,
  num.trees = 500,
  min.node.size = 50
)

# RF mean predictions
rf_mean_train <- predict(rf_fit, train)$predictions
rf_mean_cal   <- predict(rf_fit, cal)$predictions
rf_mean_test  <- predict(rf_fit, test)$predictions

# RF quantile predictions
rf_quant_cal <- predict(
  rf_fit,
  cal,
  type = "quantiles",
  quantiles = c(alpha/2, 1-alpha/2)
)$predictions

rf_quant_test <- predict(
  rf_fit,
  test,
  type = "quantiles",
  quantiles = c(alpha/2, 1-alpha/2)
)$predictions

# ============================================================
# INTERVALS
# ============================================================

z <- qnorm(1 - alpha/2)

# =============================
# Oracle
# =============================
sigma_x_test <- sigma_fun(test$x)

oracle_l <- f(test$x) - z * sigma_x_test
oracle_u <- f(test$x) + z * sigma_x_test

# =============================
# Gaussian LM
# =============================
sigma_hat <- sqrt(mean((train$y - lm_pred_train)^2))

gauss_l <- lm_pred_test - z * sigma_hat
gauss_u <- lm_pred_test + z * sigma_hat

# =============================
# Residual Quantiles (RQ RF)
# =============================
res_rf <- train$y - rf_mean_train

q_lo_rq_rf <- quantile(res_rf, alpha/2)
q_hi_rq_rf <- quantile(res_rf, 1 - alpha/2)

rq_l_rf <- rf_mean_test + q_lo_rq_rf
rq_u_rf <- rf_mean_test + q_hi_rq_rf

# =============================
# Conformal Prediction (CP RF)
# =============================
k <- ceiling((n_cal + 1)*(1 - alpha))

q_cp_rf <- sort(abs(cal$y - rf_mean_cal))[k]

cp_l_rf <- rf_mean_test - q_cp_rf
cp_u_rf <- rf_mean_test + q_cp_rf

# =============================
# Quantile Regression (QR RF)
# =============================
qr_l_rf <- rf_quant_test[,1]
qr_u_rf <- rf_quant_test[,2]

# ============================================================
# INTERVAL LENGTH AS FUNCTION OF x
# ============================================================

df_x <- data.frame(
  x = test$x,
  Oracle   = oracle_u - oracle_l,
  Gaussian = gauss_u - gauss_l,
  CP_RF    = cp_u_rf - cp_l_rf,
  RQ_RF    = rq_u_rf - rq_l_rf,
  QR_RF    = qr_u_rf - qr_l_rf
)

df_x_long <- df_x %>%
  pivot_longer(
    -x,
    names_to = "method",
    values_to = "length"
  )

df_x_long$method <- factor(
  df_x_long$method,
  levels = c(
    "Oracle",
    "Gaussian",
    "CP_RF",
    "RQ_RF",
    "QR_RF"
  )
)
method_cols <- c(
  "Oracle"   = "black",
  "Gaussian (LM)" = "#F8766D",
  "CP (RF)"    = "#53B400",
  "RQ (RF)"    = "#FB61D7",
  "QR (RF)"    = "#00B6EB"
)

df_x_long = df_x_long %>%
  mutate(
    method = recode(method,
                    "CP_RF" = "CP (RF)",
                    "Gaussian" = "Gaussian (LM)",
                    "QR_RF" = "QR (RF)",
                    "RQ_RF" = "RQ (RF)"
    )
  )

p1 <- ggplot(df_x_long,
       aes(x = x,
           y = length,
           color = method)) +
  
  geom_line(linewidth = 1) +
  
  geom_line(
    aes(y = 2*z*sigma_fun(x)),
    color = "black",
    linetype = "dashed"
  ) +
  scale_color_manual(values = method_cols) +
  labs(x = "x",
    y = "Interval length",
    color = "Method"
  ) +
  theme(legend.position = "none") +
  theme_minimal(base_size = 13)
p1 <- p1 +
  theme(legend.position = "none") +
  guides(
    color = "none",
    fill = "none"
  )
p1
# ============================================================
# INTERVAL PLOTS
# ============================================================

ord <- order(test$x)

df_plot <- bind_rows(
  
  data.frame(
    x = test$x[ord],
    y = test$y[ord],
    l = oracle_l[ord],
    u = oracle_u[ord],
    method = "Oracle"
  ),
  
  data.frame(
    x = test$x[ord],
    y = test$y[ord],
    l = gauss_l[ord],
    u = gauss_u[ord],
    method = "Gaussian"
  ),
  
  data.frame(
    x = test$x[ord],
    y = test$y[ord],
    l = cp_l_rf[ord],
    u = cp_u_rf[ord],
    method = "CP (RF)"
  ),
  
  data.frame(
    x = test$x[ord],
    y = test$y[ord],
    l = rq_l_rf[ord],
    u = rq_u_rf[ord],
    method = "RQ (RF)"
  ),
  
  data.frame(
    x = test$x[ord],
    y = test$y[ord],
    l = qr_l_rf[ord],
    u = qr_u_rf[ord],
    method = "QR (RF)"
  )
)

df_plot$method <- factor(
  df_plot$method,
  levels = c(
    "Oracle",
    "Gaussian",
    "CP (RF)",
    "RQ (RF)",
    "QR (RF)"
  )
)

ggplot(df_plot, aes(x = x)) +
  
  geom_ribbon(
    aes(ymin = l, ymax = u),
    fill = "grey70",
    alpha = 0.4
  ) +
  
  geom_point(
    aes(y = y),
    size = 0.6,
    alpha = 0.3
  ) +
  
  facet_wrap(~ method, ncol = 2) +
  
  labs(
    title = "Prediction intervals as function of x",
    x = "x",
    y = "y"
  ) +
  
  theme_minimal(base_size = 13)

# ============================================================
# CONDITIONAL COVERAGE AS FUNCTION OF x
# ============================================================

x_grid <- seq(-0.9, 0.9, length.out = 100)

B <- 50000

coverage_at_x <- function(x0) {
  
  sigma0 <- sigma_fun(x0)
  
  y_samples <- f(x0) + sigma0 * rnorm(B)
  
  # -------------------------
  # Oracle
  # -------------------------
  l_oracle <- f(x0) - z * sigma0
  u_oracle <- f(x0) + z * sigma0
  
  cov_oracle <- mean(
    y_samples >= l_oracle &
      y_samples <= u_oracle
  )
  
  # -------------------------
  # Gaussian
  # -------------------------
  mu_hat <- predict(
    lm_fit,
    data.frame(x = x0)
  )
  
  l_gauss <- mu_hat - z * sigma_hat
  u_gauss <- mu_hat + z * sigma_hat
  
  cov_gauss <- mean(
    y_samples >= l_gauss &
      y_samples <= u_gauss
  )
  
  # -------------------------
  # RF mean prediction
  # -------------------------
  pred_rf <- predict(
    rf_fit,
    data.frame(x = x0)
  )$predictions
  
  # -------------------------
  # CP RF
  # -------------------------
  l_cp <- pred_rf - q_cp_rf
  u_cp <- pred_rf + q_cp_rf
  
  cov_cp <- mean(
    y_samples >= l_cp &
      y_samples <= u_cp
  )
  
  # -------------------------
  # RQ RF
  # -------------------------
  l_rq <- pred_rf + q_lo_rq_rf
  u_rq <- pred_rf + q_hi_rq_rf
  
  cov_rq <- mean(
    y_samples >= l_rq &
      y_samples <= u_rq
  )
  
  # -------------------------
  # QR RF
  # -------------------------
  qr_pred <- predict(
    rf_fit,
    data.frame(x = x0),
    type = "quantiles",
    quantiles = c(alpha/2, 1-alpha/2)
  )$predictions
  
  l_qr <- qr_pred[1]
  u_qr <- qr_pred[2]
  
  cov_qr <- mean(
    y_samples >= l_qr &
      y_samples <= u_qr
  )
  
  data.frame(
    x = x0,
    method = c(
      "Oracle",
      "Gaussian",
      "CP_RF",
      "RQ_RF",
      "QR_RF"
    ),
    coverage = c(
      cov_oracle,
      cov_gauss,
      cov_cp,
      cov_rq,
      cov_qr
    )
  )
}

results_cond <- bind_rows(
  lapply(x_grid, coverage_at_x)
)

results_cond$method <- factor(
  results_cond$method,
  levels = c(
    "Oracle",
    "Gaussian",
    "CP_RF",
    "RQ_RF",
    "QR_RF"
  )
)

results_cond = results_cond %>%
  mutate(
    method = recode(method,
                    "CP_RF" = "CP (RF)",
                    "Gaussian" = "Gaussian (LM)",
                    "QR_RF" = "QR (RF)",
                    "RQ_RF" = "RQ (RF)"
    )
  )

ggplot(results_cond,
       aes(x = x,
           y = coverage,
           color = method)) +
  
  geom_line(linewidth = 1) +
  
  geom_hline(
    yintercept = 0.9,
    linetype = "dashed"
  ) +
  scale_color_manual(values = method_cols) +
  labs(x = "x",
    y = "Coverage"
  ) +
  
  theme_minimal(base_size = 13)

