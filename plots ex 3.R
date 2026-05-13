# ============================================================
# Nonlinear model:
# Prediction intervals
# ============================================================

library(ranger)
library(ggplot2)
library(dplyr)

# =============================
# Setup
# =============================
n <- 80

alpha <- 0.1
z <- qnorm(1 - alpha/2)

# =============================
# True regression function
# =============================
f <- function(x) {
  x*exp(x)*sin(2*pi*x)
}
x <- runif(n, -1, 1)
plot(x,f(x))

sigma <- 1

# =============================
# Data generation
# =============================
gen_data <- function(n) {
  
  x <- runif(n, -1, 1)
  
  y <- f(x) + rnorm(n, sd = sigma)
  
  data.frame(x = x, y = y)
}

# split
n_train <- n/2
n_cal   <- n/2

train <- gen_data(n_train)
cal   <- gen_data(n_cal)
test  <- gen_data(500)

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
  min.node.size = 10
)

# RF mean predictions
rf_mean_train <- predict(rf_fit, train)$predictions
rf_mean_cal   <- predict(rf_fit, cal)$predictions
rf_mean_test  <- predict(rf_fit, test)$predictions

# RF quantile predictions
rf_quant_test <- predict(
  rf_fit,
  test,
  type = "quantiles",
  quantiles = c(alpha/2, 1-alpha/2)
)$predictions

# ============================================================
# INTERVALS
# ============================================================

# =============================
# Oracle
# =============================
oracle_l <- f(test$x) - z*sigma
oracle_u <- f(test$x) + z*sigma

# =============================
# Gaussian LM
# =============================
sigma_hat <- sqrt(mean((train$y - lm_pred_train)^2))

gauss_l <- lm_pred_test - z*sigma_hat
gauss_u <- lm_pred_test + z*sigma_hat

# =============================
# CP RF
# =============================
k <- ceiling((n_cal + 1)*(1 - alpha))

q_cp <- sort(abs(cal$y - rf_mean_cal))[k]

cp_l <- rf_mean_test - q_cp
cp_u <- rf_mean_test + q_cp

# =============================
# RQ RF
# =============================
res_rf <- train$y - rf_mean_train

q_lo_rq <- quantile(res_rf, alpha/2)
q_hi_rq <- quantile(res_rf, 1 - alpha/2)

rq_l <- rf_mean_test + q_lo_rq
rq_u <- rf_mean_test + q_hi_rq

# =============================
# QR RF
# =============================
qr_l <- rf_quant_test[,1]
qr_u <- rf_quant_test[,2]

# ============================================================
# PLOT DATA
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
    method = "Gaussian (LM)"
  ),
  
  data.frame(
    x = test$x[ord],
    y = test$y[ord],
    l = cp_l[ord],
    u = cp_u[ord],
    method = "CP (RF)"
  ),
  
  data.frame(
    x = test$x[ord],
    y = test$y[ord],
    l = rq_l[ord],
    u = rq_u[ord],
    method = "RQ (RF)"
  ),
  
  data.frame(
    x = test$x[ord],
    y = test$y[ord],
    l = qr_l[ord],
    u = qr_u[ord],
    method = "QR (RF)"
  )
)

# ============================================================
# ORDER
# ============================================================

df_plot$method <- factor(
  df_plot$method,
  levels = c(
    "Oracle",
    "Gaussian (LM)",
    "CP (RF)",
    "RQ (RF)",
    "QR (RF)"
  )
)

# ============================================================
# COLORS
# ============================================================

fill_cols <- c(
  "Oracle"   = "turquoise",
  "Gaussian (LM)" = "green4",
  "CP (RF)"       = "red",
  "RQ (RF)"       = "purple",
  "QR (RF)"       = "blue"
)

# ============================================================
# PLOT
# ============================================================

ggplot(df_plot, aes(x = x)) +
  
  geom_ribbon(
    aes(ymin = l,
        ymax = u,
        fill = method),
    alpha = 0.3
  ) +
  
  geom_point(
    aes(y = y),
    size = 0.6,
    alpha = 0.35
  ) +
  
  geom_line(
    aes(y = f(x)),
    linetype = "dashed"
  ) +
  
  facet_wrap(~ method, ncol = 3) +
  
  scale_fill_manual(values = fill_cols) +
  
  labs(x = "x",
    y = "y"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "none"
  )


