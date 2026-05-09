# =============================
# Setup
# =============================
set.seed(1)

library(ggplot2)
library(dplyr)

# -----------------------------
# Parameters
# -----------------------------
n_list <- c(20, 80, 320, 1280)
M <- 200
n_test <- 100000  
alpha <- 0.1
beta  <- 2
a     <- 1
sigma <- 1

f <- function(x) beta * x + a

# -----------------------------
# Storage
# -----------------------------
results_all <- data.frame()

# =============================
# Simulation
# =============================
for (n in n_list) {
  
  cat("Running n =", n, "\n")
  
  n_train <- n/2
  n_cal   <- n/2
  
  for (m in 1:M) {
    
    # -------------------------
    # Data
    # -------------------------
    gen_data <- function(n) {
      x <- runif(n, -1, 1)
      eps <- rnorm(n, 0, sigma)
      y <- f(x) + eps
      data.frame(x = x, y = y)
    }
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # -------------------------
    # Fit model
    # -------------------------
    lm_fit <- lm(y ~ x, data = train)
    
    pred_train <- predict(lm_fit, train)
    pred_cal   <- predict(lm_fit, cal)
    pred_test  <- predict(lm_fit, test)
    
    # =========================
    # (1) Oracle
    # =========================
    z <- qnorm(1 - alpha/2)
    oracle_lower <- f(test$x) - z * sigma
    oracle_upper <- f(test$x) + z * sigma
    
    # =========================
    # (2) Gaussian 
    # =========================
    sigma_hat <- sqrt(sum((train$y - pred_train)^2) / (n_train - 2))
    
    X_mat <- model.matrix(lm_fit)
    XtX_inv <- solve(t(X_mat) %*% X_mat)
    
    x_test_mat <- model.matrix(~ x, data = test)
    
    se_pred <- sqrt(
      sigma_hat^2 * (1 + rowSums((x_test_mat %*% XtX_inv) * x_test_mat))
    )
    
    t_val <- qt(1 - alpha/2, df = n_train - 2)
    
    gauss_lower <- pred_test - t_val * se_pred
    gauss_upper <- pred_test + t_val * se_pred
    
    # =========================
    # (3) Residual quantile
    # =========================
    res_cal <- cal$y - pred_cal
    
    q_low  <- quantile(res_cal, alpha/2)
    q_high <- quantile(res_cal, 1 - alpha/2)
    
    qr_lower <- pred_test + q_low
    qr_upper <- pred_test + q_high
    
    # =========================
    # (4) Conformal
    # =========================
    abs_res_cal <- abs(cal$y - pred_cal)
    
    k <- ceiling((n_cal + 1) * (1 - alpha))
    q_cp <- sort(abs_res_cal)[k]
    
    cp_lower <- pred_test - q_cp
    cp_upper <- pred_test + q_cp
    
    # -------------------------
    # Metrics
    # -------------------------
    coverage <- function(l, u, y) mean(y >= l & y <= u)
    length   <- function(l, u) mean(u - l)
    
    # store
    results_all <- rbind(results_all, data.frame(
      n = n,
      method = c("Oracle", "Gaussian", "QR", "CP"),
      coverage = c(
        coverage(oracle_lower, oracle_upper, test$y),
        coverage(gauss_lower, gauss_upper, test$y),
        coverage(qr_lower, qr_upper, test$y),
        coverage(cp_lower, cp_upper, test$y)
      ),
      length = c(
        length(oracle_lower, oracle_upper),
        length(gauss_lower, gauss_upper),
        length(qr_lower, qr_upper),
        length(cp_lower, cp_upper)
      )
    ))
  }
}

# =============================
# Post-processing
# =============================
results_all$n <- factor(results_all$n)

# compute oracle length per n
oracle_length <- results_all |>
  filter(method == "Oracle") |>
  group_by(n) |>
  summarise(L_oracle = mean(length))

results_all <- left_join(results_all, oracle_length, by = "n")

# sqrt(n) scaling
results_all <- results_all |>
  mutate(length_scaled = sqrt(as.numeric(as.character(n))) * (length - L_oracle))

# =============================
# Plot: Coverage
# =============================
ggplot(results_all,
       aes(x = n, y = coverage, fill = method)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
  labs(y = "Coverage", x = "n") +
  theme_bw()

# =============================
# Plot: Length (raw)
# =============================
ggplot(results_all,
       aes(x = n, y = length, fill = method)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(y = "Length", x = "n") +
  theme_bw()

# =============================
# Plot: Length (scaled)
# =============================
ggplot(results_all,
       aes(x = n, y = length_scaled, fill = method)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(y = expression(sqrt(n)*(L[n] - L^"*")), x = "n") +
  theme_bw()

# (A) Raw scaling
results_all <- results_all |>
  mutate(length_scaled_raw = sqrt(as.numeric(as.character(n))) * length)

# plot
ggplot(results_all,
       aes(x = n, y = length_scaled_raw, fill = method)) +
  geom_boxplot() +
  labs(y = expression(sqrt(n)*L[n]), x = "n") +
  theme_bw()
