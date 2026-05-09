# =============================
# Setup
# =============================
set.seed(1)

library(ggplot2)
library(tidyr)

# -----------------------------
# Parameters
# -----------------------------
n_list <- c(20, 80, 320, 1280)
M <- 200
n_test <- 2000

alpha <- 0.1
beta  <- 2
a     <- 1
sigma <- 1

f <- function(x) 2*x + sin(3*x)
z <- qnorm(1 - alpha/2)


# -----------------------------
# FAST Quantile Regression
# -----------------------------
qr_fit_fast <- function(x, y, tau) {
  
  X <- cbind(1, x)
  
  loss <- function(beta) {
    r <- y - X %*% beta
    sum((tau - (r < 0)) * r)
  }
  
  # smart initialization (OLS)
  beta_init <- coef(lm(y ~ x))
  
  optim(beta_init, loss, method = "BFGS")$par
}

qr_predict <- function(beta, x) {
  cbind(1, x) %*% beta
}

# -----------------------------
# Storage
# -----------------------------
results_all <- data.frame()

# -----------------------------
# Simulation loop
# -----------------------------
for (n in n_list) {
  
  cat("Running n =", n, "\n")
  
  n_train <- n/2
  n_cal   <- n/2
  
  cov <- len <- matrix(0, M, 5)
  colnames(cov) <- colnames(len) <- c("Oracle","Linear","Naive","QR","CP")
  
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
    # Fit linear model ONCE
    # -------------------------
    lm_fit <- lm(y ~ x, data = data.frame(x = train$x, y = train$y))
    
    pred_train <- predict(lm_fit)
    pred_cal   <- predict(lm_fit, newdata = data.frame(x = cal$x))
    pred_test  <- predict(lm_fit, newdata = data.frame(x = test$x))
    
    # FIX for predict naming
    pred_cal  <- coef(lm_fit)[1] + coef(lm_fit)[2] * cal$x
    pred_test <- coef(lm_fit)[1] + coef(lm_fit)[2] * test$x
    
    # -------------------------
    # (i) Oracle
    # -------------------------
    oracle_l <- f(test$x) - z * sigma
    oracle_u <- f(test$x) + z * sigma
    
    # -------------------------
    # (ii) Linear
    # -------------------------
    sigma_hat <- sqrt(mean((train$y - pred_train)^2))
    lin_l <- pred_test - z * sigma_hat
    lin_u <- pred_test + z * sigma_hat
    
    # -------------------------
    # Residuals
    # -------------------------
    res_train <- train$y - pred_train
    
    # -------------------------
    # (iii) Naive
    # -------------------------
    ql <- quantile(res_train, alpha/2)
    qu <- quantile(res_train, 1 - alpha/2)
    
    naive_l <- pred_test + ql
    naive_u <- pred_test + qu
    
    # -------------------------
    # (iv) TRUE QR (optim)
    # -------------------------
    beta_low  <- qr_fit_fast(train$x, train$y, alpha/2)
    beta_high <- qr_fit_fast(train$x, train$y, 1 - alpha/2)
    
    qr_l <- qr_predict(beta_low, test$x)
    qr_u <- qr_predict(beta_high, test$x)
    
    # -------------------------
    # (v) Conformal
    # -------------------------
    res_cal <- abs(cal$y - pred_cal)
    k <- ceiling((n_cal + 1) * (1 - alpha))
    
    q_cp <- sort(res_cal)[k]
    
    cp_l <- pred_test - q_cp
    cp_u <- pred_test + q_cp
    
    # -------------------------
    # Metrics (vectorized)
    # -------------------------
    coverage <- function(l, u) mean(test$y >= l & test$y <= u)
    length   <- function(l, u) mean(u - l)
    
    cov[m,] <- c(
      coverage(oracle_l, oracle_u),
      coverage(lin_l, lin_u),
      coverage(naive_l, naive_u),
      coverage(qr_l, qr_u),
      coverage(cp_l, cp_u)
    )
    
    len[m,] <- c(
      length(oracle_l, oracle_u),
      length(lin_l, lin_u),
      length(naive_l, naive_u),
      length(qr_l, qr_u),
      length(cp_l, cp_u)
    )
  }
  
  results_all <- rbind(results_all, data.frame(
    n = n,
    Method = rep(colnames(cov), each = M),
    Coverage = as.vector(cov),
    Length   = as.vector(len)
  ))
}

# -----------------------------
# Order methods
# -----------------------------
results_all$Method <- factor(
  results_all$Method,
  levels = c("Oracle","Linear","Naive","QR","CP")
)

# -----------------------------
# Reshape
# -----------------------------
results_long <- pivot_longer(
  results_all,
  cols = c(Coverage, Length),
  names_to = "Metric",
  values_to = "Value"
)

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(results_long,
            aes(x = Method, y = Value, fill = factor(n))) +
  geom_boxplot(outlier.size = 0.5,
               position = position_dodge(0.75)) +
  facet_wrap(~ Metric, scales = "free_y") +
  geom_hline(
    data = data.frame(Metric = "Coverage"),
    aes(yintercept = 1 - alpha),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  theme_bw()

print(p)


curve(f(x), from=-1, to=1, col="black", lwd=2)
points(test$x, test$y, col=rgb(0,0,0,0.1))
