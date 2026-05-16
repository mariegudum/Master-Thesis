# =========================================================
# FULL SIMULATION — ALL METHODS IN ONE LOOP
# =========================================================

library(ranger)
library(ggplot2)
library(dplyr)

set.seed(1)

# =========================================================
# Setup
# =========================================================

n_list <- c(20, 80, 320, 1280)

M       <- 200
n_test  <- 100000
alpha   <- 0.1

beta  <- 2
a     <- 1
sigma <- 1

f <- function(x) beta*x + a

z <- qnorm(1 - alpha/2)

# =========================================================
# Data generator
# =========================================================

gen_data <- function(n) {
  x <- runif(n, -1, 1)
  
  data.frame(
    x = x,
    y = f(x) + rnorm(n, 0, sigma)
  )
}

# =========================================================
# Quantile regression (linear)
# =========================================================

rho <- function(u, tau) {
  u * (tau - (u < 0))
}

fit_qr <- function(X, y, tau) {
  
  optim(
    rep(0, ncol(X)),
    function(b) sum(rho(y - X %*% b, tau)),
    method = "BFGS"
  )$par
}

# =========================================================
# RF tuning depending on n
# =========================================================

get_rf_params <- function(n) {
  
  if (n <= 80) {
    
    list(
      num.trees      = 300,
      min.node.size  = 10,
      sample.fraction = 1
    )
    
  } else if (n <= 320) {
    
    list(
      num.trees      = 500,
      min.node.size  = 30,
      sample.fraction = 0.8
    )
    
  } else {
    
    list(
      num.trees      = 500,
      min.node.size  = 80,
      sample.fraction = 0.7
    )
  }
}

# =========================================================
# Metrics
# =========================================================

coverage <- function(l, u, y) {
  mean(y >= l & y <= u)
}

width <- function(l, u) {
  mean(u - l)
}

# =========================================================
# Pre-allocation
# =========================================================

total_iters <- length(n_list) * M

results <- vector("list", total_iters)

counter <- 1

start_time <- Sys.time()

# =========================================================
# MAIN LOOP
# =========================================================

for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  rf_par <- get_rf_params(n)
  
  for (m in 1:M) {
    
    # -----------------------------------------------------
    # Progress
    # -----------------------------------------------------
    
    if (counter %% 5 == 0) {
      
      elapsed <- as.numeric(
        Sys.time() - start_time,
        units = "secs"
      )
      
      eta <- elapsed/(counter-1) * (total_iters-counter+1)
      
      cat(sprintf(
        "Iter %d/%d | n=%d | ETA=%.1fs\n",
        counter,
        total_iters,
        n,
        eta
      ))
    }
    
    # -----------------------------------------------------
    # Data
    # -----------------------------------------------------
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # =====================================================
    # LM
    # =====================================================
    
    lm_fit <- lm(y ~ x, data = train)
    
    lm_pred_train <- predict(lm_fit, train)
    lm_pred_cal   <- predict(lm_fit, cal)
    lm_pred_test  <- predict(lm_fit, test)
    
    # =====================================================
    # RF
    # =====================================================
    
    rf_fit <- ranger(
      y ~ x,
      data = train,
      num.trees       = rf_par$num.trees,
      min.node.size   = rf_par$min.node.size,
      sample.fraction = rf_par$sample.fraction,
      num.threads     = parallel::detectCores()
    )
    
    rf_pred_cal <- predict(rf_fit, cal)$predictions
    rf_pred_test <- predict(rf_fit, test)$predictions
    
    # =====================================================
    # RF QUANTILE
    # =====================================================
    
    rf_qr_fit <- ranger(
      y ~ x,
      data = train,
      quantreg = TRUE,
      num.trees       = rf_par$num.trees,
      min.node.size   = rf_par$min.node.size,
      sample.fraction = rf_par$sample.fraction,
      num.threads     = parallel::detectCores()
    )
    
    qr_rf_pred <- predict(
      rf_qr_fit,
      test,
      type = "quantiles",
      quantiles = c(alpha/2, 1-alpha/2)
    )$predictions
    
    qr_rf_l <- qr_rf_pred[,1]
    qr_rf_u <- qr_rf_pred[,2]
    
    # =====================================================
    # ORACLE
    # =====================================================
    
    oracle_l <- f(test$x) - z*sigma
    oracle_u <- f(test$x) + z*sigma
    
    # =====================================================
    # GAUSSIAN LM
    # =====================================================
    
    sigma_hat <- sqrt(
      mean((train$y - lm_pred_train)^2)
    )
    
    gauss_l <- lm_pred_test - z*sigma_hat
    gauss_u <- lm_pred_test + z*sigma_hat
    
    # =====================================================
    # RQ
    # =====================================================
    
    # ---------- LM ----------
    
    rq_l_lm <- lm_pred_test +
      quantile(cal$y - lm_pred_cal, alpha/2)
    
    rq_u_lm <- lm_pred_test +
      quantile(cal$y - lm_pred_cal, 1-alpha/2)
    
    # ---------- RF ----------
    
    q_low_rf <- quantile(
      cal$y - rf_pred_cal,
      alpha/2
    )
    
    q_high_rf <- quantile(
      cal$y - rf_pred_cal,
      1-alpha/2
    )
    
    rq_l_rf <- rf_pred_test + q_low_rf
    rq_u_rf <- rf_pred_test + q_high_rf
    
    # =====================================================
    # QR
    # =====================================================
    
    # ---------- LM ----------
    
    X_cal  <- model.matrix(~x, cal)
    X_test <- model.matrix(~x, test)
    
    b_l <- fit_qr(X_cal, cal$y, alpha/2)
    b_u <- fit_qr(X_cal, cal$y, 1-alpha/2)
    
    qr_l_lm <- X_test %*% b_l
    qr_u_lm <- X_test %*% b_u
    
    # ---------- RF ----------
    
    qr_l_rf <- qr_rf_l
    qr_u_rf <- qr_rf_u
    
    # =====================================================
    # CP
    # =====================================================
    
    k <- ceiling((n_cal + 1)*(1-alpha))
    
    # ---------- LM ----------
    
    q_lm <- sort(abs(cal$y - lm_pred_cal))[k]
    
    cp_l_lm <- lm_pred_test - q_lm
    cp_u_lm <- lm_pred_test + q_lm
    
    # ---------- RF ----------
    
    q_rf <- sort(abs(cal$y - rf_pred_cal))[k]
    
    cp_l_rf <- rf_pred_test - q_rf
    cp_u_rf <- rf_pred_test + q_rf
    
    # =====================================================
    # SAVE RESULTS
    # =====================================================
    
    results[[counter]] <- data.frame(
      
      n = n,
      m = m,
      
      method = c(
        "Oracle",
        "Gaussian_LM",
        "RQ_LM",
        "RQ_RF_tuned",
        "QR_LM",
        "QR_RF_tuned",
        "CP_LM",
        "CP_RF_tuned"
      ),
      
      coverage = c(
        
        coverage(oracle_l, oracle_u, test$y),
        
        coverage(gauss_l, gauss_u, test$y),
        
        coverage(rq_l_lm, rq_u_lm, test$y),
        coverage(rq_l_rf, rq_u_rf, test$y),
        
        coverage(qr_l_lm, qr_u_lm, test$y),
        coverage(qr_l_rf, qr_u_rf, test$y),
        
        coverage(cp_l_lm, cp_u_lm, test$y),
        coverage(cp_l_rf, cp_u_rf, test$y)
      ),
      
      length = c(
        
        width(oracle_l, oracle_u),
        
        width(gauss_l, gauss_u),
        
        width(rq_l_lm, rq_u_lm),
        width(rq_l_rf, rq_u_rf),
        
        width(qr_l_lm, qr_u_lm),
        width(qr_l_rf, qr_u_rf),
        
        width(cp_l_lm, cp_u_lm),
        width(cp_l_rf, cp_u_rf)
      )
    )
    
    counter <- counter + 1
  }
}

# =========================================================
# Combine
# =========================================================

results_all <- bind_rows(results)

# =========================================================
# Relative lengths
# =========================================================

oracle_mean <- results_all %>%
  filter(method == "Oracle") %>%
  group_by(n) %>%
  summarise(
    oracle_length = mean(length),
    .groups = "drop"
  )

rel_results <- results_all %>%
  left_join(oracle_mean, by = "n") %>%
  mutate(
    rel_length = length / oracle_length
  )

# =========================================================
# COVERAGE PLOT
# =========================================================

coverage_plot <- results_all %>%
  
  filter(method != "Oracle") %>%
  
  group_by(n, method) %>%
  
  summarise(
    mean_cov = mean(coverage),
    .groups = "drop"
  ) %>%
  
  ggplot(
    aes(x = n, y = mean_cov, color = method)
  ) +
  
  geom_line() +
  geom_point() +
  
  geom_hline(
    yintercept = 1-alpha,
    linetype = "dashed"
  ) +
  
  scale_x_log10() +
  
  labs(
    title = "Coverage vs sample size",
    y = "Coverage"
  ) +
  
  theme_minimal()

print(coverage_plot)

# =========================================================
# RELATIVE LENGTH PLOT
# =========================================================

width_plot <- rel_results %>%
  
  filter(method != "Oracle") %>%
  
  group_by(n, method) %>%
  
  summarise(
    mean_rel = mean(rel_length),
    .groups = "drop"
  ) %>%
  
  ggplot(
    aes(x = n, y = mean_rel, color = method)
  ) +
  
  geom_line() +
  geom_point() +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  
  scale_x_log10() +
  
  labs(
    title = "Relative interval length",
    y = "Length / Oracle length"
  ) +
  
  theme_minimal()

print(width_plot)

# =========================================================
# COVERAGE DENSITY
# =========================================================

ggplot(
  results_all %>% filter(method != "Oracle"),
  aes(x = coverage, fill = method)
) +
  geom_density(alpha = 0.3) +
  geom_vline(
    xintercept = 1-alpha,
    linetype = "dashed"
  ) +
  labs(
    title = "Coverage distribution"
  ) +
  theme_minimal()

# =========================================================
# Runtime
# =========================================================

runtime <- Sys.time() - start_time

cat(
  "\nTotal runtime:",
  round(as.numeric(runtime, units="mins"),2),
  "minutes\n"
)