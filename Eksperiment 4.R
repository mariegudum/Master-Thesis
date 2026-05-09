set.seed(1)

library(dplyr)
library(ggplot2)
library(ranger)

# =============================
# Setup
# =============================
n_list <- c(20, 80, 320, 1280)
M      <- 200
n_test <- 100000
alpha  <- 0.1

beta <- 2; a <- 1
f <- function(x) beta*x + a

# Heteroskedastic std
sigma_fun <- function(x) 1 + x

# =============================
# Progress
# =============================
total_iters <- length(n_list) * M
counter <- 1
start_time <- Sys.time()

results_all_4 <- vector("list", total_iters)

# =============================
# Simulation
# =============================
for (n in n_list) {
  
  n_train <- n/2
  n_cal   <- n/2
  
  # per-n grid (QR tuning)
  if (n <= 80) {
    node_grid <- c(5, 10, 20)
  } else if (n <= 320) {
    node_grid <- c(20, 40, 80)
  } else {
    node_grid <- c(50, 80, 120)
  }
  
  for (m in 1:M) {
    
    # ===== Progress =====
    if (counter %% 5 == 0 || counter == 1) {
      elapsed <- as.numeric(Sys.time() - start_time, units="secs")
      progress <- counter / total_iters
      eta <- if (counter > 1) elapsed/(counter-1)*(total_iters-counter+1) else NA
      
      cat(sprintf(
        "Iter %d/%d (%.1f%%) | n=%d, m=%d | ETA=%.1fs\n",
        counter, total_iters, 100*progress, n, m, eta
      ))
    }
    
    # =========================
    # Data (heteroskedastic)
    # =========================
    gen_data <- function(n) {
      x <- runif(n, -1, 1)
      sigma_x <- sigma_fun(x)
      y <- f(x) + sigma_x * rnorm(n)
      data.frame(x = x, y = y)
    }
    
    train <- gen_data(n_train)
    cal   <- gen_data(n_cal)
    test  <- gen_data(n_test)
    
    # =========================
    # LM (misspecified variance)
    # =========================
    lm_fit <- lm(y ~ x, data = train)
    lm_pred_train <- predict(lm_fit, train)
    lm_pred_cal   <- predict(lm_fit, cal)
    lm_pred_test  <- predict(lm_fit, test)
    
    # =========================
    # RF (tuned baseline)
    # =========================
    rf_fit <- ranger(
      y ~ x,
      data = train,
      num.trees = 500,
      min.node.size = 30,
      sample.fraction = 0.8
    )
    
    rf_pred_cal  <- predict(rf_fit, cal)$predictions
    rf_pred_test <- predict(rf_fit, test)$predictions
    
    # =========================
    # Oracle (x-dependent)
    # =========================
    z <- qnorm(1 - alpha/2)
    sigma_x_test <- sigma_fun(test$x)
    
    oracle_l <- f(test$x) - z * sigma_x_test
    oracle_u <- f(test$x) + z * sigma_x_test
    
    # =========================
    # Gaussian LM
    # =========================
    sigma_hat <- sqrt(mean((train$y - lm_pred_train)^2))
    
    gauss_l <- lm_pred_test - z * sigma_hat
    gauss_u <- lm_pred_test + z * sigma_hat
    
    # =========================
    # RQ
    # =========================
    rq_l_lm <- lm_pred_test + quantile(cal$y - lm_pred_cal, alpha/2)
    rq_u_lm <- lm_pred_test + quantile(cal$y - lm_pred_cal, 1 - alpha/2)
    
    rq_l_rf <- rf_pred_test + quantile(cal$y - rf_pred_cal, alpha/2)
    rq_u_rf <- rf_pred_test + quantile(cal$y - rf_pred_cal, 1 - alpha/2)
    
    # =========================
    # CP
    # =========================
    k <- ceiling((n_cal + 1)*(1 - alpha))
    
    q_lm <- sort(abs(cal$y - lm_pred_cal))[k]
    cp_l_lm <- lm_pred_test - q_lm
    cp_u_lm <- lm_pred_test + q_lm
    
    q_rf <- sort(abs(cal$y - rf_pred_cal))[k]
    cp_l_rf <- rf_pred_test - q_rf
    cp_u_rf <- rf_pred_test + q_rf
    
    # =========================
    # QR LM
    # =========================
    X_train <- model.matrix(~x, train)
    X_test  <- model.matrix(~x, test)
    
    rho <- function(u, tau) u*(tau - (u < 0))
    
    fit_qr <- function(X, y, tau) {
      optim(rep(0, ncol(X)),
            function(b) sum(rho(y - X %*% b, tau)),
            method="BFGS")$par
    }
    
    b_l <- fit_qr(X_train, train$y, alpha/2)
    b_u <- fit_qr(X_train, train$y, 1 - alpha/2)
    
    qr_l_lm <- X_test %*% b_l
    qr_u_lm <- X_test %*% b_u
    
    # =========================
    # QR RF (per-n tuning)
    # =========================
    best_score <- Inf
    best_model <- NULL
    
    for (node_size in node_grid) {
      
      rf_tmp <- ranger(
        y ~ x,
        data = train,
        quantreg = TRUE,
        num.trees = 300,
        min.node.size = node_size,
        sample.fraction = 0.7
      )
      
      pred_cal <- predict(
        rf_tmp, cal,
        type = "quantiles",
        quantiles = c(alpha/2, 1 - alpha/2)
      )$predictions
      
      u1 <- cal$y - pred_cal[,1]
      u2 <- cal$y - pred_cal[,2]
      
      loss <- mean(
        u1 * (alpha/2 - (u1 < 0)) +
          u2 * ((1 - alpha/2) - (u2 < 0))
      )
      
      if (loss < best_score) {
        best_score <- loss
        best_model <- rf_tmp
      }
    }
    
    qr_pred <- predict(
      best_model, test,
      type="quantiles",
      quantiles=c(alpha/2, 1 - alpha/2)
    )$predictions
    
    qr_l_rf <- qr_pred[,1]
    qr_u_rf <- qr_pred[,2]
    
    # =========================
    # Metrics
    # =========================
    cov <- function(l,u) mean(test$y >= l & test$y <= u)
    wid <- function(l,u) mean(u-l)
    
    results_all_4[[counter]] <- data.frame(
      n=n, m=m,
      method=c("Oracle",
               "Gaussian_LM",
               "RQ_LM","RQ_RF_tuned",
               "QR_LM","QR_RF_tuned",
               "CP_LM","CP_RF_tuned"),
      coverage=c(
        cov(oracle_l,oracle_u),
        cov(gauss_l,gauss_u),
        cov(rq_l_lm,rq_u_lm),
        cov(rq_l_rf,rq_u_rf),
        cov(qr_l_lm,qr_u_lm),
        cov(qr_l_rf,qr_u_rf),
        cov(cp_l_lm,cp_u_lm),
        cov(cp_l_rf,cp_u_rf)
      ),
      length=c(
        wid(oracle_l,oracle_u),
        wid(gauss_l,gauss_u),
        wid(rq_l_lm,rq_u_lm),
        wid(rq_l_rf,rq_u_rf),
        wid(qr_l_lm,qr_u_lm),
        wid(qr_l_rf,qr_u_rf),
        wid(cp_l_lm,cp_u_lm),
        wid(cp_l_rf,cp_u_rf)
      )
    )
    
    counter <- counter + 1
  }
}

results_all_4 <- bind_rows(results_all_4)

total_time <- Sys.time() - start_time
cat("\nTotal runtime:", round(as.numeric(total_time, units="secs"),1), "seconds\n")


results_all_4 %>%
  group_by(n, method) %>%
  summarise(mean_cov = mean(coverage), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_cov, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  scale_x_log10() +
  labs(title = "Coverage (heteroskedastic noise)", y = "Coverage") +
  theme_minimal()


oracle_mean_4 <- results_all_4 %>%
  filter(method == "Oracle") %>%
  group_by(n) %>%
  summarise(oracle_length = mean(length), .groups="drop")

rel_results_4 <- results_all_4 %>%
  left_join(oracle_mean_4, by="n") %>%
  mutate(rel_length = length / oracle_length)

rel_results_4 %>%
  filter(method != "Oracle") %>%
  group_by(n, method) %>%
  summarise(mean_rel = mean(rel_length), .groups="drop") %>%
  ggplot(aes(x = n, y = mean_rel, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(title = "Relative interval length (heteroskedastic noise)",
       y = "Length / Oracle") +
  theme_minimal()



# længde som funktion af x


# =============================
# Setup
# =============================
n <- 1280
n_train <- n/2
n_cal   <- n/2
n_test  <- 1000
alpha   <- 0.1

beta <- 2; a <- 1
f <- function(x) beta*x + a
sigma_fun <- function(x) 1 + x

# =============================
# Data
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

# =============================
# Models
# =============================

# LM
lm_fit <- lm(y ~ x, data = train)
lm_pred_cal  <- predict(lm_fit, cal)
lm_pred_test <- predict(lm_fit, test)

# RF (quantile)
rf_fit <- ranger(
  y ~ x,
  data = train,
  quantreg = TRUE,
  num.trees = 500,
  min.node.size = 50
)

rf_pred_cal <- predict(rf_fit, cal,
                       type="quantiles",
                       quantiles=c(alpha/2, 1-alpha/2))$predictions

rf_pred_test <- predict(rf_fit, test,
                        type="quantiles",
                        quantiles=c(alpha/2, 1-alpha/2))$predictions

# =============================
# Intervals
# =============================

z <- qnorm(1 - alpha/2)

# Oracle
sigma_x_test <- sigma_fun(test$x)
oracle_l <- f(test$x) - z * sigma_x_test
oracle_u <- f(test$x) + z * sigma_x_test

# Gaussian LM
sigma_hat <- sqrt(mean((train$y - predict(lm_fit, train))^2))
gauss_l <- lm_pred_test - z * sigma_hat
gauss_u <- lm_pred_test + z * sigma_hat

# CP RF
k <- ceiling((n_cal + 1)*(1 - alpha))
q_rf <- sort(abs(cal$y - rowMeans(rf_pred_cal)))[k]

cp_l_rf <- rowMeans(rf_pred_test) - q_rf
cp_u_rf <- rowMeans(rf_pred_test) + q_rf

# QR RF
qr_l_rf <- rf_pred_test[,1]
qr_u_rf <- rf_pred_test[,2]

# =============================
# Build df_x
# =============================
df_x <- data.frame(
  x = test$x,
  Oracle = oracle_u - oracle_l,
  Gaussian = gauss_u - gauss_l,
  CP_RF = cp_u_rf - cp_l_rf,
  QR_RF = qr_u_rf - qr_l_rf
)

# =============================
# Long format
# =============================
df_x_long <- df_x %>%
  pivot_longer(-x, names_to = "method", values_to = "length")

# =============================
# Plot
# =============================
ggplot(df_x_long, aes(x = x, y = length, color = method)) +
  geom_line(linewidth = 1) +
  geom_line(aes(y = 2*z*sigma_fun(x)),
            color = "black", linetype = "dashed") +
  labs(
    x = "x",
    y = "Interval length",
    title = "Interval length as a function of x",
    color = "Method"
  ) +
  theme_minimal(base_size = 13)


# =============================
# Sort data (important for lines)
# =============================
ord <- order(test$x)

df_plot <- bind_rows(
  data.frame(
    x = test$x[ord], y = test$y[ord],
    l = oracle_l[ord], u = oracle_u[ord],
    method = "Oracle"
  ),
  data.frame(
    x = test$x[ord], y = test$y[ord],
    l = gauss_l[ord], u = gauss_u[ord],
    method = "Gaussian"
  ),
  data.frame(
    x = test$x[ord], y = test$y[ord],
    l = cp_l_rf[ord], u = cp_u_rf[ord],
    method = "CP (RF)"
  ),
  data.frame(
    x = test$x[ord], y = test$y[ord],
    l = qr_l_rf[ord], u = qr_u_rf[ord],
    method = "QR (RF)"
  )
)

# =============================
# Plot
# =============================
ggplot(df_plot, aes(x = x)) +
  
  geom_ribbon(aes(ymin = l, ymax = u),
              fill = "grey70", alpha = 0.4) +
  
  geom_point(aes(y = y), size = 0.6, alpha = 0.3) +
  
  facet_wrap(~ method, ncol = 2) +
  
  labs(
    title = "Prediction intervals as function of x",
    x = "x",
    y = "y"
  ) +
  
  theme_minimal(base_size = 13)




# =============================
# Conditional coverage test
# =============================
x_grid <- seq(-0.9, 0.9, length.out = 100)
B <- 200000  
z <- qnorm(1 - alpha/2)

# funktion til coverage ved fast x
coverage_at_x <- function(x0) {
  
  sigma0 <- sigma_fun(x0)
  y_samples <- f(x0) + sigma0 * rnorm(B)
  
  # --- Oracle ---
  l_oracle <- f(x0) - z * sigma0
  u_oracle <- f(x0) + z * sigma0
  cov_oracle <- mean(y_samples >= l_oracle & y_samples <= u_oracle)
  
  # --- Gaussian LM ---
  mu_hat <- predict(lm_fit, data.frame(x = x0))
  l_gauss <- mu_hat - z * sigma_hat
  u_gauss <- mu_hat + z * sigma_hat
  cov_gauss <- mean(y_samples >= l_gauss & y_samples <= u_gauss)
  
  # --- CP RF ---
  pred_rf <- predict(rf_fit, data.frame(x = x0))$predictions
  l_cp <- pred_rf - q_rf
  u_cp <- pred_rf + q_rf
  cov_cp <- mean(y_samples >= l_cp & y_samples <= u_cp)
  
  # --- QR RF ---
  qr_pred <- predict(
    rf_fit,
    data.frame(x = x0),
    type = "quantiles",
    quantiles = c(alpha/2, 1-alpha/2)
  )$predictions
  
  l_qr <- qr_pred[1]
  u_qr <- qr_pred[2]
  cov_qr <- mean(y_samples >= l_qr & y_samples <= u_qr)
  
  data.frame(
    x = x0,
    method = c("Oracle", "Gaussian", "CP_RF", "QR_RF"),
    coverage = c(cov_oracle, cov_gauss, cov_cp, cov_qr)
  )
}

# =============================
# Kør for flere x
# =============================
results_cond <- bind_rows(lapply(x_grid, coverage_at_x))

# =============================
# Plot
# =============================
ggplot(results_cond, aes(x = x, y = coverage, color = method)) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  labs(
    title = "Conditional coverage at fixed x",
    y = "Coverage"
  ) +
  theme_minimal()


# =============================
# Grid for x
# =============================
x_grid <- seq(-0.9, 0.9, length.out = 100)

B <- 2000
z <- qnorm(1 - alpha/2)

coverage_at_x <- function(x0) {
  
  sigma0 <- sigma_fun(x0)
  y_samples <- f(x0) + sigma0 * rnorm(B)
  
  # Oracle
  l_oracle <- f(x0) - z * sigma0
  u_oracle <- f(x0) + z * sigma0
  cov_oracle <- mean(y_samples >= l_oracle & y_samples <= u_oracle)
  
  # Gaussian
  mu_hat <- predict(lm_fit, data.frame(x = x0))
  l_gauss <- mu_hat - z * sigma_hat
  u_gauss <- mu_hat + z * sigma_hat
  cov_gauss <- mean(y_samples >= l_gauss & y_samples <= u_gauss)
  
  # CP RF
  pred_rf <- predict(rf_fit, data.frame(x = x0))$predictions
  l_cp <- pred_rf - q_rf
  u_cp <- pred_rf + q_rf
  cov_cp <- mean(y_samples >= l_cp & y_samples <= u_cp)
  
  # QR RF
  qr_pred <- predict(
    rf_fit,
    data.frame(x = x0),
    type = "quantiles",
    quantiles = c(alpha/2, 1-alpha/2)
  )$predictions
  
  l_qr <- qr_pred[1]
  u_qr <- qr_pred[2]
  cov_qr <- mean(y_samples >= l_qr & y_samples <= u_qr)
  
  data.frame(
    x = x0,
    method = c("Oracle", "Gaussian", "CP_RF", "QR_RF"),
    coverage = c(cov_oracle, cov_gauss, cov_cp, cov_qr)
  )
}

# =============================
# Compute
# =============================
results_cond <- bind_rows(lapply(x_grid, coverage_at_x))

ggplot(results_cond, aes(x = x, y = coverage, color = method)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  labs(
    title = "Conditional coverage as function of x",
    y = "Coverage"
  ) +
  theme_minimal(base_size = 13)
