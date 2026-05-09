head(ENB2012_data)

summary(ENB2012_data)
cor(ENB2012_data)
pairs(ENB2012_data)

model_Y1 <- lm(Y1 ~ ., data = ENB2012_data)
model_Y2 <- lm(Y2 ~ ., data = ENB2012_data)

plot(model_Y1)


# =========================================================
# 0. Setup
# =========================================================

library(ranger)

set.seed(1)

data <- ENB2012_data
X <- data[, paste0("X", 1:8)]

n <- nrow(data)

alpha_level <- 0.1   # 90% prediction intervals

# =========================================================
# 1. Funktion: Split Conformal med ranger
# =========================================================

split_conformal_ranger <- function(X, Y, alpha = 0.1) {
  
  n <- nrow(X)
  
  # -------------------------
  # Split data
  # -------------------------
  idx <- sample(1:n, size = floor(n/2))
  train_idx <- idx
  cal_idx   <- setdiff(1:n, train_idx)
  
  X_train <- X[train_idx, ]
  Y_train <- Y[train_idx]
  
  X_cal <- X[cal_idx, ]
  Y_cal <- Y[cal_idx]
  
  # ranger kræver data.frame med target
  train_data <- data.frame(Y = Y_train, X_train)
  cal_data   <- data.frame(X_cal)
  
  # -------------------------
  # Fit Random Forest (ranger)
  # -------------------------
  rf_model <- ranger(
    Y ~ ., 
    data = train_data,
    num.trees = 500,
    min.node.size = 5
  )
  
  # -------------------------
  # Predictions
  # -------------------------
  pred_cal <- predict(rf_model, data = cal_data)$predictions
  
  # -------------------------
  # Nonconformity scores
  # -------------------------
  alpha_i <- abs(Y_cal - pred_cal)
  
  # -------------------------
  # Quantile
  # -------------------------
  q_hat <- quantile(alpha_i, probs = 1 - alpha, type = 1)
  
  # -------------------------
  # Prediction intervals
  # -------------------------
  lower <- pred_cal - q_hat
  upper <- pred_cal + q_hat
  
  # -------------------------
  # Metrics
  # -------------------------
  coverage <- mean(Y_cal >= lower & Y_cal <= upper)
  avg_length <- mean(upper - lower)
  
  # -------------------------
  # Output
  # -------------------------
  list(
    model = rf_model,
    q_hat = q_hat,
    coverage = coverage,
    avg_length = avg_length,
    lower = lower,
    upper = upper,
    pred = pred_cal,
    Y_cal = Y_cal
  )
}

# =========================================================
# 2. Kør for Y1
# =========================================================

result_Y1 <- split_conformal_ranger(X, data$Y1, alpha = alpha_level)

cat("Y1 (Heating Load)\n")
cat("Coverage:", result_Y1$coverage, "\n")
cat("Average length:", result_Y1$avg_length, "\n\n")

# =========================================================
# 3. Kør for Y2
# =========================================================

result_Y2 <- split_conformal_ranger(X, data$Y2, alpha = alpha_level)

cat("Y2 (Cooling Load)\n")
cat("Coverage:", result_Y2$coverage, "\n")
cat("Average length:", result_Y2$avg_length, "\n\n")

# =========================================================
# 4. Plot funktion
# =========================================================

plot_conformal <- function(result, title) {
  
  plot(result$Y_cal, result$pred,
       main = title,
       xlab = "True values",
       ylab = "Predictions",
       pch = 16, col = "blue")
  
  segments(result$Y_cal,
           result$lower,
           result$Y_cal,
           result$upper,
           col = "gray")
  
  abline(0, 1, col = "red", lwd = 2)
}

# =========================================================
# 5. Visualisering
# =========================================================
par(mfrow = c(1, 1))
plot_conformal(result_Y1, "Conformal Prediction - Y1")
plot_conformal(result_Y2, "Conformal Prediction - Y2")











# =========================================================
# 0. Setup
# =========================================================

library(ranger)

set.seed(1)

data <- ENB2012_data
X <- data[, paste0("X", 1:8)]

alpha_level <- 0.1
B <- 1000   # antal gentagelser

# =========================================================
# 1. Funktion: én conformal kørsel
# =========================================================

split_conformal_ranger_once <- function(X, Y, alpha = 0.1) {
  
  n <- nrow(X)
  
  # Split
  idx <- sample(1:n, size = floor(n/2))
  train_idx <- idx
  cal_idx   <- setdiff(1:n, train_idx)
  
  X_train <- X[train_idx, ]
  Y_train <- Y[train_idx]
  
  X_cal <- X[cal_idx, ]
  Y_cal <- Y[cal_idx]
  
  # ranger kræver data.frame
  train_data <- data.frame(Y = Y_train, X_train)
  cal_data   <- data.frame(X_cal)
  
  # Fit model
  model <- ranger(
    Y ~ .,
    data = train_data,
    num.trees = 500,
    min.node.size = 5
  )
  
  # Prediction
  pred_cal <- predict(model, data = cal_data)$predictions
  
  # Nonconformity
  alpha_i <- abs(Y_cal - pred_cal)
  
  # Quantile
  q_hat <- quantile(alpha_i, probs = 1 - alpha, type = 1)
  
  # Intervals
  lower <- pred_cal - q_hat
  upper <- pred_cal + q_hat
  
  # Metrics
  coverage <- mean(Y_cal >= lower & Y_cal <= upper)
  avg_length <- mean(upper - lower)
  
  c(coverage = coverage, avg_length = avg_length)
}

# =========================================================
# 2. Monte Carlo funktion
# =========================================================

run_experiment <- function(X, Y, B = 1000, alpha = 0.1) {
  
  results <- replicate(B, split_conformal_ranger_once(X, Y, alpha))
  
  results <- t(results)
  
  data.frame(
    coverage = results[, "coverage"],
    avg_length = results[, "avg_length"]
  )
}

# =========================================================
# 3. Kør eksperiment for Y1 og Y2
# =========================================================

results_Y1 <- run_experiment(X, data$Y1, B, alpha_level)
results_Y2 <- run_experiment(X, data$Y2, B, alpha_level)

# =========================================================
# 4. Opsummering
# =========================================================

summarize_results <- function(res) {
  c(
    mean_coverage = mean(res$coverage),
    sd_coverage   = sd(res$coverage),
    mean_length   = mean(res$avg_length),
    sd_length     = sd(res$avg_length)
  )
}

summary_Y1 <- summarize_results(results_Y1)
summary_Y2 <- summarize_results(results_Y2)

cat("Y1 summary:\n")
print(summary_Y1)

cat("\nY2 summary:\n")
print(summary_Y2)

# =========================================================
# 5. Visualisering (valgfri men stærk)
# =========================================================

hist(results_Y1$coverage, main = "Coverage Y1", xlab = "Coverage")
hist(results_Y1$avg_length, main = "Interval Length Y1", xlab = "Length")

hist(results_Y2$coverage, main = "Coverage Y2", xlab = "Coverage")
hist(results_Y2$avg_length, main = "Interval Length Y2", 
     xlab = "Length")








running_mean <- function(x) {
  cumsum(x) / seq_along(x)
}
plot_convergence <- function(results, title) {
  
  lengths <- results$avg_length
  
  run_mean <- running_mean(lengths)
  
  plot(run_mean, type = "l",
       main = title,
       xlab = "Iteration (b)",
       ylab = "Running mean of interval length",
       lwd = 2)
  
  abline(h = mean(lengths), lty = 2)  # samlet gennemsnit
  
}

plot_convergence(results_Y1, "Convergence of Interval Length (Y1)")
plot_convergence(results_Y2, "Convergence of Interval Length (Y2)")

plot(results_Y1$avg_length, type = "l",
     main = "Raw Interval Lengths (Y1)",
     xlab = "Iteration",
     ylab = "Length")

running_sd <- function(x) {
  sapply(seq_along(x), function(i) {
    if (i < 2) {
      NA
    } else {
      sd(x[1:i])
    }
  })
}

plot(running_sd(results_Y1$avg_length), type = "l",
     main = "Running SD of Interval Length (Y1)",
     xlab = "Iteration",
     ylab = "SD")







# =============

# QUANTILE REGRESSION 

quantile_ranger_once <- function(X, Y, alpha = 0.1) {
  
  n <- nrow(X)
  
  # Split (samme som før for fair comparison)
  idx <- sample(1:n, size = floor(n/2))
  train_idx <- idx
  test_idx  <- setdiff(1:n, train_idx)
  
  X_train <- X[train_idx, ]
  Y_train <- Y[train_idx]
  
  X_test <- X[test_idx, ]
  Y_test <- Y[test_idx]
  
  train_data <- data.frame(Y = Y_train, X_train)
  
  # Quantile forest
  model <- ranger(
    Y ~ .,
    data = train_data,
    num.trees = 500,
    quantreg = TRUE
  )
  
  # Predict quantiles
  preds <- predict(
    model,
    data = data.frame(X_test),
    type = "quantiles",
    quantiles = c(alpha/2, 1 - alpha/2)
  )$predictions
  
  lower <- preds[, 1]
  upper <- preds[, 2]
  
  # Metrics
  coverage <- mean(Y_test >= lower & Y_test <= upper)
  avg_length <- mean(upper - lower)
  
  c(coverage = coverage, avg_length = avg_length)
}


run_quantile_experiment <- function(X, Y, B = 1000, alpha = 0.1) {
  
  results <- replicate(B, quantile_ranger_once(X, Y, alpha))
  results <- t(results)
  
  data.frame(
    coverage = results[, "coverage"],
    avg_length = results[, "avg_length"]
  )
}


results_qr_Y1 <- run_quantile_experiment(X, data$Y1, B = 1000, alpha = 0.1)
results_qr_Y2 <- run_quantile_experiment(X, data$Y2, B = 1000, alpha = 0.1)

summarize_results <- function(res) {
  c(
    mean_coverage = mean(res$coverage),
    sd_coverage   = sd(res$coverage),
    mean_length   = mean(res$avg_length),
    sd_length     = sd(res$avg_length)
  )
}

cat("Conformal RF (Y1):\n")
print(summarize_results(results_Y1))

cat("\nQuantile RF (Y1):\n")
print(summarize_results(results_qr_Y1))


plot(results_Y1$avg_length, results_qr_Y1$avg_length,
     xlab = "Conformal length",
     ylab = "Quantile length")









