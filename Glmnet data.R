# ==========================================================
# Split Conformal with glmnet on CoxExample (Fixed)
# ==========================================================

library(glmnet)

set.seed(123)

data(CoxExample, package = "glmnet")

X <- CoxExample$x
n <- nrow(X)
p <- ncol(X)

# ----------------------------------------------------------
# 1. Create sparse high-dimensional beta
# ----------------------------------------------------------

beta_true <- rep(0, p)
beta_true[1:10] <- 2    # sparse signal in first 10 features

# Generate regression outcome
Y <- as.vector(X %*% beta_true + rnorm(n))

# ----------------------------------------------------------
# 2. Train/Test split
# ----------------------------------------------------------

train_idx <- sample(1:n, size = floor(0.7*n))
test_idx  <- setdiff(1:n, train_idx)

X_train <- X[train_idx, ]
Y_train <- Y[train_idx]

X_test  <- X[test_idx, ]
Y_test  <- Y[test_idx]

# ----------------------------------------------------------
# 3. Split Conformal Function
# ----------------------------------------------------------

split_conformal_glmnet <- function(X, Y, X_new, alpha = 0.05){
  
  n <- nrow(X)
  
  # Split into train/calibration
  idx <- sample(1:n, size = floor(n/2))
  
  X_train <- X[idx, ]
  Y_train <- Y[idx]
  
  X_calib <- X[-idx, ]
  Y_calib <- Y[-idx]
  
  m <- length(Y_calib)
  
  # Lasso fit
  cvfit <- cv.glmnet(X_train, Y_train, alpha = 1)
  model <- glmnet(X_train, Y_train, alpha = 1, lambda = cvfit$lambda.min)
  
  # Calibration residuals
  Y_calib_hat <- predict(model, X_calib)
  R <- abs(Y_calib - Y_calib_hat)
  
  # Quantile
  k <- ceiling((m + 1) * (1 - alpha))
  q_hat <- sort(R)[k]
  
  # Predictions on new data
  Y_new_hat <- as.vector(predict(model, X_new))
  
  lower <- Y_new_hat - q_hat
  upper <- Y_new_hat + q_hat
  
  return(list(lower=lower, upper=upper))
}

# ----------------------------------------------------------
# 4. Apply conformal prediction
# ----------------------------------------------------------

alpha <- 0.05

conf <- split_conformal_glmnet(
  X = X_train,
  Y = Y_train,
  X_new = X_test,
  alpha = alpha
)

covered <- (Y_test >= conf$lower) & (Y_test <= conf$upper)

cat("Marginal coverage:", mean(covered), "\n")
cat("Average interval width:", mean(conf$upper - conf$lower), "\n")


# ----------------------------------------------------------
# 6. Approximate Conditional Coverage (Binning)
# ----------------------------------------------------------

# Use first feature for conditioning proxy
X1_test <- X_test[, 1]

# Bin into quintiles
bins <- cut(
  X1_test,
  breaks = quantile(X1_test, probs = seq(0, 1, 0.2)),
  include.lowest = TRUE
)

conditional_coverage <- tapply(covered, bins, mean)

cat("\nConditional coverage by X1 bins:\n")
print(conditional_coverage)


plot_df <- data.frame(
  X1 = X_test[,1],
  Y = Y_test,
  lower = conf$lower,
  upper = conf$upper,
  covered = covered
)

ggplot(plot_df, aes(x = X1, y = Y)) +
  geom_point(aes(color = covered), alpha = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  scale_color_manual(values = c("red", "black")) +
  labs(
    title = "Split Conformal Prediction Intervals",
    subtitle = "Red points indicate miscoverage",
    x = expression(X[1]),
    y = "Y"
  ) +
  theme_minimal()

plot_df$width <- plot_df$upper - plot_df$lower

ggplot(plot_df, aes(x = width)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Distribution of Conformal Interval Widths",
    x = "Interval Width",
    y = "Frequency"
  ) +
  theme_minimal()


library(dplyr)

plot_df$bin <- cut(
  plot_df$X1,
  breaks = quantile(plot_df$X1, probs = seq(0, 1, 0.2)),
  include.lowest = TRUE
)

cond_cov <- plot_df %>%
  group_by(bin) %>%
  summarise(coverage = mean(covered),
            n = n())

ggplot(cond_cov, aes(x = bin, y = coverage)) +
  geom_col(fill = "gray70") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(0,1) +
  labs(
    title = "Approximate Conditional Coverage",
    subtitle = "Dashed line = nominal 95% level",
    x = expression(X[1]~"bins"),
    y = "Coverage"
  ) +
  theme_minimal()

ggplot(cond_cov, aes(x = bin, y = coverage, fill = bin)) +
  geom_col() +
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             color = "red") +
  ylim(0, 1) +
  labs(
    title = "Approximate Conditional Coverage",
    subtitle = "Dashed line = nominal 95% level",
    x = expression(X[1]~"bins"),
    y = "Coverage"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

cond_cov$deviation <- ifelse(cond_cov$coverage < 0.95,
                             "Under-coverage",
                             "Over-coverage")

ggplot(cond_cov, aes(x = bin, y = coverage, fill = deviation)) +
  geom_col() +
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             color = "black") +
  scale_fill_manual(values = c("Under-coverage" = "red",
                               "Over-coverage" = "steelblue")) +
  ylim(0, 1) +
  labs(
    title = "Approximate Conditional Coverage",
    subtitle = "Red bars fall below nominal 95% level",
    x = expression(X[1]~"bins"),
    y = "Coverage"
  ) +
  theme_minimal()


ggplot(plot_df, aes(x = X1, y = as.numeric(covered))) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(
    title = "Smoothed Conditional Coverage",
    x = expression(X[1]),
    y = "Estimated Conditional Coverage"
  ) +
  theme_minimal()

plot_df$pred <- (plot_df$upper + plot_df$lower)/2

ggplot(plot_df, aes(x = pred, y = as.numeric(covered))) +
  geom_smooth(method = "loess", se = FALSE, color = "darkgreen") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(
    title = "Coverage vs Predicted Value",
    x = "Predicted Y",
    y = "Estimated Coverage"
  ) +
  theme_minimal()

ggplot(plot_df, aes(x = X1, y = as.numeric(!covered))) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(
    title = "Smoothed Miscoverage Probability",
    y = "Estimated P(Y not in C(X))"
  ) +
  theme_minimal()

ggplot(plot_df, aes(x = X1, y = pred)) +
  geom_point(alpha = .6) +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = .3) +
  theme_minimal()
