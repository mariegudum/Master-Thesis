# SImulate data

set.seed(1)

n <- 1000
X <- runif(n, -2, 2)

sigma <- 0.5 + abs(X)
Y <- X^2 + rnorm(n, 0, sigma)

data <- data.frame(X, Y)

# Split data
train_frac <- 0.6
calib_frac <- 0.2

n_train <- floor(train_frac * n)
n_calib <- floor(calib_frac * n)

idx <- sample(1:n)

train <- data[idx[1:n_train], ]
calib <- data[idx[(n_train+1):(n_train+n_calib)], ]
test  <- data[idx[(n_train+n_calib+1):n], ]

# Split conformal

# Mean model (polynomial to capture nonlinearity)
model <- lm(Y ~ poly(X, 2), data = train)

# Calibration residuals
res_calib <- abs(calib$Y - predict(model, calib))

alpha <- 0.1
q_hat <- quantile(res_calib, 1 - alpha)

# Test predictions
pred_test <- predict(model, test)

lower_cp <- pred_test - q_hat
upper_cp <- pred_test + q_hat

# KNN quantile estimator

knn_quantile <- function(x0, X_train, Y_train, alpha, k = 50) {
  # distances
  d <- abs(X_train - x0)
  
  # indices of k nearest neighbors
  idx <- order(d)[1:k]
  
  # empirical quantile
  quantile(Y_train[idx], probs = alpha)
}

# CQR

alpha_lo <- alpha / 2
alpha_hi <- 1 - alpha / 2

q_lo_cal <- sapply(calib$X, function(x)
  knn_quantile(x, train$X, train$Y, alpha_lo))

q_hi_cal <- sapply(calib$X, function(x)
  knn_quantile(x, train$X, train$Y, alpha_hi))

E_cal <- pmax(q_lo_cal - calib$Y,
              calib$Y - q_hi_cal)

q_cqr <- quantile(E_cal, 1 - alpha)

q_lo_test <- sapply(test$X, function(x)
  knn_quantile(x, train$X, train$Y, alpha_lo))

q_hi_test <- sapply(test$X, function(x)
  knn_quantile(x, train$X, train$Y, alpha_hi))

lower_cqr <- q_lo_test - q_cqr
upper_cqr <- q_hi_test + q_cqr

# Resilts

coverage_cp  <- mean(test$Y >= lower_cp & test$Y <= upper_cp)
coverage_cqr <- mean(test$Y >= lower_cqr & test$Y <= upper_cqr)

width_cp  <- mean(upper_cp - lower_cp)
width_cqr <- mean(upper_cqr - lower_cqr)

cat("CP Coverage:", coverage_cp, "\n")
cat("CQR Coverage:", coverage_cqr, "\n")

cat("CP Width:", width_cp, "\n")
cat("CQR Width:", width_cqr, "\n")

# Repeat

run_experiment <- function(seed) {
  set.seed(seed)
  
  idx <- sample(1:n)
  train <- data[idx[1:n_train], ]
  calib <- data[idx[(n_train+1):(n_train+n_calib)], ]
  test  <- data[idx[(n_train+n_calib+1):n], ]
  
  # CP
  model <- lm(Y ~ poly(X, 2), data = train)
  res_calib <- abs(calib$Y - predict(model, calib))
  q_hat <- quantile(res_calib, 1 - alpha)
  
  pred_test <- predict(model, test)
  lower_cp <- pred_test - q_hat
  upper_cp <- pred_test + q_hat
  
  # CQR
  q_lo_cal <- sapply(calib$X, function(x)
    knn_quantile(x, train$X, train$Y, alpha_lo))
  
  q_hi_cal <- sapply(calib$X, function(x)
    knn_quantile(x, train$X, train$Y, alpha_hi))
  
  E_cal <- pmax(q_lo_cal - calib$Y,
                calib$Y - q_hi_cal)
  
  q_cqr <- quantile(E_cal, 1 - alpha)
  
  q_lo_test <- sapply(test$X, function(x)
    knn_quantile(x, train$X, train$Y, alpha_lo))
  
  q_hi_test <- sapply(test$X, function(x)
    knn_quantile(x, train$X, train$Y, alpha_hi))
  
  lower_cqr <- q_lo_test - q_cqr
  upper_cqr <- q_hi_test + q_cqr
  
  c(
    coverage_cp  = mean(test$Y >= lower_cp & test$Y <= upper_cp),
    coverage_cqr = mean(test$Y >= lower_cqr & test$Y <= upper_cqr),
    width_cp  = mean(upper_cp - lower_cp),
    width_cqr = mean(upper_cqr - lower_cqr)
  )
}

results <- replicate(100, run_experiment(sample(1:10000,1)))
apply(results, 1, mean)

# Plots

x_grid <- seq(-2, 2, length.out = 200)
grid_df <- data.frame(X = x_grid)

# CP model
model <- lm(Y ~ poly(X, 2), data = train)
res_calib <- abs(calib$Y - predict(model, calib))
q_hat <- quantile(res_calib, 1 - alpha)

pred_grid <- predict(model, grid_df)
lower_cp_grid <- pred_grid - q_hat
upper_cp_grid <- pred_grid + q_hat

q_lo_grid <- sapply(x_grid, function(x)
  knn_quantile(x, train$X, train$Y, alpha_lo))

q_hi_grid <- sapply(x_grid, function(x)
  knn_quantile(x, train$X, train$Y, alpha_hi))

lower_cqr_grid <- q_lo_grid - q_cqr
upper_cqr_grid <- q_hi_grid + q_cqr

plot(test$X, test$Y,
     pch = 16, col = rgb(0,0,0,0.3),
     main = "Conformalized Quantile Regression",
     xlab = "X", ylab = "Y")

lines(x_grid, q_lo_grid, col = "darkgreen", lwd = 2)
lines(x_grid, q_hi_grid, col = "darkgreen", lwd = 2)

lines(x_grid, lower_cqr_grid, col = "blue", lwd = 2)
lines(x_grid, upper_cqr_grid, col = "blue", lwd = 2)


plot(test$X, test$Y,
     pch = 16, col = rgb(0,0,0,0.3),
     main = "CP vs CQR",
     xlab = "X", ylab = "Y")

# CP
lines(x_grid, lower_cp_grid, col = "red", lwd = 2)
lines(x_grid, upper_cp_grid, col = "red", lwd = 2)

# CQR
lines(x_grid, lower_cqr_grid, col = "blue", lwd = 2)
lines(x_grid, upper_cqr_grid, col = "blue", lwd = 2)

legend("topright",
       legend = c("CP", "CQR"),
       col = c("red", "blue"),
       lwd = 2)

