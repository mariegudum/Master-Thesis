set.seed(42)

# -----------------------------
# Step 1: Simulate data
# -----------------------------
n_calib <- 100
n_test  <- 50

X_calib <- runif(n_calib, 0, 1)
X_test  <- seq(0, 1, length.out = n_test)

Y_calib <- sin(2 * pi * X_calib) + 0.1 * rnorm(n_calib)
Y_test  <- sin(2 * pi * X_test)  + 0.1 * rnorm(n_test)

# -----------------------------
# Step 2: Fit linear regression
# -----------------------------
model <- lm(Y_calib ~ X_calib)

Y_pred_calib <- predict(model, newdata = data.frame(X_calib = X_calib))
Y_pred_test  <- predict(model, newdata = data.frame(X_calib = X_test))

# -----------------------------
# Step 3: Nonconformity scores
# -----------------------------
scores <- abs(Y_calib - Y_pred_calib)

# -----------------------------
# Step 4: Stratification
# -----------------------------
n_strata <- 5
strata_edges <- seq(0, 1, length.out = n_strata + 1)

strata_scores <- vector("list", n_strata)

for (i in 1:n_strata) {
  mask <- X_calib >= strata_edges[i] & X_calib < strata_edges[i+1]
  strata_scores[[i]] <- scores[mask]
}

# -----------------------------
# Step 5: Compute intervals
# -----------------------------
alpha <- 0.1
intervals_lower <- numeric(n_test)
intervals_upper <- numeric(n_test)

for (j in 1:n_test) {
  x <- X_test[j]
  
  # Find stratum
  stratum_index <- findInterval(x, strata_edges, rightmost.closed = TRUE)
  stratum_index <- min(stratum_index, n_strata)
  
  calib_scores <- strata_scores[[stratum_index]]
  
  if (length(calib_scores) == 0) {
    q <- quantile(scores, 1 - alpha)
  } else {
    k <- ceiling((length(calib_scores) + 1) * (1 - alpha))
    q <- sort(calib_scores)[min(k, length(calib_scores))]
  }
  
  y_hat <- Y_pred_test[j]
  intervals_lower[j] <- y_hat - q
  intervals_upper[j] <- y_hat + q
}

# -----------------------------
# Step 6: Conditional coverage
# -----------------------------
for (i in 1:n_strata) {
  mask <- X_test >= strata_edges[i] & X_test < strata_edges[i+1]
  
  if (sum(mask) == 0) next
  
  covered <- Y_test[mask] >= intervals_lower[mask] &
    Y_test[mask] <= intervals_upper[mask]
  
  coverage <- mean(covered)
  
  cat(sprintf("Stratum %d: coverage = %.2f, target = %.2f, #calib pts = %d\n",
              i, coverage, 1 - alpha, length(strata_scores[[i]])))
}

# -----------------------------
# Plot
# -----------------------------
plot(X_test, Y_test, col="red", pch=16,
     ylim=range(c(intervals_lower, intervals_upper, Y_test)),
     main="Stratified Conformal Prediction",
     xlab="X", ylab="Y")

points(X_calib, Y_calib, col="blue", pch=16)
lines(X_test, Y_pred_test, col="black", lwd=2)
lines(X_test, intervals_lower, col="orange", lty=2)
lines(X_test, intervals_upper, col="orange", lty=2)
