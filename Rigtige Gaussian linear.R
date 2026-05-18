# =========================
# Classical linear prediction interval
# =========================

# Number of estimated parameters (including intercept)
p <- length(coef(lm_fit))

# Residual standard error
rss <- sum((train$y - lm_pred_train)^2)
sigma_hat <- sqrt(rss / (n_train - p))

# t-quantile
t_quant <- qt(1 - alpha/2, df = n_train - p)

# Design matrices
X_train <- model.matrix(lm_fit)
X_test  <- model.matrix(~ x, data = test)

# (X'X)^(-1)
XtX_inv <- solve(t(X_train) %*% X_train)

# Leverage terms for test observations
h <- rowSums((X_test %*% XtX_inv) * X_test)

# Prediction interval widths
pi_width <- t_quant * sigma_hat * sqrt(1 + h)

# Lower and upper bounds
gauss_l <- lm_pred_test - pi_width
gauss_u <- lm_pred_test + pi_width