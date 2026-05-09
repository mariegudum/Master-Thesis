# --- 9. Split Conformal Prediction (Elastic Net) ----------------------------
ls()
set.seed(123)
data("riboflavin", package="hdi")

x <- as.matrix(riboflavin[,-1])
y <- as.numeric(riboflavin[,1])

x_scaled <- scale(x)

n <- nrow(x_scaled)

alpha_en <- 0.5

# Train/test split
train_idx <- sample(1:n, size = floor(0.7*n))
test_idx  <- setdiff(1:n, train_idx)

X_train <- x_scaled[train_idx,]
Y_train <- y[train_idx]

X_test <- x_scaled[test_idx,]
Y_test <- y[test_idx]

# Split conformal function using glmnet
split_conformal_glmnet <- function(X, Y, X_new, alpha = 0.05, alpha_en = 0.5){
  
  n <- nrow(X)
  
  # Split train/calibration
  idx <- sample(1:n, floor(n/2))
  
  X_train <- X[idx,]
  Y_train <- Y[idx]
  
  X_cal <- X[-idx,]
  Y_cal <- Y[-idx]
  
  m <- length(Y_cal)
  
  # Elastic net model
  cvfit <- cv.glmnet(X_train, Y_train, alpha = alpha_en)
  model <- glmnet(X_train, Y_train, alpha = alpha_en,
                  lambda = cvfit$lambda.min)
  
  # Calibration residuals
  pred_cal <- predict(model, X_cal)
  R <- abs(Y_cal - pred_cal)
  
  # Quantile
  k <- ceiling((m + 1)*(1 - alpha))
  q_hat <- sort(R)[k]
  
  # Predictions for test data
  pred_test <- as.vector(predict(model, X_new))
  
  lower <- pred_test - q_hat
  upper <- pred_test + q_hat
  
  return(list(lower = lower,
              upper = upper,
              pred = pred_test,
              q_hat = q_hat))
}

# Apply conformal prediction
conf <- split_conformal_glmnet(
  X = X_train,
  Y = Y_train,
  X_new = X_test,
  alpha = 0.05,
  alpha_en = alpha_en
)

covered <- (Y_test >= conf$lower) & (Y_test <= conf$upper)

cat("Marginal coverage:", mean(covered), "\n")
cat("Average interval width:", mean(conf$upper - conf$lower), "\n")

# Use first gene as conditioning variable
X1_test <- X_test[,1]

bins <- cut(
  X1_test,
  breaks = quantile(X1_test, probs = seq(0,1,0.2)),
  include.lowest = TRUE
)

conditional_coverage <- tapply(covered, bins, mean)

print(conditional_coverage)

library(ggplot2)

cond_df <- data.frame(
  bin = names(conditional_coverage),
  coverage = as.numeric(conditional_coverage)
)

cond_df$deviation <- ifelse(cond_df$coverage < 0.95,
                            "Under-coverage",
                            "Over-coverage")

ggplot(cond_df, aes(x = bin, y = coverage, fill = deviation)) +
  geom_col(color="black") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             color = "black") +
  scale_fill_manual(values=c(
    "Under-coverage"="red",
    "Over-coverage"="steelblue"
  )) +
  ylim(0,1) +
  labs(
    title="Conditional Coverage (Riboflavin Data)",
    x="Gene Expression Bin",
    y="Coverage"
  ) +
  theme_minimal()

plot_df <- data.frame(
  X1 = X_test[,1],
  Y = Y_test,
  lower = conf$lower,
  upper = conf$upper,
  covered = covered
)

ggplot(plot_df, aes(x=X1, y=Y)) +
  geom_point(aes(color=covered), alpha=0.7) +
  geom_errorbar(aes(ymin=lower, ymax=upper), alpha=0.3) +
  scale_color_manual(values=c("red","black")) +
  labs(
    title="Conformal Prediction Intervals — Riboflavin",
    x="Gene Expression (Gene 1)",
    y="Riboflavin Production"
  ) +
  theme_minimal()
