install.packages("gbm")
library(mlbench)
library(lpSolve) # to solve quantile regression optimization problem
library(gbm) # nonlinear


# from Meinhausen
data("Ozone", package = "mlbench") # perhaps remove missing data
data("BostonHousing", package = "mlbench")

# Others

data("BreastCancer", package = "mlbench")
data("DNA", package = "mlbench")

head(Ozone)
head(BostonHousing)

df <- BostonHousing

# Response
Y <- df$medv

# Design matrix (include intercept!)
X <- model.matrix(medv ~ ., data = df)

quantile_lp <- function(X, Y, tau) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Objective: (beta, u, v)
  f.obj <- c(rep(0, p), tau * rep(1, n), (1 - tau) * rep(1, n))
  
  # Constraints: X beta + u - v = Y
  A <- cbind(X, diag(n), -diag(n))
  
  result <- lp(direction = "min",
               objective.in = f.obj,
               const.mat = A,
               const.dir = rep("=", n),
               const.rhs = Y)
  
  beta_hat <- result$solution[1:p]
  return(beta_hat)
}

tau <- 0.5

beta_hat <- quantile_lp(X, Y, tau)
beta_hat

Y_pred <- X %*% beta_hat

pinball_loss <- function(y, q, tau) {
  ifelse(y >= q,
         tau * (y - q),
         (1 - tau) * (q - y))
}

mean(pinball_loss(Y, Y_pred, tau))

beta_low  <- quantile_lp(X, Y, tau = 0.025)
beta_high <- quantile_lp(X, Y, tau = 0.975)

Q_low  <- X %*% beta_low
Q_high <- X %*% beta_high

coverage <- mean(Y >= Q_low & Y <= Q_high)
coverage

# 1. Fitted vs Observed (Median)
plot(Y_pred, Y,
     xlab = "Predicted median",
     ylab = "Observed",
     main = "Quantile Regression Fit")
abline(0,1)

# 2. Prediction Intervals Plot
plot(Y, pch=16, col="black")
points(Q_low, col="blue", pch=16)
points(Q_high, col="red", pch=16)
legend("topright", legend=c("Y","Q0.025","Q0.975"),
       col=c("black","blue","red"), pch=16)

# 3. Coverage Visualization
inside <- (Y >= Q_low & Y <= Q_high)

plot(Y, col=ifelse(inside, "green", "red"), pch=16)
legend("topright", legend=c("Covered","Not covered"),
       col=c("green","red"), pch=16)

# 4. Interval Width vs Covariate
width <- Q_high - Q_low

plot(df$lstat, width,
     xlab="lstat",
     ylab="Interval width",
     main="Heteroskedasticity")

# 5. Pinball Loss Distribution
loss <- pinball_loss(Y, Y_pred, 0.5)

hist(loss, breaks=30, main="Pinball Loss Distribution")


# 6. Coefficient Plot 
barplot(beta_hat,
        main="Estimated Quantile Regression Coefficients")



