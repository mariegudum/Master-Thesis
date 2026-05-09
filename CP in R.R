library(ggplot2)
library(insight)
library(patchwork)
library(marginaleffects)

conformal <- function(X, alpha = 0.05) {
  # split the original data randomly in two roughly equal parts
  idx <- sample(c(rep(FALSE, floor(length(X) / 2)), rep(TRUE, ceiling(length(X) / 2))))
  X_train <- X[idx]
  X_calib <- X[!idx]
  # fit a new model in the training set
  X_bar <- mean(X_train)
  # conformity score in calibration set: squared residual
  # d: absolute residual
  conf <- data.frame(
    d = abs(X_calib - X_bar),
    score = (X_calib - X_bar)^2)
  # nth conformity score and d
  conf <- conf[order(conf$score), ]
  idx <- ceiling(((length(X) / 2) + 1) * (1 - alpha))
  d <- conf$d[idx]
  # bounds
  out <- c(X_bar - d, X_bar + d)
  return(out)
}
conformal(rnorm(100, mean = 0, sd = 1))

simulate <- function(mean = 0, sd = 1, N = 1000, alpha = .05) {
  X <- rnorm(N, mean = mean, sd = sd)
  Z <- rnorm(N, mean = mean, sd = sd)
  bounds <- conformal(X, alpha = alpha)
  out <- Z > bounds[1] & Z < bounds[2]
  return(out)
}

# Different DGP parameters
mean(replicate(1000, simulate()))
conformal <- function(model, newdata, alpha = 0.05) {
  # compute bounds using the original data and model
  dat <- insight::get_data(model)
  response <- insight::find_response(model)
  # split the original data randomly in two roughly equal parts
  idx <- sample(c(rep(FALSE, floor(nrow(dat) / 2)), rep(TRUE, ceiling(nrow(dat) / 2))))
  dat_train <- dat[idx, ]
  dat_calib <- dat[!idx, ]
  # fit a new model in training split
  mod <- update(model, data = dat_train)
  # conformity score in calibration split: squared residuals
  # d: absolute residuals
  conf <- predictions(
    mod,
    vcov = FALSE,
    newdata = dat_calib)
  conf$score <- residuals(mod)^2
  conf$d <- abs(residuals(mod))
  # nth conformity score and d
  conf <- conf[order(conf$score), ]
  idx <- ceiling(((nrow(dat) / 2) + 1) * (1 - alpha))
  d <- conf$d[idx]
  # add bounds around predictions made on new data
  out <- predictions(
    mod,
    vcov = FALSE,
    newdata = newdata)
  out$conf.low <- out$estimate - d
  out$conf.high <- out$estimate + d
  return(out)
}
set.seed(1024)

# simulate two datasets
N <- 100

X <- rnorm(N)
past <- data.frame(
  X = X,
  Y = X + X^2 + X^3 + rnorm(N))

X <- rnorm(N)
future <- data.frame(
  X = X,
  Y = X + X^2 + X^3 + rnorm(N))

# fit model
model <- lm(Y ~ X + I(X^2) + I(X^3), data = past)

# predictions on a grid of predictions, using the original model
pred <- predictions(model)

# use conformal inference to build prediction sets
conf <- conformal(model, newdata = future)

# display results
ggplot(pred) +
  geom_ribbon(data = conf, aes(x = X, ymin = conf.low, ymax = conf.high), alpha = .1, fill = "black") +
  geom_ribbon(aes(x = X, ymin = conf.low, ymax = conf.high), fill = "orange") +
  geom_point(aes(X, Y), alpha = .5) +
  geom_line(aes(X, estimate)) +
  xlim(-2, 2)

simulate <- function(misspecified = FALSE, N = 1000) {
  # past: data used to compute bounds 
  X <- rnorm(N)
  past <- data.frame(
    X = X,
    Y = X + X^2 + X^3 + rnorm(N))
  # future: data we want to predict on
  X <- rnorm(N)
  future <- data.frame(
    X = X,
    Y = X + X^2 + X^3 + rnorm(N))
  # misspecified model excludes polynomial terms
  if (isTRUE(misspecified)) {
    model <- lm(Y ~ X + I(X^2) + I(X^3), data = past)
  } else {
    model <- lm(Y ~ X, data = past)
  }
  pred <- conformal(
    model = model,
    newdata = future,
    alpha = .05)
  coverage <- pred$Y > pred$conf.low & pred$Y < pred$conf.high
  return(coverage)
}

# correct model
coverage <- replicate(1000, simulate())
mean(coverage)

# misspecified model
coverage <- replicate(1000, simulate(misspecified = TRUE))
mean(coverage) # 0.9482 
# The coverage rate is about correct, even when the model is misspecified.
