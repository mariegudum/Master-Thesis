library(ciTools)
library(trending)
library(ggplot2)
library(patchwork)
library(MASS)

# generate data
x <- rnorm(100, mean = 0)
y <- rpois(n = 100, lambda = exp(1.5 + 0.5*x))
dat <- data.frame(x = x, y = y)
fit <- glm(y ~ x , family = poisson(link = "log"))

# use ciTools to add prediction interval
dat1 <- add_pi(dat, fit, names = c("lpb", "upb"), alpha = 0.1, nsims = 20000)
#> Warning in add_pi.glm(dat, fit, names = c("lpb", "upb"), alpha = 0.1, nsims =
#> 20000): The response is not continuous, so Prediction Intervals are approximate
head(dat1)
#>            x  y     pred lpb upb
#> 1 -1.7051837  2 1.978464   0   5
#> 2 -0.7533523  3 3.232369   1   6
#> 3  0.4280385  9 5.944711   2  10
#> 4 -0.1654847  2 4.377160   1   8
#> 5  0.4824726  4 6.113966   2  11
#> 6  0.9575575 12 7.811480   4  13

# add intervals with trending (no uncertainty in parameters)
poisson_model <- glm_model(y ~ x, family = "poisson")
fitted_model <- fit(poisson_model, dat)
dat2 <- predict(fitted_model, simulate_pi = FALSE, uncertain = FALSE, alpha = 0.1)
dat2 <- get_result(dat2)
head(dat2[[1]])
#> <trending_prediction> 6 x 7
#>       y      x estimate lower_ci upper_ci lower_pi upper_pi
#>   <int>  <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
#> 1     2 -1.71      1.98     1.66     2.36        0        5
#> 2     3 -0.753     3.23     2.88     3.63        1        6
#> 3     9  0.428     5.94     5.54     6.38        2       10
#> 4     2 -0.165     4.38     4.01     4.77        1        8
#> 5     4  0.482     6.11     5.70     6.56        2       10
#> 6    12  0.958     7.81     7.24     8.43        4       13

# add intervals with trending (uncertainty in parameters)
dat3 <- predict(fitted_model, simulate_pi = FALSE, alpha = 0.1)
dat3 <- get_result(dat3)
head(dat3[[1]])
#> <trending_prediction> 6 x 7
#>       y      x estimate lower_ci upper_ci lower_pi upper_pi
#>   <int>  <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
#> 1     2 -1.71      1.98     1.66     2.36        0        5
#> 2     3 -0.753     3.23     2.88     3.63        0        7
#> 3     9  0.428     5.94     5.54     6.38        2       11
#> 4     2 -0.165     4.38     4.01     4.77        1        9
#> 5     4  0.482     6.11     5.70     6.56        2       11
#> 6    12  0.958     7.81     7.24     8.43        3       13

# plots
p1 <- ggplot(dat1, aes(x, y)) +
  geom_point(size = 1) +
  geom_line(aes(y = pred), size = 1.2) +
  geom_ribbon(aes(ymin = lpb, ymax = upb), alpha = 0.2) +
  geom_ribbon(aes(ymin = `lower_pi`, ymax = `upper_pi`), data = dat2[[1]], alpha = 0.4) +
  ggtitle(
    "Poisson regression with prediction intervals and no uncertainty in parameters", 
    subtitle = "Model fit (black line), with bootstrap intervals (gray), parametric intervals (dark gray)"
  ) +
  coord_cartesian(ylim=c(0, 30))
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.

p2 <- ggplot(dat1, aes(x, y)) +
  geom_point(size = 1) +
  geom_line(aes(y = pred), size = 1.2) +
  geom_ribbon(aes(ymin = lpb, ymax = upb), alpha = 0.4) +
  geom_ribbon(aes(ymin = `lower_pi`, ymax = `upper_pi`), data = dat3[[1]], alpha = 0.2) +
  ggtitle(
    "Poisson regression with prediction intervals and uncertainty in parameters", 
    subtitle = "Model fit (black line), with parametric intervals (gray), bootstrap intervals (dark gray)"
  ) +
  coord_cartesian(ylim=c(0, 30))

p1 / p2

