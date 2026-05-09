run_experiment <- function(){
  
  n <- nrow(x_scaled)
  
  # 70/30 split
  train_idx <- sample(1:n, floor(0.7*n))
  
  X_train <- x_scaled[train_idx,]
  Y_train <- y[train_idx]
  
  X_test <- x_scaled[-train_idx,]
  Y_test <- y[-train_idx]
  
  # Ridge
  cv_ridge <- cv.glmnet(X_train,Y_train,alpha=0,nfolds=5)
  ridge_model <- glmnet(X_train,Y_train,alpha=0,lambda=cv_ridge$lambda.min)
  
  ridge_pred <- predict(ridge_model,X_test)
  
  ridge_train_pred <- predict(ridge_model,X_train)
  sigma_ridge <- sqrt(mean((Y_train-ridge_train_pred)^2))
  
  z <- qnorm(.975)
  
  ridge_lower <- ridge_pred - z*sigma_ridge
  ridge_upper <- ridge_pred + z*sigma_ridge
  
  ridge_cov <- mean(Y_test>=ridge_lower & Y_test<=ridge_upper)
  ridge_width <- mean(ridge_upper-ridge_lower)
  
  # Lasso
  cv_lasso <- cv.glmnet(X_train,Y_train,alpha=1,nfolds=5)
  lasso_model <- glmnet(X_train,Y_train,alpha=1,lambda=cv_lasso$lambda.min)
  
  lasso_pred <- predict(lasso_model,X_test)
  
  lasso_train_pred <- predict(lasso_model,X_train)
  sigma_lasso <- sqrt(mean((Y_train-lasso_train_pred)^2))
  
  lasso_lower <- lasso_pred - z*sigma_lasso
  lasso_upper <- lasso_pred + z*sigma_lasso
  
  lasso_cov <- mean(Y_test>=lasso_lower & Y_test<=lasso_upper)
  lasso_width <- mean(lasso_upper-lasso_lower)
  
  # Split conformal ridge
  conf_ridge <- split_conformal(X_train,Y_train,X_test,penalty=0)
  
  conf_ridge_cov <- mean(Y_test>=conf_ridge$lower &
                           Y_test<=conf_ridge$upper)
  
  conf_ridge_width <- mean(conf_ridge$upper-conf_ridge$lower)
  
  # Split conformal lasso
  conf_lasso <- split_conformal(X_train,Y_train,X_test,penalty=1)
  
  conf_lasso_cov <- mean(Y_test>=conf_lasso$lower &
                           Y_test<=conf_lasso$upper)
  
  conf_lasso_width <- mean(conf_lasso$upper-conf_lasso$lower)
  
  c(ridge_cov,ridge_width,
    lasso_cov,lasso_width,
    conf_ridge_cov,conf_ridge_width,
    conf_lasso_cov,conf_lasso_width)
}

set.seed(123)

results <- replicate(200, run_experiment())
results <- t(results)


summary_table <- data.frame(
  Method = c("Ridge","Lasso","Conformal Ridge","Conformal Lasso"),
  Coverage = c(mean(results[,1]),
               mean(results[,3]),
               mean(results[,5]),
               mean(results[,7])),
  Width = c(mean(results[,2]),
            mean(results[,4]),
            mean(results[,6]),
            mean(results[,8]))
)

print(summary_table)

# standard deviations
summary_table$Coverage_sd <- c(sd(results[,1]),
                               sd(results[,3]),
                               sd(results[,5]),
                               sd(results[,7]))

summary_table$Width_sd <- c(sd(results[,2]),
                            sd(results[,4]),
                            sd(results[,6]),
                            sd(results[,8]))

# Distribution

coverage_df <- data.frame(
  Ridge = results[,1],
  Lasso = results[,3],
  Conformal_Ridge = results[,5],
  Conformal_Lasso = results[,7]
)

coverage_long <- reshape2::melt(coverage_df)
colnames(coverage_long) <- c("Method","Coverage")


coverage_plot <- ggplot(coverage_long,
                        aes(x=Method,y=Coverage,fill=Method)) +
  geom_boxplot(alpha=0.7) +
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color="red",
             linewidth=1) +
  labs(
    title="Distribution of Empirical Coverage (200 repetitions)",
    x="Method",
    y="Coverage"
  ) +
  theme_minimal() +
  theme(legend.position="none")

coverage_plot


width_df <- data.frame(
  Ridge = results[,2],
  Lasso = results[,4],
  Conformal_Ridge = results[,6],
  Conformal_Lasso = results[,8]
)

width_long <- reshape2::melt(width_df)
colnames(width_long) <- c("Method","Width")


width_plot <- ggplot(width_long,
                     aes(x = Method, y = Width, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Distribution of Prediction Interval Width (200 repetitions)",
    x = "Method",
    y = "Interval Width"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

width_plot

library(patchwork)

combined_plot <- coverage_plot + width_plot
combined_plot


scatter_df <- data.frame(
  Coverage = c(results[,1], results[,3], results[,5], results[,7]),
  Width = c(results[,2], results[,4], results[,6], results[,8]),
  Method = rep(c("Ridge","Lasso","Conformal Ridge","Conformal Lasso"),
               each = nrow(results))
)

scatter_plot <- ggplot(scatter_df,
                       aes(x = Width, y = Coverage, color = Method)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  labs(
    title = "Coverage–Width Tradeoff (200 repetitions)",
    x = "Average Interval Width",
    y = "Empirical Coverage"
  ) +
  theme_minimal()

scatter_plot

