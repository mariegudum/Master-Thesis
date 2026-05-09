library(glmnet)

full_conformal_glmnet <- function(X, Y, x_new,
                                  alpha = 0.05,
                                  penalty = 1,
                                  grid_size = 100){
  
  n <- nrow(X)
  
  # grid of candidate y values
  y_grid <- seq(min(Y), max(Y), length = grid_size)
  
  accept <- rep(FALSE, grid_size)
  
  for(i in 1:grid_size){
    
    y_candidate <- y_grid[i]
    
    # augmented dataset
    X_aug <- rbind(X, x_new)
    Y_aug <- c(Y, y_candidate)
    
    # fit glmnet path once
    fit <- glmnet(X_aug, Y_aug, alpha = penalty)
    
    # choose lambda using CV on original data
    cv <- cv.glmnet(X, Y, alpha = penalty)
    
    lambda <- cv$lambda.min
    
    # predictions
    preds <- predict(fit, X_aug, s = lambda)
    
    # residuals
    R <- abs(Y_aug - preds)
    
    R_new <- R[n+1]
    
    p_val <- mean(R >= R_new)
    
    accept[i] <- (p_val > alpha)
  }
  
  interval <- range(y_grid[accept])
  
  return(interval)
}

x_new <- X_test[1,,drop=FALSE]

full_ridge <- full_conformal_glmnet(
  X_train,
  Y_train,
  x_new,
  penalty = 0
)

full_lasso <- full_conformal_glmnet(
  X_train,
  Y_train,
  x_new,
  penalty = 1
)


# Coverage and width

ridge_cov <- mean(Y_test >= ridge_lower & Y_test <= ridge_upper)
lasso_cov <- mean(Y_test >= lasso_lower & Y_test <= lasso_upper)

conf_ridge_cov <- mean(Y_test >= full_ridge$lower & Y_test <= full_ridge$upper)
conf_lasso_cov <- mean(Y_test >= full_lasso$lower & Y_test <= full_lasso$upper)

ridge_width <- mean(ridge_upper - ridge_lower)
lasso_width <- mean(lasso_upper - lasso_lower)

conf_ridge_width <- mean(full_ridge$upper - full_ridge$lower)
conf_lasso_width <- mean(full_lasso$upper - full_lasso$lower)

comparison <- data.frame(
  Method=c("Ridge","Lasso","Conformal Ridge","Conformal Lasso"),
  Coverage=c(ridge_cov,lasso_cov,conf_ridge_cov,conf_lasso_cov),
  Width=c(ridge_width,lasso_width,conf_ridge_width,conf_lasso_width)
)
comparison

# Coverage-width tradeoff plot 

coverage_plot <- ggplot(comparison,
                        aes(x=Width,y=Coverage,color=Method))+
  geom_point(size=4)+
  geom_text(aes(label=Method),vjust=-1)+
  geom_hline(yintercept=.95,linetype="dashed",color="red")+
  labs(
    title="Coverage–Width Tradeoff",
    x="Average Interval Width",
    y="Empirical Coverage"
  )+
  theme_minimal()

ggsave("coverage_width_tradeoff.png",
       coverage_plot,
       width=6,
       height=4)

# Conditional coverage 

covered_conf <- (Y_test >= conf_lasso$lower &
                   Y_test <= conf_lasso$upper)

bins <- cut(
  X_test[,1],
  breaks=quantile(X_test[,1],probs=seq(0,1,0.2)),
  include.lowest=TRUE
)

conditional_coverage <- tapply(covered_conf,bins,mean)

cond_df <- data.frame(
  bin=names(conditional_coverage),
  coverage=as.numeric(conditional_coverage)
)

cond_plot <- ggplot(cond_df,
                    aes(x=bin,y=coverage))+
  geom_col(fill="steelblue")+
  geom_hline(yintercept=.95,
             linetype="dashed",
             color="red")+
  ylim(0,1)+
  labs(
    title="Conditional Coverage",
    x="Feature Bin",
    y="Coverage"
  )+
  theme_minimal()

ggsave("conditional_coverage.png",
       cond_plot,
       width=6,
       height=4)
