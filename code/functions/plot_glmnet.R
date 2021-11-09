# Description: create a trace plot of lasso or ridge regression, highlighting 
#               certain features
# 
# Arguments:
# glmnet_fit:       fit object returned either by glmnet or cv.glmnet
# lambda:           a value of the penalty parameter; if null and cv.glmnet 
#                    was called, then automatically set to the value chosen by CV
# features_to_plot: either number or names of features to highlight in the plot;
#                    if lambda defined, then for lasso it defaults to plotting 
#                    variables with nonzero coefficients
# 
# Dependencies: The scales and glmnetUtils packages must be installed.
plot_glmnet = function(glmnet_fit, data, lambda = NULL, features_to_plot = NULL){
  
  # extract coefficients
  if(!is.null(glmnet_fit$beta)){
    beta_hat = glmnet_fit$beta
  } else{
    beta_hat = glmnet_fit$glmnet.fit$beta
  }
  
  # extract feature names
  features = rownames(beta_hat)
  
  # extract alpha
  alpha = glmnet_fit$call$alpha
  
  # extract formula
  formula = glmnet_fit$call$formula
  
  # get model matrix
  if(is.null(glmnet_fit$use.model.frame)){
    print("Error: glmnet_fit must be obtained from glmnetUtils rather than glmnet!")
    return()
  } else{
    if(glmnet_fit$use.model.frame){
      X = glmnetUtils:::makeModelComponentsMF(formula, data)$x
    } else{
      X = glmnetUtils:::makeModelComponents(formula, data)$x
    } 
  }
  
  # translate coefficients to standardized scale
  X_centered <- apply(X, 2, function(x) x - mean(x))
  X_sd <- apply(X_centered, 2, function(x) sqrt(sum(x^2) / nrow(X)))
  beta_hat_std = apply(beta_hat, 2, function(x) x*X_sd)
  
  # set lambda using one-standard-error-rule if cv.glmnet() was called
  if(is.null(lambda) & !is.null(glmnet_fit$lambda.1se)){
    lambda = glmnet_fit$lambda.1se
  }
  
  if(!is.null(lambda)){
    lambda = glmnet_fit$lambda[which.min(abs(glmnet_fit$lambda-lambda))]
  }
  
  # set features_to_plot
  if(is.null(features_to_plot)){
    # if not specified but lambda is specified, plot coefficients with nonzero
    # coefficients at lambda = lambda
    if(!is.null(lambda)){
      coefs = beta_hat_std[,glmnet_fit$lambda == lambda]
      if(alpha > 0){
        features_to_plot = features[coefs != 0]
      }
    } 
  } else if(is.numeric(features_to_plot)){
    # if given as integer m
    if(length(features_to_plot) > 1){
      print("Error: if numeric, features_to_plot should have length 1")
      return(NULL)
    } else{
      num_features_to_plot = features_to_plot
      if(alpha > 0){
        # lasso and elastic net:
        #  first num_features_to_plot variables to enter lasso path
        features_to_plot = apply(beta_hat_std, 
                                 1, 
                                 function(row)(min(c(Inf,which(row != 0))))) %>% 
          sort() %>% 
          head(num_features_to_plot) %>% 
          names()
      } else{
        # ridge:
        #  largest num_features_to_plot variables at lambda = lambda
        if(!is.null(lambda)){
          coefs = beta_hat_std[,glmnet_fit$lambda == lambda]
        } else{
          coefs = beta_hat_std[,ncol(beta_hat)]
        }
        features_to_plot = abs(coefs) %>% sort(decreasing = TRUE) %>% 
          head(num_features_to_plot) %>% names()
      }
    }
  } else if(is.character(features_to_plot)){
    # features to plot specified explicitly
  } else{
    print("Error: features_to_plot should either be null, numeric, or a character")
    return(NULL)
  }
  
  # data frame to plot
  df_to_plot = t(beta_hat_std) %>% 
    as.matrix() %>% 
    as_tibble() %>% 
    bind_cols(lambda = glmnet_fit$lambda) %>%
    pivot_longer(-lambda, names_to = "Feature", values_to = "beta_hat_std") %>%
    mutate(feature_label = ifelse(Feature %in% features_to_plot, Feature, ""))
  
  # split features into those that will be highlighted and those that won't
  df_to_plot_highlight = df_to_plot %>% filter(Feature %in% features_to_plot)
  df_to_plot_other = df_to_plot %>% filter(!(Feature %in% features_to_plot))
  
  # define reverse-log transformation
  reverselog_trans <- function(base = 10){
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
              scales::log_breaks(base = base), 
              domain = c(1e-100, Inf))
  }
  
  # produce final plot
  df_to_plot_highlight %>% 
    ggplot(aes(x = lambda, y = beta_hat_std, group = Feature)) + 
    geom_line(data = df_to_plot_other, color = "darkgray") +
    geom_line(aes(color = Feature)) + 
    geom_vline(xintercept = lambda, linetype = "dashed") +
    scale_x_continuous(trans = reverselog_trans()) +  
    xlab(expr(lambda)) +
    ylab("Standardized Coefficient") + 
    theme_bw() + 
    theme(legend.position = "bottom",
          legend.title = element_blank())
}

extract_std_coefs = function(glmnet_fit, data, lambda = NULL){
  
  # extract coefficients
  if(!is.null(glmnet_fit$beta)){
    beta_hat = glmnet_fit$beta
  } else{
    beta_hat = glmnet_fit$glmnet.fit$beta
  }
  
  # extract feature names
  features = rownames(beta_hat)
  
  # extract alpha
  alpha = glmnet_fit$call$alpha
  
  # extract formula
  formula = glmnet_fit$call$formula
  
  # get model matrix
  if(is.null(glmnet_fit$use.model.frame)){
    print("Error: glmnet_fit must be obtained from glmnetUtils rather than glmnet!")
    return()
  } else{
    if(glmnet_fit$use.model.frame){
      X = glmnetUtils:::makeModelComponentsMF(formula, data)$x
    } else{
      X = glmnetUtils:::makeModelComponents(formula, data)$x
    } 
  }
  
  # translate coefficients to standardized scale
  X_centered <- apply(X, 2, function(x) x - mean(x))
  X_sd <- apply(X_centered, 2, function(x) sqrt(sum(x^2) / nrow(X)))
  beta_hat_std = diag(X_sd) %*% beta_hat
  
  # set lambda using one-standard-error-rule if cv.glmnet() was called
  if(is.null(lambda) & !is.null(glmnet_fit$lambda.1se)){
    lambda = glmnet_fit$lambda.1se
  }
  
  tibble(feature = features,
        coefficient = beta_hat_std[,which.min(abs(glmnet_fit$lambda-lambda))])
}

# Description: Plot CV as a function of alpha for cva.glmnet
#               
# 
# Arguments:
# glmnet_fit:       fit object returned by cva.glmnet
# 
# Dependencies: The glmnetUtils package must be installed.
plot_cva_glmnet = function(elnet_fit){
  # extract the list of alpha values used
  alpha = elnet_fit$alpha
  
  # extract the list of minimum CV errors
  error = sapply(elnet_fit$modlist, function(mod) {min(mod$cvm)})
  
  # find the best alpha
  best_alpha =  alpha[which.min(error)]
  
  # plot
  tibble(alpha, error) %>%
    ggplot(aes(x = alpha, y = error)) + 
    geom_point() + 
    geom_line() +
    geom_vline(xintercept = best_alpha, linetype = "dashed") +
    theme_bw() + 
    labs(
      x = expr(alpha),
      y = "Minimum CV error"
    )
}

# Description: Extract best elastic net fit over alpha
#               
# 
# Arguments:
# glmnet_fit:       fit object returned by cva.glmnet
# 
# Dependencies: The glmnetUtils package must be installed.
extract_best_elnet = function(elnet_fit){
  # extract the list of alpha values used
  alpha = elnet_fit$alpha
  
  # extract the list of minimum CV errors
  error = sapply(elnet_fit$modlist, function(mod) {min(mod$cvm)})
  
  # extract the best fit
  out = elnet_fit$modlist[[which.min(error)]]
  
  # append the corresponding value of alpha
  out$alpha = alpha[which.min(error)]
  
  # append use.model.frame and formula
  out$use.model.frame = elnet_fit$use.model.frame
  out$call$formula = elnet_fit$call$formula
  
  # return
  out
}
