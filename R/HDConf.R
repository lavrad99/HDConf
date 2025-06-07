pinv_SVD <- function(X , numerical_zero = 1e-7){
  #' Function for performing the Moore-Penrose pseudo-inverse using SVD
  #' @param X is a matrix
  #' @param numerical_zero is the level at which an eigenvalue is considered zero
  #' @return The Moore-Penrose pseudo-inverse of X
  my_helpf <- function(z){
    if(abs(z) < numerical_zero){
      return(0)
    }else{
      return(1/z)
    }
  }
  my_SVD <- propack.svd(X)
  return(my_SVD$v%*%diag(sapply(my_SVD$d , my_helpf))%*%t(my_SVD$u))
}


Ridge_SVD <- function(X,y,lambda_list){
  #' Function for Ridge regression using SVD
  #' @param X is the design matrix
  #' @param y is the response vector
  #' @param lambda_list is the list of lambdas for which we want to perform Ridge regression
  #' @return the ridge coefficients when regressing y on X, using lambda_list as penalties

  my_SVD <- propack.svd(X)
  U_Ty <- t(my_SVD$u)%*%y
  d <- my_SVD$d
  d_sq <- d^2
  V <- my_SVD$v
  my_coefs <- sapply(lambda_list , function(lambda){return(V%*%diag(d/(lambda+d_sq))%*%U_Ty)})
  return(my_coefs)
}

create_lambda_list_Ridge <- function(X , n_lambdas = 200){
  #' Function for creating a set of candidate lambdas for ridge regression
  #' @param X is the design matrix
  #' @param n_lambdas is the number of lambdas we want
  #' @return vector of candidate lambdas

  sq_F_norm <- norm(X , type = 'F')^2
  min_log_value <- -10
  return(sq_F_norm*exp(min_log_value+(-min_log_value)*(1:n_lambdas)/(n_lambdas+1)))

}


CV_select_lambda_Ridge <- function(X,y,n_lambdas = 200 , n_folds = 5){
  #' Function for K-fold cross-validation using Ridge regression
  #' @param X is the design matrix
  #' @param y is the response vector
  #' @param n_lambdas is the bumber of lambdas we want to test
  #' @param n_folds is the number of folds
  #' @return the average MSE across folds for each lambda and the lambda with the lowest average MSE


  my_ret <- list()

  n <- dim(X)[1]

  my_rearangement <- sample(x = 1:n , size = n , replace = F)

  X <- X[my_rearangement , ] #Permute columns to destroy any row-wise structure
  y <- y[my_rearangement]

  my_lambdas <- create_lambda_list_Ridge(X = X , n_lambdas = n_lambdas)
  my_ret$lambdas <- my_lambdas
  my_folds <- split(1:n, cut(1:n, breaks = n_folds, labels = FALSE))

  MSE_fold <- function(fold){
    X_train <- X[-fold ,]
    y_train <- y[-fold]

    X_test <- X[fold ,]
    y_test <- y[fold]

    my_coefs <- Ridge_SVD(X = X_train , y = y_train , lambda_list = my_lambdas)

    my_MSEs <- apply(X = X_test%*%my_coefs-y_test ,MARGIN = 2, FUN = function(x){return(mean(x^2))})

    return(my_MSEs)
  }

  my_MSEs_all_folds <- sapply(my_folds, MSE_fold)


  avg_MSE <- apply(X = my_MSEs_all_folds , MARGIN = 1 , mean)

  my_ret$avg_MSEs <- avg_MSE

  my_ret$lambda.min <- my_lambdas[which.min(avg_MSE)]

  return(my_ret)
}



HD_OLS <- function(x,y,conf, min_DF = 20){
  #' Function for high-dimensional OLS
  #' @param x is the exposure vector
  #' @param y is the response vector
  #' @param conf is the confound matrix
  #' @param min_DF is the minimal number of degrees of freedom in the PC-OLS approach
  #' @return the p-value and the sign of the coefficient for x
  #' @details If the full model has degrees of freedom above min_DF, the function defaults to regular OLS

  my_ret <- list()

  if(is.null(dim(conf)[2])){

    warning("Confound matrix has no second dimension. Returned p-value without deconfounding.")

    my_test <- cor.test(x,y)

    my_ret$sign_beta_x <- sign(my_test[['estimate']])

    my_ret$pval <- my_test[['p.value']]

    return(my_ret)

  }

  if(dim(conf)[1]-dim(conf)[2]-2 >= min_DF){
    my_df <- as.data.frame(unname(conf))
    my_df$x <- x
    my_df$y <- y

  }else{
    n <- dim(conf)[1]
    n_PCs <- n-min_DF-2
    conf_PC <- prcomp(x = conf , rank. = n_PCs, scale. = T , center = T)$x
    my_df <- as.data.frame(conf_PC)
    my_df$x <- x
    my_df$y <- y
  }

  my_HD_OLS_summary <- summary(lm(formula = y ~ . , data = my_df))

  my_ret$sign_beta_x <- sign(my_HD_OLS_summary$coefficients['x' , 'Estimate'])

  my_ret$pval <- my_HD_OLS_summary$coefficients['x' , 'Pr(>|t|)']

  return(my_ret)

}


PCA_OLS <- function(x,y,conf,prop_VE = 0.95, min_DF = 20){
  #' Function for high-dimensional "OLS-PCA"
  #' @param x is the exposure vector
  #' @param y is the response vector
  #' @param conf is the confound matrix
  #' @param prop_VE is the % of Variance Explained criterion for including PCs
  #' @param min_DF is the minimal number of degrees of freedom in the OLS-PCA approach
  #' @return the p-value and the sign of the coefficient for x
  #' @details If the VE criterion yields degrees of freedom below min_DF, the function defaults to HD_OLS

  my_ret <- list()

  my_PCA <- prcomp(x = conf , scale. = T, center = T)

  n_PCs_VE <- which.max(cumsum(summary(my_PCA)$importance[2,]) > prop_VE)

  n_PCs_min_DF <- dim(conf)[1] - min_DF - 2

  n_PCs <- min(c(n_PCs_VE , n_PCs_min_DF))

  my_ret$n_PCs <- n_PCs

  my_df <- as.data.frame(my_PCA$x[ , 1:n_PCs])

  my_df$x <- x

  my_df$y <- y

  my_PCA_OLS_summary <- summary(lm(y ~ . , data = my_df))

  my_ret$sign_beta_x <- sign(my_PCA_OLS_summary$coefficients['x' , 'Estimate'])

  my_ret$pval <- my_PCA_OLS_summary$coefficients['x' , 'Pr(>|t|)']

  return(my_ret)

}



Ridge_pvalue <- function(x,y,conf, method = 'Cule-DeIorio', n_lambdas = 200, t_test = T){
  #' Function for extracting p-values from Ridge regression using the method described in Cule, Vineis and De Iorio (2011)
  #' @param x is the exposure vector
  #' @param y is the response vector
  #' @param conf is the confound matrix
  #' @param method refers to method for selecting ridge penalty ("Cule-DeIorio" (default) or "CV")
  #' @param n_lambdas refers is the number of lambdas used in cross-validation (only used when method = "CV")
  #' @param t_test is a boolean determining whether to use t-tests (default) or z-tests
  #' @return the p-value and the sign of the coefficient for x
  my_ret <- list()
  if(method == 'Cule-DeIorio'){
    my_df <- as.data.frame(unname(conf))

    my_df$x <- x

    my_df$y <- y

    my_Ridge <- linearRidge(y ~ . , data = my_df , lambda = 'automatic')

    scaled_X <- my_Ridge$x

    y <- my_Ridge$y

    lambda_Ridge <- my_Ridge$lambda[my_Ridge$chosen.nPCs]

    XTX <- t(scaled_X)%*%scaled_X

    XTX_plus_ldiag_inv <- solve(XTX+lambda_Ridge*diag(dim(my_df)[2]-1))

    H <- scaled_X%*%XTX_plus_ldiag_inv%*%t(scaled_X)

    df_noise <- dim(my_df)[1]-2*sum(diag(H))+sum(H^2)

    df_tstat <- dim(my_df)[1]-sum(diag(H))

    y_hat <- H%*%y

    sigma_sq_hat <- t(y-y_hat)%*%(y-y_hat)/df_noise

    var_beta <- as.numeric(sigma_sq_hat)*XTX_plus_ldiag_inv%*%XTX%*%XTX_plus_ldiag_inv

    SE_x <- sqrt(var_beta['x' , 'x'])

    beta_hat_x <- my_Ridge$coef['x' , my_Ridge$chosen.nPCs]

    my_ret$sign_beta_x <- sign(beta_hat_x)

    z_score_x <- beta_hat_x/SE_x
    if(t_test == F){
      my_ret$pval <- 2*min(pnorm(q = z_score_x) , 1-pnorm(q = z_score_x))
    }else{
      my_ret$pval <- 2*min(pt(q = z_score_x , df = df_tstat) , 1-pt(q = z_score_x , df = df_tstat))
    }

    return(my_ret)
  }else if(method == 'CV'){
    my_df <- as.data.frame(unname(conf))

    my_df$x <- x

    my_df$y <- y

    X <- scale(as.matrix(base::subset(my_df , select = -c(y))))

    y <- scale(as.matrix(base::subset(my_df , select = c(y))))

    my_Ridge_CV <- CV_select_lambda_Ridge(X=X,y=y , n_lambdas = n_lambdas)

    lambda_Ridge <- my_Ridge_CV$lambda.min

    XTX <- t(X)%*%X

    XTX_plus_ldiag_inv <- solve(XTX+lambda_Ridge*diag(dim(my_df)[2]-1))

    H <- X%*%XTX_plus_ldiag_inv%*%t(X)

    df_noise <- dim(my_df)[1]-sum(diag(2*H))+sum(H^2)

    df_tstat <- dim(my_df)[1]-sum(diag(H))

    y_hat <- H%*%y

    sigma_sq_hat <- t(y-y_hat)%*%(y-y_hat)/df_noise

    var_beta <- as.numeric(sigma_sq_hat)*XTX_plus_ldiag_inv%*%XTX%*%XTX_plus_ldiag_inv

    SE_x <- sqrt(var_beta['x' , 'x'])

    beta_hat_x <- Ridge_SVD(X=X,y=y, lambda_list = c(lambda_Ridge))[[dim(X)[2]]]

    my_ret$sign_beta_x <- sign(beta_hat_x)

    z_score_x <- beta_hat_x/SE_x

    if(t_test == F){
      my_ret$pval <- 2*min(pnorm(q = z_score_x) , 1-pnorm(q = z_score_x))
    }else{
      my_ret$pval <- 2*min(pt(q = z_score_x , df = df_tstat) , 1-pt(q = z_score_x , df = df_tstat))
    }
    return(my_ret)
  }
  stop("invalid method name")
}

Ridge_PT <- function(x,y,conf,test_method = 'FLH1', fit_method = 'svd' ,n_perms = 2000 , n_lambdas = 200){
  #' Function performing permutation testing using ridge regression as described in Hemerik, Thoresen and Finos 2020
  #' @param x is the exposure vector
  #' @param y is the response vector
  #' @param conf is the confound matrix
  #' @param test_method is the method used for testing,"FLH1" (default), "FLH2" or "DR"
  #' @param fit_method is the method used for fitting, "svd" (default) or "glmnet"
  #' @param n_perms is the number of permutations in the test
  #' @param n_lambdas is the number of lambdas tested for selection using 5-fold cross-validation
  #' @return the p-value and the sign of the correlation between x and y after adjusting for confounds

  my_ret <- list()

  Z <- scale(conf)

  x <- scale(x)

  y <- scale(y)

  if(fit_method == 'svd'){
    if(test_method == 'FLH1'){

      my_Ridge_CV_x <- CV_select_lambda_Ridge(X=Z,y=x,n_lambdas = n_lambdas)

      lambda_x <- my_Ridge_CV_x$lambda.min

      my_Ridge_CV_y <- CV_select_lambda_Ridge(X=Z,y=y,n_lambdas = n_lambdas)

      lambda_y <- my_Ridge_CV_y$lambda.min

      beta_maker_x <- solve(t(Z)%*%Z+lambda_x*diag(dim(Z)[2]))%*%t(Z)

      beta_maker_y <- solve(t(Z)%*%Z+lambda_y*diag(dim(Z)[2]))%*%t(Z)

      my_signal <- Z%*%beta_maker_y%*%y

      my_resids <- y - my_signal

      my_resids_x <- x - Z%*%beta_maker_x%*%x

      my_ret$sign_beta_x <- sign(cor(my_resids , my_resids_x)[[1]])

      T_1 <- abs(cor(my_resids , my_resids_x)[[1]])

      my_new_ys <- replicate(n = n_perms ,
                             expr = matrix(my_signal+sample(my_resids, replace=FALSE), ncol =1),
                             simplify = 'matrix')

      my_new_resids <- my_new_ys - Z%*%beta_maker_y%*%my_new_ys

      T_list <- c(T_1 , apply(X = my_new_resids, MARGIN = 2, FUN = function(a){return(abs(cor(a,my_resids_x)[[1]]))}))

      my_ret$pval <-  mean(T_list >= T_1)

    }else if(test_method == 'FLH2'){

      my_Ridge_CV_y <- CV_select_lambda_Ridge(X=Z,y=y,n_lambdas = n_lambdas)

      lambda_y <- my_Ridge_CV_y$lambda.min

      beta_maker_y <- solve(t(Z)%*%Z+lambda_y*diag(dim(Z)[2]))%*%t(Z)

      my_signal <- Z%*%beta_maker_y%*%y

      my_resids <- y - my_signal

      my_ret$sign_beta_x <- sign(cor(my_resids , x)[[1]])

      T_1 <- abs(cor(my_resids , x)[[1]])

      my_new_ys <- replicate(n = n_perms ,
                             expr = matrix(my_signal+sample(my_resids, replace=FALSE), ncol =1),
                             simplify = 'matrix')

      my_new_resids <- my_new_ys - Z%*%beta_maker_y%*%my_new_ys

      T_list <- c(T_1 , apply(X = my_new_resids, MARGIN = 2, FUN = function(a){return(abs(cor(a,x)[[1]]))}))

      my_ret$pval <- mean(T_list >= T_1)

    }else if(test_method == 'DR'){

      my_Ridge_CV_x <- CV_select_lambda_Ridge(X=Z,y=x,n_lambdas = n_lambdas)

      lambda_x <- my_Ridge_CV_x$lambda.min

      my_Ridge_CV_y <- CV_select_lambda_Ridge(X=Z,y=x,n_lambdas = n_lambdas)

      lambda_y <- my_Ridge_CV_y$lambda.min

      beta_maker_x <- solve(t(Z)%*%Z+lambda_x*diag(dim(Z)[2]))%*%t(Z)

      beta_maker_y <- solve(t(Z)%*%Z+lambda_y*diag(dim(Z)[2]))%*%t(Z)

      my_signal <- Z%*%beta_maker_y%*%y

      my_resids <- y - my_signal

      my_resids_x <- x - Z%*%beta_maker_x%*%x

      my_ret$sign_beta_x <- sign(cor(y , my_resids_x)[[1]])

      T_1 <- abs(cor(y , my_resids_x)[[1]])


      my_new_ys <- replicate(n = n_perms ,
                             expr = matrix(my_signal+sample(my_resids, replace=FALSE), ncol =1),
                             simplify = 'matrix')

      T_list <- c(T_1 , apply(X = my_new_ys, MARGIN = 2, FUN = function(a){return(abs(cor(a,my_resids_x)[[1]]))}))

      my_ret$pval <- mean(T_list >= T_1)

    }else{
      stop("invalid option for test_method")
    }

  }else if(fit_method == 'glmnet'){

    if(test_method == 'FLH1'){

      my_Ridge_CV_x <- cv.glmnet(x = Z , y = x, alpha = 0, intercept = F,
                                 nlambda = n_lambdas, nfolds = 5)

      lambda_x <- my_Ridge_CV_x$lambda.min

      my_Ridge_CV_y <- cv.glmnet(x = Z , y = y, alpha = 0, intercept = F,
                                 nlambda = n_lambdas, nfolds = 5)

      lambda_y <- my_Ridge_CV_y$lambda.min

      my_Ridge_x <- glmnet(x = Z ,
                           y = x, alpha = 0, intercept = F, lambda = lambda_x)

      my_Ridge_y <- glmnet(x = Z ,
                           y = y, alpha = 0, intercept = F, lambda = lambda_y)

      my_signal <- predict(my_Ridge_y , Z)

      my_resids <- y - my_signal

      my_resids_x <- x - predict(my_Ridge_x , Z)

      T_1 <- abs(cor(my_resids , my_resids_x)[[1]])

      my_ret$sign_beta_x <- sign(cor(my_resids , my_resids_x)[[1]])

      my_ret$sign_beta_x <- sign(cor(my_resids , my_resids_x)[[1]])

      my_null_draw <- function(){
        y_new <- my_signal+sample(my_resids, replace=FALSE)
        my_Ridge_new <- glmnet(x = Z , y = y_new,
                               alpha = 0, intercept = F, lambda = lambda_y)

        return(abs(cor(y_new - predict(my_Ridge_new , Z) , my_resids_x)[[1]]))
      }


      T_list <- c(T_1 , replicate( n = n_perms , expr = my_null_draw()))

      my_ret$pval <- mean(T_list >= T_1)

    }else if(test_method == 'FLH2'){

      my_Ridge_CV_y <- cv.glmnet(x = Z , y = y, alpha = 0, intercept = F,
                                 nlambda = n_lambdas, nfolds = 5)

      lambda_y <- my_Ridge_CV_y$lambda.min

      my_Ridge_y <- glmnet(x = Z ,
                           y = y, alpha = 0, intercept = F, lambda = lambda_y)

      my_signal <- predict(my_Ridge_y , Z)

      my_resids <- y - my_signal

      T_1 <- abs(cor(my_resids , x)[[1]])

      my_null_draw <- function(){
        y_new <- my_signal+sample(my_resids, replace=FALSE)
        my_Ridge_new <- glmnet(x = Z , y = y_new,
                               alpha = 0, intercept = F, lambda = lambda_y)

        return(abs(cor(y_new - predict(my_Ridge_new , Z) , x)[[1]]))
      }

      T_list <- c(T_1 , replicate( n = n_perms , expr = my_null_draw()))

      pvalue_Ridge_PT <- mean(T_list >= T_1)

    }else if(test_method == 'DR'){

      my_Ridge_CV_x <- cv.glmnet(x = Z , y = x, alpha = 0, intercept = F,
                                 nlambda = n_lambdas, nfolds = 5)

      lambda_x <- my_Ridge_CV_x$lambda.min

      my_Ridge_CV_y <- cv.glmnet(x = Z , y = y, alpha = 0, intercept = F,
                                 nlambda = n_lambdas, nfolds = 5)

      lambda_y <- my_Ridge_CV_y$lambda.min

      beta_maker_x <- solve(t(Z)%*%Z+dim(Z)[1]*lambda_x*diag(dim(Z)[2]))%*%t(Z)

      beta_maker_y <- solve(t(Z)%*%Z+dim(Z)[1]*lambda_y*diag(dim(Z)[2]))%*%t(Z)

      my_signal <- Z%*%beta_maker_y%*%y

      my_resids <- y - my_signal

      my_resids_x <- x - Z%*%beta_maker_y%*%x

      my_ret$sign_beta_x <- sign(cor(y , my_resids_x)[[1]])

      T_1 <- abs(cor(y , my_resids_x)[[1]])


      my_new_ys <- replicate(n = n_perms ,
                             expr = matrix(my_signal+sample(my_resids, replace=FALSE), ncol =1),
                             simplify = 'matrix')

      T_list <- c(T_1 , apply(X = my_new_ys, MARGIN = 2, FUN = function(a){return(abs(cor(a,my_resids_x)[[1]]))}))

      my_ret$pval <- mean(T_list >= T_1)

    }else{
      stop("invalid option for test_method")
    }



  }else{
    stop("invalid option for fit_method")
  }

  return(my_ret)

}



Naive_LASSO <- function(x,y,conf,criterion = 'min',n_lambdas = 200, min_DF = 20){
  #' Function performing high-dimensional estimation via Naive LASSO selection (see Zhao, Witten and Shojaie (2021))
  #' @param x is the exposure vector
  #' @param y is the response vector
  #' @param conf is the confound matrix
  #' @param n_lambdas is the number of lambdas tested for selection using 5-fold cross-validation
  #' @param criterion refers to the criterion for choosing lambda, "min" (default) or "1se"
  #' @param min_DF is the minimal number of degrees of freedom in the final OLS model
  #' @return the p-value and the sign of the coefficient for x
  #' @details If the full model has degrees of freedom above min_DF, the function uses HD_OLS at the last step
  my_ret <- list()

  C <- scale(as.matrix(unname(conf)))

  x <- scale(as.matrix(x))

  y <- scale(as.matrix(y))

  n <- dim(C)[1]

  my_LASSO_CV <- cv.glmnet(x = C, y = y, nlambda = n_lambdas, alpha = 1, intercept = F, nfolds = 5)

  if(criterion == 'min'){
    my_lambda <- my_LASSO_CV$lambda.min
  }else if(criterion == '1se'){
    my_lambda <- my_LASSO_CV$lambda.1se
  }else{
    stop("invalid name of criterion")
  }
  my_LASSO <- glmnet(x = C , y = y,
                     family = 'gaussian', alpha = 1,
                     lambda = my_lambda,
                     intercept = F)

  my_coefs_indic <- coef(my_LASSO)[-1] != 0

  my_ret$selected_vars <- my_coefs_indic

  my_HD_OLS<- HD_OLS(x = x, y = y,
                     conf = C[ , my_coefs_indic],
                     min_DF = min_DF)

  my_ret$sign_beta_x <- my_HD_OLS$sign_beta_x

  my_ret$pval <- my_HD_OLS$pval

  return(my_ret)

}



Naive_ElasticNet <- function(x,y,conf,criterion = 'min',n_lambdas = 200, n_alphas = 10, min_DF = 20){
  #' Function performing high-dimensional estimation via Naive ElasticNet selection
  #' @param x is the exposure vector
  #' @param y is the response vector
  #' @param conf is the confound matrix
  #' @param n_lambdas is the number of lambdas tested for selection using 5-fold cross-validation
  #' @param n_alphas is the number of alphas tested for selection using 5-fold cross-validation
  #' @param criterion refers to the criterion for choosing lambda, "min" (default) or "1se"
  #' @param min_DF is the minimal number of degrees of freedom in the final OLS model
  #' @return the p-value and the sign of the coefficient for x
  #' @details If the full model has degrees of freedom above min_DF, the function uses HD_OLS at the last step

  my_ret <- list()

  C <- scale(as.matrix(unname(conf)))

  x <- scale(as.matrix(x))

  y <- scale(as.matrix(y))

  n <- dim(C)[1]

  alpha_list <- (1:n_alphas)/(1+n_alphas)

  my_EN_CV <- cva.glmnet(x = C , y = y , intercept = F , alpha = alpha_list ,
                         nlambda = n_lambdas, nfolds = 5)

  my_min_index <- which.min(sapply(my_EN_CV$modlist , function(g){return(min(g$cvm))}))

  my_alpha <- alpha_list[my_min_index]

  if(criterion == 'min'){
    my_lambda <- my_EN_CV$modlist[[my_min_index]]$lambda.min
  }else if(criterion == '1se'){
    my_lambda <- my_EN_CV$modlist[[my_min_index]]$lambda.1se
  }else{
    stop("invalid name of criterion")
  }
  my_EN <- glmnet(x = C , y = y,
                  family = 'gaussian', alpha = my_alpha,
                  lambda = my_lambda,
                  intercept = F)

  my_coefs_indic <- coef(my_EN)[-1] != 0

  my_ret$selected_vars <- my_coefs_indic

  my_HD_OLS<- HD_OLS(x = x, y = y,
                     conf = C[ , my_coefs_indic],
                     min_DF = min_DF)

  my_ret$sign_beta_x <- my_HD_OLS$sign_beta_x

  my_ret$pval <- my_HD_OLS$pval

  return(my_ret)

}



LADS <- function(x,y,conf,criterion = 'min',n_lambdas = 200, min_DF = 20){
  #' Function performing high-dimensional estimation via LASSO Double Selection (Belloni, Cherozhukov and Hansen (2014))
  #' @param x is the exposure vector
  #' @param y is the response vector
  #' @param conf is the confound matrix
  #' @param n_lambdas is the number of lambdas tested for selection using 5-fold cross-validation
  #' @param criterion refers to the criterion for choosing lambda, "min" (default) or "1se"
  #' @param min_DF is the minimal number of degrees of freedom in the final OLS model
  #' @return the p-value and the sign of the coefficient for x
  #' @details If the full model has degrees of freedom above min_DF, the function uses HD_OLS at the last step
  my_ret <- list()


  C <- scale(as.matrix(unname(conf)))

  x <- scale(as.matrix(x))

  y <- scale(as.matrix(y))

  n <- dim(C)[1]

  my_LASSO_CV_y <- cv.glmnet(x = C, y = y, nlambda = n_lambdas, alpha = 1, intercept = F, nfolds = 5)

  if(criterion == 'min'){
    my_lambda_y <- my_LASSO_CV_y$lambda.min
  }else if(criterion == '1se'){
    my_lambda_y <- my_LASSO_CV_y$lambda.1se
  }else{
    stop("invalid name of criterion")
  }

  my_LASSO_y <- glmnet(x = C , y = y,
                       family = 'gaussian', alpha = 1,
                       lambda = my_lambda_y,
                       intercept = F)

  my_coefs_indic_y <- coef(my_LASSO_y)[-1] != 0


  my_LASSO_CV_x <- cv.glmnet(x = C, y = x, nlambda = n_lambdas, alpha = 1, nfolds = 5)

  if(criterion == 'min'){
    my_lambda_x <- my_LASSO_CV_x$lambda.min
  }else if(criterion == '1se'){
    my_lambda_x <- my_LASSO_CV_x$lambda.1se
  }else{
    stop("invalid name of criterion")
  }

  my_LASSO_x <- glmnet(x = C , y = x,
                       family = 'gaussian', alpha = 1,
                       lambda = my_lambda_x,
                       intercept = F)

  my_coefs_indic_x <- coef(my_LASSO_x)[-1] != 0

  my_coefs_indic <- my_coefs_indic_y | my_coefs_indic_x

  my_ret$selected_vars <- my_coefs_indic


  my_HD_OLS <- HD_OLS(x = x, y = y,
                      conf = C[ , my_coefs_indic],
                      min_DF = min_DF)

  my_ret$sign_beta_x <- my_HD_OLS$sign_beta_x

  my_ret$pval <- my_HD_OLS$pval
  return(my_ret)
}




ENDS <- function(x,y,conf,criterion = 'min',n_lambdas = 200, n_alphas = 10, min_DF = 20){
  #' Function performing high-dimensional estimation via LASSO Double Selection (Belloni, Cherozhukov and Hansen (2014))
  #' @param x is the exposure vector
  #' @param y is the response vector
  #' @param conf is the confound matrix
  #' @param n_lambdas is the number of lambdas tested for selection using 5-fold cross-validation
  #' @param n_alphas is the number of alphas tested for selection using 5-fold cross-validation
  #' @param criterion refers to the criterion for choosing lambda, "min" (default) or "1se"
  #' @param min_DF is the minimal number of degrees of freedom in the final OLS model
  #' @return the p-value and the sign of the coefficient for x
  #' @details If the full model has degrees of freedom above min_DF, the function uses HD_OLS at the last step

  my_ret <- list()


  C <- scale(as.matrix(unname(conf)))

  x <- scale(as.matrix(x))

  y <- scale(as.matrix(y))

  n <- dim(C)[1]

  alpha_list <- (1:n_alphas)/(1+n_alphas)

  my_EN_CV_y <- cva.glmnet(x = C , y = y , intercept = F , alpha = alpha_list ,
                           nlambda = n_lambdas, nfolds = 5)

  my_min_index_y <- which.min(sapply(my_EN_CV_y$modlist , function(g){return(min(g$cvm))}))

  my_alpha_y <- alpha_list[my_min_index_y]

  if(criterion == 'min'){
    my_lambda_y <- my_EN_CV_y$modlist[[my_min_index_y]]$lambda.min
  }else if(criterion == '1se'){
    my_lambda_y <- my_EN_CV_y$modlist[[my_min_index_y]]$lambda.1se
  }else{
    stop("invalid name of criterion")
  }

  my_EN_y <- glmnet(x = C , y = y,
                    family = 'gaussian', alpha = my_alpha_y,
                    lambda = my_lambda_y,
                    intercept = F)

  my_coefs_indic_y <- coef(my_EN_y)[-1] != 0


  my_EN_CV_x <- cva.glmnet(x = C , y = x , intercept = F , alpha = alpha_list ,
                           nlambda = n_lambdas, nfolds = 5)

  my_min_index_x <- which.min(sapply(my_EN_CV_x$modlist , function(g){return(min(g$cvm))}))

  my_alpha_x <- alpha_list[my_min_index_x]

  if(criterion == 'min'){
    my_lambda_x <- my_EN_CV_x$modlist[[my_min_index_x]]$lambda.min
  }else if(criterion == '1se'){
    my_lambda_x <- my_EN_CV_x$modlist[[my_min_index_x]]$lambda.1se
  }else{
    stop("invalid name of criterion")
  }

  my_EN_x <- glmnet(x = C , y = x,
                    family = 'gaussian', alpha = my_alpha_x,
                    lambda = my_lambda_x,
                    intercept = F)

  my_coefs_indic_x <- coef(my_EN_x)[-1] != 0

  my_coefs_indic <- my_coefs_indic_y | my_coefs_indic_x

  my_ret$selected_vars <- my_coefs_indic


  my_HD_OLS <- HD_OLS(x = x, y = y,
                      conf = C[ , my_coefs_indic],
                      min_DF = min_DF)

  my_ret$sign_beta_x <- my_HD_OLS$sign_beta_x

  my_ret$pval <- my_HD_OLS$pval
  return(my_ret)
}



DESLA <- function(x,y,conf,inv_method = 'JM', n_lambdas = 200){
  #' Function performing de-sparsified LASSO estimation
  #' @param x is the exposure vector
  #' @param y is the response vector
  #' @param conf is the confound matrix
  #' @param inv_method is the method for estimating the precision matrix, either "JM" (default) (see Javanmard and Montanari (2014)) or "MP" (faster) (see Boot and Nibbering (2017))
  #' @param n_lambdas is the number of lambdas tested for selection using 5-fold cross-validation
  #' @return the p-value and the sign of the coefficient for x

  my_ret <- list()

  my_df <- as.data.frame(unname(conf))

  my_df$x <- x

  my_df$y <- y

  if(inv_method == 'JM'){

    X <- scale(as.matrix(base::subset(my_df , select = -c(y))))

    y <- scale(my_df$y)

    my_SSLasso <- SSLasso(X = X , y = y,verbose = F)

    my_ret$sign_beta_x <- sign(my_SSLasso$unb.coef[['x']])

    my_ret$pval <- my_SSLasso$pvals[['x']]

  }else if(inv_method == 'MP'){
    X <- as.matrix(scale(subset(my_df , select = -c(y))))

    y <- as.matrix(scale(my_df[ , 'y']))

    n <- dim(X)[1]

    p <- dim(X)[2]


    my_LASSO_CV <- cv.glmnet(x = X, y = y, nlambda = n_lambdas, intercept = F, alpha = 1, nfolds = 5)

    my_lambda <- my_LASSO_CV$lambda.min

    my_LASSO <- glmnet(x = X , y = y,
                       family = 'gaussian', alpha = 1,
                       lambda = my_lambda,
                       intercept = F)

    beta_LASSO <- coef(my_LASSO)

    n_vars <- sum(abs(beta_LASSO)!= 0)

    sigma_sq_hat <- sum((y-predict(my_LASSO , X))^2)/(n-n_vars)

    M <- pinv_SVD(t(X)%*%X/(n-1))

    beta_DSL <- as.matrix(beta_LASSO[-1 , 1])+M%*%t(X)%*%(y-predict(my_LASSO , X))/n

    my_ret$sign_beta_x <- sign(beta_DSL['x',1])

    var_DSL <- sigma_sq_hat*M%*%t(X)%*%X%*%t(M)/n^2

    z_score_x <- beta_DSL[p,1]/sqrt(var_DSL[p,p])

    my_ret$pval <- 2*min(pnorm(q = z_score_x) , 1-pnorm(q = z_score_x))

  }else{
    stop("invalid name of inv_method")
  }
  return(my_ret)
}
