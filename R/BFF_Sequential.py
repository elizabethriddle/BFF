#' BFF Change Point Detector Sequential Initial
#'
#' Detects change points within data stream sequentially using a Bayesian adaptive estimation
#' procedure and posterior p-values.
#'
#' @param data data stream to perform change point detection on
#' @param threshold_val vector of thresholds to use to assess whether new data point is a change point
#' @param burnin the intial period to use to build the model on which no change points are detected
#' @param grace_period period after a change in which changes are not detected
#' @param sw sliding window size for p-value calibration
#' @param alpha alpha parameter for lambda beta prior
#' @param beta beta parameter for lambda beta prior
#' @param post_p to calculate posterior p-values
#' @param post_p to calculate predictive posterior p-values
#' @return return the change points detected by the algorithm for each of the dection procedures, for the different thresholds.
#' @export
BFF_gaussian_seq_initial <- function(data,threshold_val = c(0.05),burnin=200,grace_period=20,sw=2000,alpha=39,beta=1.8,post_p = TRUE,pred_p = TRUE){
  N_values_old <- 0
  D_values_old <- 0
  M_values_old <- 0
  N_values_old_old <- 0
  D_values_old_old <- 0
  M_values_old_old <- 0
  mu_values <- c()
  current_sigma2 <- 0
  p_values_change_mu<-c()
  p_values_change_sigma2<-c()
  p_values_change_lambda<-c()
  p_values_predictive<-c()
  mu_0_prior_value <- 5
  sigma2_0_prior_value <- 1
  sigma2_alpha_0_prior_value <- 3.36
  sigma2_beta_0_prior_value <- 8.72
  mu_0_prior_value_old <- 5
  sigma2_0_prior_value_old <- 1
  sigma2_alpha_0_prior_value_old <- 3.36
  sigma2_beta_0_prior_value_old <- 8.72


  anomalous_mu_threshold_sw <- list()
  anomalous_sigma2_threshold_sw <- list()
  anomalous_lambda_threshold_sw <- list()
  anomalous_mu_threshold_sw_uncal <- list()
  anomalous_sigma2_threshold_sw_uncal <- list()
  anomalous_lambda_threshold_sw_uncal <- list()


  anomalous_predictive_threshold <- list()
  anomalous_predictive_threshold_cal <- list()

  for(j in c(1:length(threshold_val))){
    if(post_p == TRUE){
      anomalous_mu_threshold_sw[[j]] <- c(0)
      anomalous_sigma2_threshold_sw[[j]] <- c(0)
      anomalous_lambda_threshold_sw[[j]] <- c(0)

      anomalous_mu_threshold_sw_uncal[[j]] <- c(0)
      anomalous_sigma2_threshold_sw_uncal[[j]] <- c(0)
      anomalous_lambda_threshold_sw_uncal[[j]] <- c(0)
    }
    if(pred_p == TRUE){
      anomalous_predictive_threshold[[j]] <- c(0)
      anomalous_predictive_threshold_cal[[j]] <- c(0)
    }
  }


  mu_MAE <- 0

  p_values_calibrated_mu_sw <- c()
  p_values_calibrated_sigma2_sw <- c()
  p_values_calibrated_lambda_sw<- c()
  p_values_calibrated_predictive <- c()

  for(i in c(1:length(data))){
    model_results <- gaussian_unknown_mean_var_BFF_sequential(data[i],N_values_old,D_values_old,M_values_old,mu_0_prior_value,sigma2_0_prior_value,sigma2_alpha_0_prior_value,sigma2_beta_0_prior_value,alpha,beta)

    if(i>1){
      if(post_p == TRUE){
        if(model_results$mu_estimate<current_mu){
          p_values_change_mu<-c(p_values_change_mu,2*pnorm(model_results$mu_estimate,current_mu,sqrt(current_sigma2/(D_values_old+(1/sigma2_0_prior_value)))))
        }else{
          p_values_change_mu<-c(p_values_change_mu,2*pnorm(model_results$mu_estimate,current_mu,sqrt(current_sigma2/(D_values_old+(1/sigma2_0_prior_value))),lower.tail = FALSE))
        }
        # find sigma posterior parameters:
        sigma_2_alpha_0 <- (0.5*D_values_old)+sigma2_alpha_0_prior_value
        sigma_2_beta_0 <- sigma2_beta_0_prior_value + 0.5*((mu_0_prior_value^2/sigma2_0_prior_value)+M_values_old-((N_values_old+(mu_0_prior_value/sigma2_0_prior_value))^2/(D_values_old+(1/sigma2_0_prior_value))))
        p_values_change_lambda<-c(p_values_change_lambda,p_value_calculator_lambda(lambda_posterior_unk_mu_sigma2_sequential_2_beta_prior_unnorm,model_results$lambda_estimate,data_new=data[i-1],N_prev=N_values_old_old,D_prev=D_values_old_old,M_prev=M_values_old_old,mu_0=mu_0_prior_value_old,sigma2_0=sigma2_0_prior_value_old,alpha_0=sigma2_alpha_0_prior_value_old,beta_0=sigma2_beta_0_prior_value_old,alpha=alpha,beta=beta))


      }
      if(pred_p == TRUE){
        p_values_predictive <- c(p_values_predictive,BFF_predictive_p_value_calculator_approx(data[i],current_mu,N_values_old,D_values_old,M_values_old,mu_0_prior_value_old,sigma2_0_prior_value_old,sigma2_alpha_0_prior_value_old,sigma2_beta_0_prior_value_old))
        # Detect Anom
        for(j in c(1:length(threshold_val))){
          anomalous_predictive_threshold[[j]] <- threshold_p_value_detector_seq(i,tail(p_values_predictive,1),anomalous_predictive_threshold[[j]],threshold_val[j],grace_period)
        }
      }
    }
    if(i>=(sw+2)){
      if(post_p == TRUE){
        p_values_ecdf_mu <- ecdf(p_values_change_mu[c(1:sw)+(length(p_values_change_mu)-(sw+1))])
        for(j in c(1:length(threshold_val))){
          anomalous_mu_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_mu(p_values_change_mu[length(p_values_change_mu)]),anomalous_mu_threshold_sw[[j]],threshold_val[j],grace_period)
        }


        p_values_ecdf_lambda <- ecdf(p_values_change_lambda[c(1:sw)+(length(p_values_change_lambda)-(sw+1))])
        for(j in c(1:length(threshold_val))){
          anomalous_lambda_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_lambda(p_values_change_lambda[length(p_values_change_lambda)]),anomalous_lambda_threshold_sw[[j]],threshold_val[j],grace_period)
        }

        # # Anom detection uncalibrated p-values
        for(j in c(1:length(threshold_val))){
          anomalous_mu_threshold_sw_uncal[[j]] <- threshold_p_value_detector_seq(i,tail(p_values_change_mu,1),anomalous_mu_threshold_sw_uncal[[j]],threshold_val,grace_period)
        }
        for(j in c(1:length(threshold_val))){
          anomalous_lambda_threshold_sw_uncal[[j]] <- threshold_p_value_detector_seq(i,tail(p_values_change_lambda,1),anomalous_lambda_threshold_sw_uncal[[j]],threshold_val,grace_period)
        }


        # remove a value from p_values
        p_values_change_mu <- p_values_change_mu[-1]
        p_values_change_lambda <- p_values_change_lambda[-1]
      }
      if(pred_p == TRUE){
        p_values_ecdf <- ecdf(p_values_predictive[c(1:sw)+(length(p_values_predictive)-(sw+1))])
        # Detect Anom
        for(j in c(1:length(threshold_val))){
          anomalous_predictive_threshold_cal[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf(p_values_predictive[length(p_values_predictive)]),anomalous_predictive_threshold_cal[[j]],threshold_val[j],grace_period)
        }

        p_values_predictive <- p_values_predictive[-1]
      }
    }
    current_mu <- model_results$mu_estimate
    current_sigma2 <- model_results$sigma2_estimate
    N_values_old_old <- N_values_old
    D_values_old_old <- D_values_old
    M_values_old_old <- M_values_old
    N_values_old <- model_results$N_new
    D_values_old <- model_results$D_new
    M_values_old <- model_results$M_new
    mu_0_prior_value_old <- mu_0_prior_value
    sigma2_0_prior_value_old <- sigma2_0_prior_value
    sigma2_alpha_0_prior_value_old <- sigma2_alpha_0_prior_value
    sigma2_beta_0_prior_value_old <- sigma2_beta_0_prior_value
    mu_0_prior_value <- model_results$mu_estimate
    sigma2_alpha_0_prior_value <- 1/2
    sigma2_beta_0_prior_value <- (2/3)*model_results$sigma2_estimate

  }

  if((post_p == TRUE) & (pred_p == FALSE)){
    return(list(
      anomalous_mu_threshold_sw = anomalous_mu_threshold_sw,
      anomalous_lambda_threshold_sw = anomalous_lambda_threshold_sw,
      p_values_change_mu=p_values_change_mu,
      p_values_change_lambda=p_values_change_lambda,
      current_mu = current_mu,
      current_sigma2 = current_sigma2,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      M_values_old_old = M_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      M_values_old = M_values_old,
      mu_0_prior_value_old = mu_0_prior_value_old,
      sigma2_0_prior_value_old = sigma2_0_prior_value_old,
      sigma2_alpha_0_prior_value_old = sigma2_alpha_0_prior_value_old,
      sigma2_beta_0_prior_value_old = sigma2_beta_0_prior_value_old,
      mu_0_prior_value = mu_0_prior_value,
      sigma2_alpha_0_prior_value = sigma2_alpha_0_prior_value,
      sigma2_beta_0_prior_value = sigma2_beta_0_prior_value
    ))
  }
  else if((post_p == FALSE) & (pred_p == TRUE)){
    return(list(
      anomalous_predictive_threshold_cal = anomalous_predictive_threshold_cal,
      p_values_predictive = p_values_predictive,
      current_mu = current_mu,
      current_sigma2 = current_sigma2,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      M_values_old_old = M_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      M_values_old = M_values_old,
      mu_0_prior_value_old = mu_0_prior_value_old,
      sigma2_0_prior_value_old = sigma2_0_prior_value_old,
      sigma2_alpha_0_prior_value_old = sigma2_alpha_0_prior_value_old,
      sigma2_beta_0_prior_value_old = sigma2_beta_0_prior_value_old,
      mu_0_prior_value = mu_0_prior_value,
      sigma2_alpha_0_prior_value = sigma2_alpha_0_prior_value,
      sigma2_beta_0_prior_value = sigma2_beta_0_prior_value
    ))
  }
  else if((post_p == TRUE) & (pred_p == TRUE)){
    return(list(
      anomalous_mu_threshold_sw = anomalous_mu_threshold_sw,
      anomalous_lambda_threshold_sw = anomalous_lambda_threshold_sw,
      p_values_change_mu=p_values_change_mu,
      p_values_change_lambda=p_values_change_lambda,
      anomalous_predictive_threshold_cal = anomalous_predictive_threshold_cal,
      p_values_predictive = p_values_predictive,
      current_mu = current_mu,
      current_sigma2 = current_sigma2,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      M_values_old_old = M_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      M_values_old = M_values_old,
      mu_0_prior_value_old = mu_0_prior_value_old,
      sigma2_0_prior_value_old = sigma2_0_prior_value_old,
      sigma2_alpha_0_prior_value_old = sigma2_alpha_0_prior_value_old,
      sigma2_beta_0_prior_value_old = sigma2_beta_0_prior_value_old,
      mu_0_prior_value = mu_0_prior_value,
      sigma2_alpha_0_prior_value = sigma2_alpha_0_prior_value,
      sigma2_beta_0_prior_value = sigma2_beta_0_prior_value
    ))

  }




}




#' BFF Change Point Detector Sequential Update
#'
#' Detects change points within data stream sequentially using a Bayesian adaptive estimation
#' procedure and posterior p-values.
#'
#' @param data new data point
#' @param i data point number
#' @param threshold_val vector of thresholds to use to assess whether new data point is a change point
#' @param grace_period period after a change in which changes are not detected
#' @param sw sliding window size for p-value calibration
#' @param alpha alpha parameter for lambda beta prior
#' @param beta beta parameter for lambda beta prior
#' @param post_p to calculate posterior p-values
#' @param post_p to calculate predictive posterior p-values
#' @param current_mu estimated mu from previous time
#' @param current_sigma2 estimated sigma2 from previous time
#' @param N_values_old_old old old N value
#' @param D_values_old_old old old D value
#' @param M_values_old_old old old M value
#' @param N_values_old old N value
#' @param D_values_old old D value
#' @param M_values_old old M value
#' @param mu_0_prior_value_old old prior mean for mu
#' @param sigma2_0_prior_value_old old prior var for mu
#' @param sigma2_alpha_0_prior_value_old old prior alpha for sigma2
#' @param sigma2_beta_0_prior_value_old old prior beta for sigma2
#' @param mu_0_prior_value prior mean for mu
#' @param sigma2_alpha_0_prior_value prior alpha for sigma2
#' @param sigma2_beta_0_prior_value prior beta for sigma2
#' @param anomalous_mu_threshold_sw current anomalous mu posterior p-values list
#' @param anomalous_lambda_threshold_sw current anomalous lambda posterior p-values list
#' @param p_values_change_mu uncalibrated mu posterior p-values
#' @param p_values_change_lambda uncalibrated lambda posterior p-values
#' @param anomalous_predictive_threshold_cal current anomalous predictive posterior p-values list
#' @param p_values_predictive uncalibrated predictive posterior p-values
#' @return return the change points detected by the algorithm for each of the dection procedures, for the different thresholds.
#' @export
BFF_gaussian_seq_update <- function(data,i,threshold_val = c(0.05),grace_period=20,sw=2000,alpha=39,beta=1.8,post_p = TRUE,pred_p = TRUE,
                                    current_mu,current_sigma2,N_values_old_old,D_values_old_old,M_values_old_old,N_values_old,D_values_old,M_values_old,
                                    mu_0_prior_value_old,sigma2_0_prior_value_old,sigma2_alpha_0_prior_value_old,sigma2_beta_0_prior_value_old,
                                    mu_0_prior_value,sigma2_alpha_0_prior_value,sigma2_beta_0_prior_value,
                                    anomalous_mu_threshold_sw,anomalous_lambda_threshold_sw,p_values_change_mu,
                                    p_values_change_lambda,anomalous_predictive_threshold_cal,p_values_predictive){

  sigma2_0_prior_value <- 1

  model_results <- gaussian_unknown_mean_var_BFF_sequential(data[i],N_values_old,D_values_old,M_values_old,mu_0_prior_value,sigma2_0_prior_value,sigma2_alpha_0_prior_value,sigma2_beta_0_prior_value,alpha,beta)


  sigma2est <- model_results$sigma2_estimate/(model_results$D_new+(1/sigma2_0_prior_value))

  if(post_p == TRUE){
    p_values_ecdf_mu <- ecdf(p_values_change_mu[c(1:sw)+(length(p_values_change_mu)-(sw+1))])
    for(j in c(1:length(threshold_val))){
      anomalous_mu_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_mu(p_values_change_mu[length(p_values_change_mu)]),anomalous_mu_threshold_sw[[j]],threshold_val[j],grace_period)
    }


    p_values_ecdf_lambda <- ecdf(p_values_change_lambda[c(1:sw)+(length(p_values_change_lambda)-(sw+1))])
    for(j in c(1:length(threshold_val))){
      anomalous_lambda_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_lambda(p_values_change_lambda[length(p_values_change_lambda)]),anomalous_lambda_threshold_sw[[j]],threshold_val[j],grace_period)
    }


    # remove a value from p_values
    p_values_change_mu <- p_values_change_mu[-1]
    p_values_change_lambda <- p_values_change_lambda[-1]
  }
  if(pred_p == TRUE){
    p_values_ecdf <- ecdf(p_values_predictive[c(1:sw)+(length(p_values_predictive)-(sw+1))])
    # Detect Anom
    for(j in c(1:length(threshold_val))){
      anomalous_predictive_threshold_cal[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf(p_values_predictive[length(p_values_predictive)]),anomalous_predictive_threshold_cal[[j]],threshold_val[j],grace_period)
    }

    p_values_predictive <- p_values_predictive[-1]
  }

  current_mu <- model_results$mu_estimate
  current_sigma2 <- model_results$sigma2_estimate
  N_values_old_old <- N_values_old
  D_values_old_old <- D_values_old
  M_values_old_old <- M_values_old
  N_values_old <- model_results$N_new
  D_values_old <- model_results$D_new
  M_values_old <- model_results$M_new
  mu_0_prior_value_old <- mu_0_prior_value
  sigma2_0_prior_value_old <- sigma2_0_prior_value
  sigma2_alpha_0_prior_value_old <- sigma2_alpha_0_prior_value
  sigma2_beta_0_prior_value_old <- sigma2_beta_0_prior_value
  mu_0_prior_value <- model_results$mu_estimate
  sigma2_alpha_0_prior_value <- 1/2
  sigma2_beta_0_prior_value <- (2/3)*model_results$sigma2_estimate



  if((post_p == TRUE) & (pred_p == FALSE)){
    return(list(
      anomalous_mu_threshold_sw = anomalous_mu_threshold_sw,
      anomalous_lambda_threshold_sw = anomalous_lambda_threshold_sw,
      p_values_change_mu=p_values_change_mu,
      p_values_change_lambda=p_values_change_lambda,
      current_mu = current_mu,
      current_sigma2 = current_sigma2,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      M_values_old_old = M_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      M_values_old = M_values_old,
      mu_0_prior_value_old = mu_0_prior_value_old,
      sigma2_0_prior_value_old = sigma2_0_prior_value_old,
      sigma2_alpha_0_prior_value_old = sigma2_alpha_0_prior_value_old,
      sigma2_beta_0_prior_value_old = sigma2_beta_0_prior_value_old,
      mu_0_prior_value = mu_0_prior_value,
      sigma2_alpha_0_prior_value = sigma2_alpha_0_prior_value,
      sigma2_beta_0_prior_value = sigma2_beta_0_prior_value
    ))
  }
  else if((post_p == FALSE) & (pred_p == TRUE)){
    return(list(
      anomalous_predictive_threshold_cal = anomalous_predictive_threshold_cal,
      p_values_predictive = p_values_predictive,
      current_mu = current_mu,
      current_sigma2 = current_sigma2,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      M_values_old_old = M_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      M_values_old = M_values_old,
      mu_0_prior_value_old = mu_0_prior_value_old,
      sigma2_0_prior_value_old = sigma2_0_prior_value_old,
      sigma2_alpha_0_prior_value_old = sigma2_alpha_0_prior_value_old,
      sigma2_beta_0_prior_value_old = sigma2_beta_0_prior_value_old,
      mu_0_prior_value = mu_0_prior_value,
      sigma2_alpha_0_prior_value = sigma2_alpha_0_prior_value,
      sigma2_beta_0_prior_value = sigma2_beta_0_prior_value
    ))
  }
  else if((post_p == TRUE) & (pred_p == TRUE)){
    return(list(
      anomalous_mu_threshold_sw = anomalous_mu_threshold_sw,
      anomalous_lambda_threshold_sw = anomalous_lambda_threshold_sw,
      p_values_change_mu=p_values_change_mu,
      p_values_change_lambda=p_values_change_lambda,
      anomalous_predictive_threshold_cal = anomalous_predictive_threshold_cal,
      p_values_predictive = p_values_predictive,
      current_mu = current_mu,
      current_sigma2 = current_sigma2,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      M_values_old_old = M_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      M_values_old = M_values_old,
      mu_0_prior_value_old = mu_0_prior_value_old,
      sigma2_0_prior_value_old = sigma2_0_prior_value_old,
      sigma2_alpha_0_prior_value_old = sigma2_alpha_0_prior_value_old,
      sigma2_beta_0_prior_value_old = sigma2_beta_0_prior_value_old,
      mu_0_prior_value = mu_0_prior_value,
      sigma2_alpha_0_prior_value = sigma2_alpha_0_prior_value,
      sigma2_beta_0_prior_value = sigma2_beta_0_prior_value
    ))

  }




}








#' BFF Change Point Detector Sequential Poisson Initial
#'
#' Detects change points within data stream sequentially using a Bayesian adaptive estimation
#' procedure and predictive posterior p-values. Additionally estimates parameters
#'
#' @param data data stream to perform change point detection on
#' @param threshold_val vector of thresholds to use to assess whether new data point is a change point
#' @param burnin the intial period to use to build the model on which no change points are detected
#' @param grace_period period after a change in which changes are not detected
#' @param alpha alpha parameter for lambda beta prior
#' @param beta beta parameter for lambda beta prior
#' @param sw sliding window size for p-value calibration
#' @param param_est whether to estimate paramters
#' @return return the change points detected by the algorithm for each of the dection procedures, for the different thresholds. Collects p-values and model parameter estimates.
#' @export
BFF_poisson_seq_initial <- function(data,threshold_val = c(0.05),burnin=200,grace_period=20,sw=2000,alpha=39,beta=1.8,post_p = TRUE,pred_p = TRUE){
  N_values_old <- 0
  D_values_old <- 0
  F_values_old <- 0
  N_values_old_old <- 0
  D_values_old_old <- 0
  F_values_old_old <- 0

  current_gamma <- 0
  p_values_change_gamma<-c()
  p_values_change_lambda<-c()
  p_values_change_pred <- c()
  alpha_0_prior_value <- 2
  beta_0_prior_value <- 1
  alpha_0_prior_value_old <- 2
  beta_0_prior_value_old <- 1
  alpha_prior_value <- alpha
  beta_prior_value <- beta


  anomalous_gamma_threshold_sw <- list()
  anomalous_lambda_threshold_sw <- list()
  anomalous_pred_threshold_sw <- list()


  for(j in c(1:length(threshold_val))){
    anomalous_gamma_threshold_sw[[j]] <- c(0)
    anomalous_lambda_threshold_sw[[j]] <- c(0)
    anomalous_pred_threshold_sw[[j]] <- c(0)
  }

  anomalous_gamma_threshold_sw_uncal <- list()
  anomalous_lambda_threshold_sw_uncal <- list()
  anomalous_pred_threshold_sw_uncal <- list()

  for(j in c(1:length(threshold_val))){
    anomalous_gamma_threshold_sw_uncal[[j]] <- c(0)
    anomalous_lambda_threshold_sw_uncal[[j]] <- c(0)
    anomalous_pred_threshold_sw_uncal[[j]] <- c(0)
  }

  gamma_values <- c()
  lambda_values <- c()


  for(i in c(1:length(data))){
    model_results <- poisson_unknown_rate_BFF_sequential(data[i],N_values_old,D_values_old,F_values_old,alpha_0_prior_value,beta_0_prior_value,alpha_prior_value,beta_prior_value)
    #print(c(i,tail(N_values,1),tail(D_values,1),tail(M_values,1),tail(lambda_values,1)))
    gamma_values <- c(gamma_values,model_results$gamma_estimate)
    lambda_values <- c(lambda_values,model_results$lambda_estimate)

    if(i>1){
      gamma_alpha_0 <- (N_values_old)+alpha_0_prior_value
      gamma_beta_0 <- beta_0_prior_value + D_values_old
      if(post_p == TRUE){
        if(model_results$gamma_estimate>((gamma_alpha_0-1)/(gamma_beta_0))){
          p_values_change_gamma<-c(p_values_change_gamma,pgamma(model_results$gamma_estimate,gamma_alpha_0,gamma_beta_0,lower.tail = F))
        }else{
          p_values_change_gamma<-c(p_values_change_gamma,pgamma(model_results$gamma_estimate,gamma_alpha_0,gamma_beta_0,lower.tail = T))
        }
        p_values_change_lambda<-c(p_values_change_lambda,p_value_calculator_lambda(lambda_posterior_poisson_unk_rate_sequential_beta_prior_unnorm,model_results$lambda_estimate,data_new=data[i-1],N_prev=N_values_old_old,D_prev=D_values_old_old,F_prev=F_values_old_old,alpha_0=alpha_0_prior_value_old,beta_0=beta_0_prior_value_old,alpha=39,beta=1.8))
      }
      if(pred_p == TRUE){
        p_values_change_pred <- c(p_values_change_pred,BFF_predictive_p_value_calculator_approx_poisson(x_new=data[i],N_value=N_values_old,D_value=D_values_old,alpha_0_prior_value=alpha_0_prior_value_old,beta_0_prior_value=beta_0_prior_value_old))
      }
    }

    if(i>=(sw+2)){
      if(post_p == TRUE){
        p_values_ecdf_gamma <- ecdf(p_values_change_gamma[c(1:sw)+(length(p_values_change_gamma)-(sw+1))])
        for(j in c(1:length(threshold_val))){
          anomalous_gamma_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_gamma(p_values_change_gamma[length(p_values_change_gamma)]),anomalous_gamma_threshold_sw[[j]],threshold_val[j],grace_period)
        }

        p_values_ecdf_lambda <- ecdf(p_values_change_lambda[c(1:sw)+(length(p_values_change_lambda)-(sw+1))])
        for(j in c(1:length(threshold_val))){
          anomalous_lambda_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_lambda(p_values_change_lambda[length(p_values_change_lambda)]),anomalous_lambda_threshold_sw[[j]],threshold_val[j],grace_period)
        }
      }
      if(pred_p == TRUE){
        p_values_ecdf_pred <- ecdf(p_values_change_pred[c(1:sw)+(length(p_values_change_pred)-(sw+1))])
        for(j in c(1:length(threshold_val))){
          anomalous_pred_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_pred(p_values_change_pred[length(p_values_change_pred)]),anomalous_pred_threshold_sw[[j]],threshold_val[j],grace_period)
        }
      }

      # remove a value from p_values
      if(post_p == TRUE){
        p_values_change_gamma <- p_values_change_gamma[-1]
        p_values_change_lambda <- p_values_change_lambda[-1]
      }
      if(pred_p == TRUE){
        p_values_change_pred <- p_values_change_pred[-1]
      }
    }

    N_values_old_old <- N_values_old
    D_values_old_old <- D_values_old
    F_values_old_old <- F_values_old
    N_values_old <- model_results$N_new
    D_values_old <- model_results$D_new
    F_values_old <- model_results$F_new
    alpha_0_prior_value_old <- alpha_0_prior_value
    beta_0_prior_value_old <- beta_0_prior_value
    alpha_0_prior_value <- model_results$gamma_estimate+1
    beta_0_prior_value <- 1

  }
  if((post_p == TRUE) & (pred_p == TRUE)){
    return(list(
      anomalous_gamma_threshold_sw = anomalous_gamma_threshold_sw,
      anomalous_lambda_threshold_sw = anomalous_lambda_threshold_sw,
      p_values_change_gamma = p_values_change_gamma,
      p_values_change_lambda =  p_values_change_lambda,
      anomalous_pred_threshold_sw = anomalous_pred_threshold_sw,
      p_values_change_pred = p_values_change_pred,
      gamma_values = gamma_values,
      lambda_values = lambda_values,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      F_values_old_old = F_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      F_values_old = F_values_old,
      alpha_0_prior_value_old = alpha_0_prior_value_old,
      beta_0_prior_value_old = beta_0_prior_value_old,
      alpha_0_prior_value = alpha_0_prior_value,
      beta_0_prior_value = beta_0_prior_value
    ))
  }
  if((post_p == TRUE) & (pred_p == FALSE)){
    return(list(
      anomalous_gamma_threshold_sw = anomalous_gamma_threshold_sw,
      anomalous_lambda_threshold_sw = anomalous_lambda_threshold_sw,
      p_values_change_gamma = p_values_change_gamma,
      p_values_change_lambda =  p_values_change_lambda,
      gamma_values = gamma_values,
      lambda_values = lambda_values,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      F_values_old_old = F_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      F_values_old = F_values_old,
      alpha_0_prior_value_old = alpha_0_prior_value_old,
      beta_0_prior_value_old = beta_0_prior_value_old,
      alpha_0_prior_value = alpha_0_prior_value,
      beta_0_prior_value = beta_0_prior_value
    ))
  }
  if((post_p == FALSE) & (pred_p == TRUE)){
    return(list(
      anomalous_pred_threshold_sw = anomalous_pred_threshold_sw,
      p_values_change_pred = p_values_change_pred,
      gamma_values = gamma_values,
      lambda_values = lambda_values,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      F_values_old_old = F_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      F_values_old = F_values_old,
      alpha_0_prior_value_old = alpha_0_prior_value_old,
      beta_0_prior_value_old = beta_0_prior_value_old,
      alpha_0_prior_value = alpha_0_prior_value,
      beta_0_prior_value = beta_0_prior_value
    ))
  }


}



#' BFF Change Point Detector Sequential Poisson Update
#'
#' Detects change points within data stream sequentially using a Bayesian adaptive estimation
#' procedure and predictive posterior p-values. Additionally estimates parameters
#'
#' @param data new data point
#' @param i data point number
#' @param threshold_val vector of thresholds to use to assess whether new data point is a change point
#' @param grace_period period after a change in which changes are not detected
#' @param sw sliding window size for p-value calibration
#' @param alpha alpha parameter for lambda beta prior
#' @param beta beta parameter for lambda beta prior
#' @param post_p to calculate posterior p-values
#' @param post_p to calculate predictive posterior p-values
#' @param current_gamma estimated gamma from previous time
#' @param current_lambda estimated lambda from previous time
#' @param N_values_old_old old old N value
#' @param D_values_old_old old old D value
#' @param F_values_old_old old old F value
#' @param N_values_old old N value
#' @param D_values_old old D value
#' @param F_values_old old F value
#' @param alpha_0_prior_value_old old prior alpha for gamma
#' @param beta_0_prior_value_old old prior beta for gamma
#' @param alpha_0_prior_value prior alpha for gamma
#' @param beta_0_prior_value prior beta for gamma
#' @param anomalous_gamma_threshold_sw current anomalous gamma posterior p-values list
#' @param anomalous_lambda_threshold_sw current anomalous lambda posterior p-values list
#' @param p_values_change_gamma uncalibrated gamma posterior p-values
#' @param p_values_change_lambda uncalibrated lambda posterior p-values
#' @param anomalous_pred_threshold_sw current anomalous predictive posterior p-values list
#' @param p_values_change_pred uncalibrated predictive posterior p-values
#' @return return the change points detected by the algorithm for each of the dection procedures, for the different thresholds.
#' @export


BFF_poisson_seq_update <- function(data,i,threshold_val = c(0.05),grace_period=20,sw=2000,alpha=39,beta=1.8,post_p = TRUE,pred_p = TRUE,
                                    current_gamma,current_lambda,N_values_old_old,D_values_old_old,F_values_old_old,N_values_old,D_values_old,F_values_old,
                                    alpha_0_prior_value_old,beta_0_prior_value_old,
                                    alpha_0_prior_value,beta_0_prior_value,
                                    anomalous_gamma_threshold_sw,anomalous_lambda_threshold_sw,p_values_change_gamma,
                                    p_values_change_lambda,anomalous_pred_threshold_sw,p_values_change_pred){

  model_results <- poisson_unknown_rate_BFF_sequential(data[i],N_values_old,D_values_old,F_values_old,alpha_0_prior_value,beta_0_prior_value,alpha_prior_value,beta_prior_value)
  gamma_values <- c(gamma_values,model_results$gamma_estimate)
  lambda_values <- c(lambda_values,model_results$lambda_estimate)

  gamma_alpha_0 <- (N_values_old)+alpha_0_prior_value
  gamma_beta_0 <- beta_0_prior_value + D_values_old
  if(post_p == TRUE){
    if(model_results$gamma_estimate>((gamma_alpha_0-1)/(gamma_beta_0))){
      p_values_change_gamma<-c(p_values_change_gamma,pgamma(model_results$gamma_estimate,gamma_alpha_0,gamma_beta_0,lower.tail = F))
    }else{
      p_values_change_gamma<-c(p_values_change_gamma,pgamma(model_results$gamma_estimate,gamma_alpha_0,gamma_beta_0,lower.tail = T))
    }
    p_values_change_lambda<-c(p_values_change_lambda,p_value_calculator_lambda(lambda_posterior_poisson_unk_rate_sequential_beta_prior_unnorm,model_results$lambda_estimate,data_new=data[i-1],N_prev=N_values_old_old,D_prev=D_values_old_old,F_prev=F_values_old_old,alpha_0=alpha_0_prior_value_old,beta_0=beta_0_prior_value_old,alpha=39,beta=1.8))
  }
  if(pred_p == TRUE){
    p_values_change_pred <- c(p_values_change_pred,BFF_predictive_p_value_calculator_approx_poisson(x_new=data[i],N_value=N_values_old,D_value=D_values_old,alpha_0_prior_value=alpha_0_prior_value_old,beta_0_prior_value=beta_0_prior_value_old))
  }


  if(post_p == TRUE){
    p_values_ecdf_gamma <- ecdf(p_values_change_gamma[c(1:sw)+(length(p_values_change_gamma)-(sw+1))])
    for(j in c(1:length(threshold_val))){
      anomalous_gamma_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_gamma(p_values_change_gamma[length(p_values_change_gamma)]),anomalous_gamma_threshold_sw[[j]],threshold_val[j],grace_period)
    }

    p_values_ecdf_lambda <- ecdf(p_values_change_lambda[c(1:sw)+(length(p_values_change_lambda)-(sw+1))])
    for(j in c(1:length(threshold_val))){
      anomalous_lambda_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_lambda(p_values_change_lambda[length(p_values_change_lambda)]),anomalous_lambda_threshold_sw[[j]],threshold_val[j],grace_period)
    }
  }
  if(pred_p == TRUE){
    p_values_ecdf_pred <- ecdf(p_values_change_pred[c(1:sw)+(length(p_values_change_pred)-(sw+1))])
    for(j in c(1:length(threshold_val))){
      anomalous_pred_threshold_sw[[j]] <- threshold_p_value_detector_seq(i,p_values_ecdf_pred(p_values_change_pred[length(p_values_change_pred)]),anomalous_pred_threshold_sw[[j]],threshold_val[j],grace_period)
    }
  }

  # remove a value from p_values
  if(post_p == TRUE){
    p_values_change_gamma <- p_values_change_gamma[-1]
    p_values_change_lambda <- p_values_change_lambda[-1]
  }
  if(pred_p == TRUE){
    p_values_change_pred <- p_values_change_pred[-1]
  }


  N_values_old_old <- N_values_old
  D_values_old_old <- D_values_old
  F_values_old_old <- F_values_old
  N_values_old <- model_results$N_new
  D_values_old <- model_results$D_new
  F_values_old <- model_results$F_new
  alpha_0_prior_value_old <- alpha_0_prior_value
  beta_0_prior_value_old <- beta_0_prior_value
  alpha_0_prior_value <- model_results$gamma_estimate+1
  beta_0_prior_value <- 1


  if((post_p == TRUE) & (pred_p == TRUE)){
    return(list(
      anomalous_gamma_threshold_sw = anomalous_gamma_threshold_sw,
      anomalous_lambda_threshold_sw = anomalous_lambda_threshold_sw,
      p_values_change_gamma = p_values_change_gamma,
      p_values_change_lambda =  p_values_change_lambda,
      anomalous_pred_threshold_sw = anomalous_pred_threshold_sw,
      p_values_change_pred = p_values_change_pred,
      gamma_values = gamma_values,
      lambda_values = lambda_values,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      F_values_old_old = F_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      F_values_old = F_values_old,
      mu_0_prior_value_old = mu_0_prior_value_old,
      alpha_0_prior_value_old = alpha_0_prior_value_old,
      beta_0_prior_value_old = beta_0_prior_value_old,
      alpha_0_prior_value = alpha_0_prior_value,
      beta_0_prior_value = beta_0_prior_value
    ))
  }
  if((post_p == TRUE) & (pred_p == FALSE)){
    return(list(
      anomalous_gamma_threshold_sw = anomalous_gamma_threshold_sw,
      anomalous_lambda_threshold_sw = anomalous_lambda_threshold_sw,
      p_values_change_gamma = p_values_change_gamma,
      p_values_change_lambda =  p_values_change_lambda,
      gamma_values = gamma_values,
      lambda_values = lambda_values,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      F_values_old_old = F_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      F_values_old = F_values_old,
      mu_0_prior_value_old = mu_0_prior_value_old,
      alpha_0_prior_value_old = alpha_0_prior_value_old,
      beta_0_prior_value_old = beta_0_prior_value_old,
      alpha_0_prior_value = alpha_0_prior_value,
      beta_0_prior_value = beta_0_prior_value
    ))
  }
  if((post_p == FALSE) & (pred_p == TRUE)){
    return(list(
      anomalous_pred_threshold_sw = anomalous_pred_threshold_sw,
      p_values_change_pred = p_values_change_pred,
      gamma_values = gamma_values,
      lambda_values = lambda_values,
      N_values_old_old = N_values_old_old,
      D_values_old_old = D_values_old_old,
      F_values_old_old = F_values_old_old,
      N_values_old = N_values_old,
      D_values_old = D_values_old,
      F_values_old = F_values_old,
      mu_0_prior_value_old = mu_0_prior_value_old,
      alpha_0_prior_value_old = alpha_0_prior_value_old,
      beta_0_prior_value_old = beta_0_prior_value_old,
      alpha_0_prior_value = alpha_0_prior_value,
      beta_0_prior_value = beta_0_prior_value
    ))
  }


}
