    devtools::install_github("elizabethriddle/BFF")

    ## Error : 'format_warning' is not an exported object from 'namespace:cli'
    ## Error : 'format_warning' is not an exported object from 'namespace:cli'
    ## Error : 'format_warning' is not an exported object from 'namespace:cli'

    ## Downloading GitHub repo elizabethriddle/BFF@master

    ## Error in utils::download.file(url, path, method = method, quiet = quiet,  : 
    ##   cannot open URL 'https://api.github.com/repos/elizabethriddle/BFF/tarball/master'

    #devtools::install("~/Documents/PhD Year 1/RCode/Bayesian Forgetting Factor/BFF")
    library(BFF)

    ## Loading required package: Matrix

    ## Loading required package: ggplot2

    ## Loading required package: grid

    ## Loading required package: gridExtra

    ## Loading required package: matrixStats

    ## Loading required package: plyr

    ## 
    ## Attaching package: 'plyr'

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## Loading required package: invgamma

    ## Loading required package: rootSolve

    ## Loading required package: gmm

    ## Loading required package: sandwich

    ## Loading required package: ks

    ## Loading required package: ModelMetrics

    ## 
    ## Attaching package: 'ModelMetrics'

    ## The following object is masked from 'package:base':
    ## 
    ##     kappa

    ## Loading required package: changepoint

    ## Loading required package: zoo

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

    ## Successfully loaded changepoint package version 2.2.2
    ##  NOTE: Predefined penalty values changed in version 2.2.  Previous penalty values with a postfix 1 i.e. SIC1 are now without i.e. SIC and previous penalties without a postfix i.e. SIC are now with a postfix 0 i.e. SIC0. See NEWS and help files for further details.

    ## Loading required package: wbs

    ## Loading required package: ocp

    ## Loading required package: ffstream

    ## Loading required package: Rcpp

    ## Loading required package: DeCAFS

    ## Loading required package: LaplacesDemon

    ## 
    ## Attaching package: 'LaplacesDemon'

    ## The following objects are masked from 'package:invgamma':
    ## 
    ##     dinvchisq, dinvgamma, rinvchisq, rinvgamma

Illustrative Batch Example
==========================

This simulation is found in the thesis with a single change point and a
trend period.

    set.seed(10)
    noisy_lambda_illustrative_data <- c(rnorm(200,0,1),rnorm(100,5,1),rnorm(50,seq(5,-5,length.out = 50),1),rnorm(150,-5,1))
    results_BFF_post_param_illus <- BFF_gaussian_predictive_paramest(noisy_lambda_illustrative_data,threshold_val = c(0.05),burnin=100,grace_period = 10,alpha=39,beta=1.8,sw=100)


    g1<-ggplot()+ theme_bw()+geom_line(data=data.frame(x=c(1:500),y=noisy_lambda_illustrative_data),aes(x,y))+geom_point(data=data.frame(x=c(200),y=min(noisy_lambda_illustrative_data[101:500])-0.1),shape=4,size=4,aes(x,y))+geom_point(data=data.frame(x=300,y=min(noisy_lambda_illustrative_data[101:500])-0.1),shape=1,size=4,aes(x,y))+xlab("Observation")+ylim(min(noisy_lambda_illustrative_data[101:500])-0.1,max(noisy_lambda_illustrative_data[101:500]))+ylab("Value")+scale_x_continuous(limits=c(101,500),minor_breaks = seq(150,500,by=50))+theme(text = element_text(size=15))
    g2<-ggplot()+ theme_bw()+geom_line(data=data.frame(x=c(1:500),y=results_BFF_post_param_illus$lambda_values),aes(x,y))+geom_point(data=data.frame(x=c(200),y=min(results_BFF_post_param_illus$lambda_values[101:500])-0.01),shape=4,size=4,aes(x,y))+geom_point(data=data.frame(x=300,y=min(results_BFF_post_param_illus$lambda_values[101:500])-0.01),shape=1,size=4,aes(x,y))+xlab("Observation")+ylim(min(results_BFF_post_param_illus$lambda_values[101:500])-0.01,max(results_BFF_post_param_illus$lambda_values[101:500]))+ylab(expression(lambda~Estimate))+scale_x_continuous(limits=c(101,500),minor_breaks = seq(150,500,by=50))+theme(text = element_text(size=15))
    grid.arrange(g1,g2,ncol=1)

    ## Warning: Removed 100 row(s) containing missing values (geom_path).

    ## Warning: Removed 100 row(s) containing missing values (geom_path).

![](Usage_files/figure-markdown_strict/unnamed-chunk-2-1.png)

Sequential Example
==================

Will now apply sequentially to larger data example.

Data Generation:
----------------

    n<-250000

    # Implementing same data as thesis large scale examples:

    mean_signal_vals <- 0
    mean_signal <- c()


    bkp <- sort(c(sample(c(2:(n-2)),n/200)))
    if((n/200)>1){
      possible_values <- setdiff(c(2:n-2),(rep(bkp,each=31)+rep(c(0:30),times=length(bkp))))
      while(any(diff(c(0,bkp))<30)){
        not_valid <- which(diff(c(0,bkp))<30)
        bkp<-bkp[-not_valid]
        bkp<-sort(c(bkp,sample(possible_values,n/200-length(bkp))))
        possible_values <- setdiff(possible_values,(rep(bkp,each=31)+rep(c(0:30),times=length(bkp))))
      }
    }
    lg <- diff(c(0, bkp, n)) ### rep(n/K, K)
    trend_points <- sample(which((lg[-1]>81) &  (lg[-1]<300)),n/1000)+1

    for(k in 1:length(lg)){
      if(k>1){
        mean_signal_vals <- c(mean_signal_vals,mean_signal_vals[k-1]+(sample(c(-1,1),1)*runif(1,3,6)))
      }
      if((k) %in% trend_points){
        # Randomly choose trend length
        trend_length <- sample(c(50:(lg[k]-30)),1)
        # Randomly choose gradient
        chosen_grad <- sample(c(0.05,0.06,0.07,0.08),1)
        mean_signal <- c(mean_signal,seq(from=(mean_signal_vals[k-1]),by=(sample(c(-1,1),1)*chosen_grad),length.out=trend_length))
        mean_signal <- c(mean_signal,rep(tail(mean_signal,1),lg[k]-trend_length))
        mean_signal_vals[k] <- tail(mean_signal,1)
      }
      else{
        mean_signal <- c(mean_signal,rep(mean_signal_vals[k],lg[k]))
      }
      
    }
    trend_points <- trend_points-1
    bkp_wo_trend <- bkp[-trend_points]
    noisy_signal <- rnorm(length(mean_signal),mean_signal,sqrt(1))

Apply initial model:
--------------------

    # Apply initial model first to 10,000 of data:

    initial_BFF_results <- BFF_gaussian_seq_initial(noisy_signal[1:10000],threshold_val = c(0.005),burnin=100,grace_period=20,sw=2000,alpha=39,beta=1.8,post_p = TRUE,pred_p = TRUE)

Update Model:
-------------

    current_BFF_results <- initial_BFF_results 
    for (i in c(10001:250000)) {
      new_BFF_results <- BFF_gaussian_seq_update(noisy_signal[i],noisy_signal[i-1],i,threshold_val = c(0.005),grace_period=20,sw=2000,alpha=39,beta=1.8,post_p = TRUE,pred_p = TRUE,
                                                 current_BFF_results$current_mu,
                                                 current_BFF_results$current_sigma2,
                                                 current_BFF_results$N_values_old_old,
                                                 current_BFF_results$D_values_old_old,
                                                 current_BFF_results$M_values_old_old,
                                                 current_BFF_results$N_values_old,
                                                 current_BFF_results$D_values_old,
                                                 current_BFF_results$M_values_old,
                                                 current_BFF_results$mu_0_prior_value_old,
                                                 current_BFF_results$sigma2_0_prior_value_old,
                                                 current_BFF_results$sigma2_alpha_0_prior_value_old,
                                                 current_BFF_results$sigma2_beta_0_prior_value_old,
                                                 current_BFF_results$mu_0_prior_value,
                                                 current_BFF_results$sigma2_alpha_0_prior_value,
                                                 current_BFF_results$sigma2_beta_0_prior_value,
                                                 current_BFF_results$anomalous_mu_threshold_sw,
                                                 current_BFF_results$anomalous_lambda_threshold_sw,
                                                 current_BFF_results$p_values_change_mu,
                                                 current_BFF_results$p_values_change_lambda,
                                                 current_BFF_results$anomalous_predictive_threshold_cal,
                                                 current_BFF_results$p_values_predictive)
      # Do something here with estimates if desired such as store them
      current_BFF_results <- new_BFF_results 
    }


    # Performance measures
    bff_results_post_lambda_cal_fast_0.005 <- c(unlist(performance_large_changepoint(current_BFF_results$anomalous_lambda_threshold_sw[[1]][-1]-1,bkp_wo_trend[bkp_wo_trend>2000],length(noisy_signal),20)))
      
    bff_results_predictive_cal_fast_0.005 <- c(unlist(performance_large_changepoint(current_BFF_results$anomalous_predictive_threshold_cal[[1]][-1]-1,bkp_wo_trend[bkp_wo_trend>2000],length(noisy_signal),20)))

    results_df <- rbind(bff_results_post_lambda_cal_fast_0.005,bff_results_predictive_cal_fast_0.005)
    colnames(results_df) <- c("F1","ARL0","ARL1","CCD","DNF","FP","Positives")
    rownames(results_df) <- c("BFF Post Lambda","BFF Pred")
    results_df

    ##                        F1      ARL0      ARL1       CCD       DNF  FP Positives
    ## BFF Post Lambda 0.8350056 4818.4314 1.0280374 0.7542800 0.9350811  52       801
    ## BFF Pred        0.7348148  862.0174 0.1209677 0.7492447 0.7209302 288      1032
