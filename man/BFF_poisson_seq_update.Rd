% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BFF_Sequential.R
\name{BFF_poisson_seq_update}
\alias{BFF_poisson_seq_update}
\title{BFF Change Point Detector Sequential Poisson Update}
\usage{
BFF_poisson_seq_update(
  data,
  i,
  threshold_val = c(0.05),
  grace_period = 20,
  sw = 2000,
  alpha = 39,
  beta = 1.8,
  post_p = TRUE,
  pred_p = TRUE,
  current_gamma,
  current_lambda,
  N_values_old_old,
  D_values_old_old,
  F_values_old_old,
  N_values_old,
  D_values_old,
  F_values_old,
  alpha_0_prior_value_old,
  beta_0_prior_value_old,
  alpha_0_prior_value,
  beta_0_prior_value,
  anomalous_gamma_threshold_sw,
  anomalous_lambda_threshold_sw,
  p_values_change_gamma,
  p_values_change_lambda,
  anomalous_pred_threshold_sw,
  p_values_change_pred
)
}
\arguments{
\item{data}{new data point}

\item{i}{data point number}

\item{threshold_val}{vector of thresholds to use to assess whether new data point is a change point}

\item{grace_period}{period after a change in which changes are not detected}

\item{sw}{sliding window size for p-value calibration}

\item{alpha}{alpha parameter for lambda beta prior}

\item{beta}{beta parameter for lambda beta prior}

\item{post_p}{to calculate predictive posterior p-values}

\item{current_gamma}{estimated gamma from previous time}

\item{current_lambda}{estimated lambda from previous time}

\item{N_values_old_old}{old old N value}

\item{D_values_old_old}{old old D value}

\item{F_values_old_old}{old old F value}

\item{N_values_old}{old N value}

\item{D_values_old}{old D value}

\item{F_values_old}{old F value}

\item{alpha_0_prior_value_old}{old prior alpha for gamma}

\item{beta_0_prior_value_old}{old prior beta for gamma}

\item{alpha_0_prior_value}{prior alpha for gamma}

\item{beta_0_prior_value}{prior beta for gamma}

\item{anomalous_gamma_threshold_sw}{current anomalous gamma posterior p-values list}

\item{anomalous_lambda_threshold_sw}{current anomalous lambda posterior p-values list}

\item{p_values_change_gamma}{uncalibrated gamma posterior p-values}

\item{p_values_change_lambda}{uncalibrated lambda posterior p-values}

\item{anomalous_pred_threshold_sw}{current anomalous predictive posterior p-values list}

\item{p_values_change_pred}{uncalibrated predictive posterior p-values}

\item{data_old}{previous data point}
}
\value{
return the change points detected by the algorithm for each of the dection procedures, for the different thresholds.
}
\description{
Detects change points within data stream sequentially using a Bayesian adaptive estimation
procedure and predictive posterior p-values. Additionally estimates parameters
}
