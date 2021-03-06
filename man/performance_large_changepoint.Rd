% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Background_functions.R
\name{performance_large_changepoint}
\alias{performance_large_changepoint}
\title{Calculate performance of CP procedure on data}
\usage{
performance_large_changepoint(estimated_cp, true_cp, datasize, D_period)
}
\arguments{
\item{estimated_cp}{discovered change points}

\item{true_cp}{true change points}

\item{datasize}{size of data}

\item{D_period}{grace period}
}
\value{
return all performance measures F1, ARL0, ARL1, recall, precision, no positives, FP etc
}
\description{
Calculate performance of CP procedure on data
}
