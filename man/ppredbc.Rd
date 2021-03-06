% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_model.R
\name{ppredbc}
\alias{ppredbc}
\title{Calculates the PIT values}
\usage{
ppredbc(obs, mean, sd, lambda)
}
\arguments{
\item{mean}{The mean of the fitted distribution}

\item{sd}{The standard deviation of the fitted distribution}

\item{lambda}{The lambda used in the Box Cox transformation}

\item{obs.bc}{Box Cox transformed observations}
}
\value{
The PIT values normalized to the zero one intervall
}
\description{
Calculates the PIT values
}
