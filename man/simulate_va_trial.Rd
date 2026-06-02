\name{simulate_va_trial}
\alias{simulate_va_trial}
\title{Simulate a visit-aligned trial}
\usage{
simulate_va_trial(n = 1000, scenario = "S1", intervention = NULL, seed = NULL)
}
\arguments{
\item{n}{Number of subjects.}
\item{scenario}{Scenario name.}
\item{intervention}{Optional static intervention, 0 or 1. Leave \code{NULL}
for the observed treatment mechanism.}
\item{seed}{Optional random seed.}
}
\value{A \code{va_data} object.}
\description{
Generates example scheduled-visit trial data for package examples and tests.
}
