\name{run_vacttmle_toy_checks}
\alias{run_vacttmle_toy_checks}
\title{Run package toy checks}
\usage{
run_vacttmle_toy_checks(n = 5000, seed = 7001)
}
\arguments{
\item{n}{Analytic toy-simulation sample size.}
\item{seed}{Random seed.}
}
\value{A data.table with estimates, analytic truths, errors, and pass flags.}
\description{
Runs lightweight analytic validation checks for the composite survival target.
}
