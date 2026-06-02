\name{validate_va_data}
\alias{validate_va_data}
\title{Validate visit-aligned data}
\usage{
validate_va_data(obs_data, endpoint = "composite")
}
\arguments{
\item{obs_data}{Object produced by \code{as_va_data()} or \code{simulate_va_trial()}.}
\item{endpoint}{Endpoint type. Only \code{"composite"} is supported.}
}
\value{The validated object, invisibly.}
\description{
Checks required components, endpoint coding, and basic person-interval
consistency for the composite endpoint.
}
