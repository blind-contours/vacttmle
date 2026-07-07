\name{as_va_data}
\alias{as_va_data}
\title{Build a visit-aligned analysis object}
\usage{
as_va_data(subject_data, interval_data, visit_times, tau, visit_L = NULL,
  visit_A = NULL, visit_g = NULL, params = list())
}
\arguments{
\item{subject_data}{Subject-level data.}
\item{interval_data}{Person-interval data.}
\item{visit_times}{Scheduled visit times.}
\item{tau}{Analysis horizon.}
\item{visit_L, visit_A, visit_g}{Optional visit matrices.}
\item{params}{Optional metadata list.}
}
\value{A \code{va_data} object.}
\description{
Converts subject-level and person-interval tables into the analysis object used
by the package estimators.
}
