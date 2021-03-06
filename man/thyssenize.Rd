% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thyssenize.R
\name{thyssenize}
\alias{thyssenize}
\title{Calculate area weighted values for regions based on station time series data using Thyssen ploygons}
\usage{
thyssenize(station_data, location, shape, shape_id, dates, cores = NULL)
}
\arguments{
\item{station_data}{Tibble providing the station time series data of a measured variable. The data must be organized in three columns where the first column provides the date, the second a unique station ID, and the third column the variable value. See examples for example table.}

\item{location}{sf object that only holds the corresponding station identifiers and the station location.}

\item{shape}{sf object providing polygons for which the area weighted means should be calculated.}

\item{shape_id}{Polygon attribute that should be used as identifier for value averaging}

\item{dates}{Vector of dates for which weighted averages should be calculated.}

\item{cores}{Parallel computing is implemented to reduce computation time.}
}
\value{
Returns a tibble holding dates, shape_id, and averaged values
}
\description{
Calculate area weighted values for regions based on station time series data using Thyssen ploygons
}
