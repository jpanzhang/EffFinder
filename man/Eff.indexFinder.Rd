% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Eff.indexFinder.R
\name{Eff.indexFinder}
\alias{Eff.indexFinder}
\title{Eff.indexFinder, to detect the Normalization efficiency of any normalization factor (single reference gene or possible multi-RG combination)
according to its serial number from Index file}
\usage{
Eff.indexFinder(NF_index)
}
\arguments{
\item{NF_index}{The serial number of any NF, you can find from the Index file}
}
\value{
NEdata_order The normalization efficiency (NE) values of normalization factors you want to detect
}
\description{
Eff.indexFinder, to detect the Normalization efficiency of any normalization factor (single reference gene or possible multi-RG combination)
according to its serial number from Index file
}
\examples{
Eff.indexFinder(NF_index=c(1:12))
Eff.indexFinder(NF_index=12:20)

}
