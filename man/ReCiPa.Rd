\name{ReCiPa}
\alias{ReCiPa}
\alias{check.null.dim}
\docType{package}
\title{
Redundancy Control in Pathways databases
}
\description{
Pathways in a database could have many redundancies among them. This package allows the user to set a maximum value for the proportion of these redundancies.
}
\details{
\tabular{ll}{
Package: \tab ReCiPa\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2012-12-21\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Juan C. Vivar

Maintainer: Juan C. Vivar <jvivar@nccu.edu>
}
\examples{ \dontrun{
ReCiPa("original_pathway_db", "new_db_name", max_overlapping, min_overlapping)
}
}
\keyword{pathway}

