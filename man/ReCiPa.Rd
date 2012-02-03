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
Version: \tab 2.0\cr
Date: \tab 2012-02-03\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Juan C. Vivar

Maintainer: Juan C. Vivar <jvivar@nccu.edu>
}
\examples{ \dontrun{
ReCiPa("original_pathway_db.txt", "new_name", max_overlapping)
}
}
\keyword{pathway}

