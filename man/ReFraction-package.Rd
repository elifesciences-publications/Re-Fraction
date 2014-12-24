\name{ReFraction-package}
\alias{ReFraction-package}
\alias{ReFraction}
\docType{package}
\title{
Re-Fraction
}
\description{
A machine learning approach for deterministic identification of protein homologs and splice variants in large-scale MS-based proteomics. 
This algorithm relies on the fractionation performed on protein level by SDS-PAGE for example and use such information to reassign peptides that were assigned non-uniquely to multiple proteins. 
Please refer to reference [1] for details.
}
\details{
\tabular{ll}{
Package: \tab ReFraction\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2012-04-04\cr
License: \tab GPL-2\cr
}

}
\author{
Pengyi Yang

Maintainer: Pengyi Yang <yangpy7@gmail.com>
}
\references{
[1] Pengyi Yang, Sean J. Humphrey, Daniel J. Fazakerley, Matthew J. Prior, Guang Yang, David E. James, and Jean Yee-Hwa Yang, 
Re-Fraction: a machine learning approach for deterministic identification of protein homologs and splice variants in large-scale MS-based proteomics, 
Journal of Proteome Research, http://dx.doi.org/10.1021/pr300072j
}

\keyword{ package }
\seealso{

}
\examples{
# Use following for details on how to apply ReFraction with example data.
?applyReFraction 

# You need to provide protein database and peptide table. You can prepare your data according to the example data. 
# There is also a function for extracting protein properties from a fasta database file.
?extractDatabase

}
