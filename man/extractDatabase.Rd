\name{extractDatabase}
\alias{extractDatabase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   protein database extraction
}
\description{
   Extracting protein properties from a given fasta protein database. This information is required by ReFraction.
}
\usage{
extractDatabase(database.fasta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{database.fasta}{
   A protein database used for peptide identification in fasta format. You need to provide the complete path to this file.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{protein.database}{
 A matrix that containing protein properties extracted from protein sequences. This matrix is essential for ReFraction learning.
 }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
