\name{proteinBarPlot}
\alias{proteinBarPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
working on...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
proteinBarPlot(protein.original, protein.fraction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein.original}{
%%     ~~Describe \code{protein.original} here~~
}
  \item{protein.fraction}{
%%     ~~Describe \code{protein.fraction} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (protein.original, protein.fraction) 
{
    dev.new(width = 3.5, height = 5)
    barplot(c(length(protein.original), length(protein.fraction)), 
        names = c("Original", "Re-Fraction"), ylab = "Number of Deterministic Protein Identifications", 
        col = c("darkblue", "red"), border = NA)
    text(0.7, length(protein.original)/2, length(protein.original), 
        cex = 1, col = "white")
    text(1.9, length(protein.fraction)/2, length(protein.fraction), 
        cex = 1, col = "white")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
