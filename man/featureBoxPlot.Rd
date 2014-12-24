\name{featureBoxPlot}
\alias{featureBoxPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
working on...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
featureBoxPlot(fraction, numFraction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fraction}{
%%     ~~Describe \code{fraction} here~~
}
  \item{numFraction}{
%%     ~~Describe \code{numFraction} here~~
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
function (fraction, numFraction) 
{
    dev.new(width = 9, height = 7)
    par(mfrow = c(2, 2), mai = c(0.45, 0.65, 0.45, 0.25), cex = 0.75)
    boxplot(fraction$fraction.weight, col = rainbow(numFraction + 
        1), names = paste("F", 1:numFraction, sep = ""), xlab = "Fraction", 
        ylab = "Mass (kDa)", main = "(a)")
    boxplot(fraction$fraction.length, col = rainbow(numFraction + 
        1), names = paste("F", 1:numFraction, sep = ""), ylab = "Length (aa)", 
        main = "(b)")
    boxplot(lapply(fraction$fraction.tryptic, log2), col = rainbow(numFraction + 
        1), names = paste("F", 1:numFraction, sep = ""), ylab = "Log2(Number of Peptides)", 
        main = "(c)")
    boxplot(fraction$fraction.pI, col = rainbow(numFraction + 
        1), names = paste("F", 1:numFraction, sep = ""), ylab = "Isoelectric Point", 
        main = "(d)")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
