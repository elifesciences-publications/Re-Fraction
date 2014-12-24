\name{peptideBarPlot}
\alias{peptideBarPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
working on...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
peptideBarPlot(peptide.assignment.fraction, peptide.assignment.original, prefix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peptide.assignment.fraction}{
%%     ~~Describe \code{peptide.assignment.fraction} here~~
}
  \item{peptide.assignment.original}{
%%     ~~Describe \code{peptide.assignment.original} here~~
}
  \item{prefix}{
%%     ~~Describe \code{prefix} here~~
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
function (peptide.assignment.fraction, peptide.assignment.original, 
    prefix) 
{
    peptide.slice.unique.idx <- lapply(strsplit(peptide.assignment.fraction, 
        ";"), length) == 1
    peptide.original.unique.idx <- lapply(strsplit(peptide.assignment.original, 
        ";"), length) == 1
    unique.count <- matrix(nrow = 2, ncol = 2)
    unique.count[1, 1] <- sum(peptide.original.unique.idx)
    unique.count[1, 2] <- sum(peptide.slice.unique.idx)
    unique.count[2, 1] <- length(peptide.original.unique.idx) - 
        sum(peptide.original.unique.idx)
    unique.count[2, 2] <- length(peptide.slice.unique.idx) - 
        sum(peptide.slice.unique.idx)
    dev.new(width = 3.5, height = 5)
    barplot(unique.count, ylim = c(0, 20000), ylab = "Number of Peptides", 
        names = c("Original", "Re-Fraction"), col = c("darkblue", 
            "red"), legend = c("Unique Peptide", "Shared Peptide"), 
        horiz = F, border = NA)
    text(0.7, unique.count[1, 1]/2, unique.count[1, 1], cex = 1, 
        col = "white")
    text(0.7, unique.count[1, 1] + unique.count[2, 1]/2, unique.count[2, 
        1], cex = 1, col = "white")
    text(1.9, unique.count[1, 2]/2, unique.count[1, 2], cex = 1, 
        col = "white")
    text(1.9, unique.count[1, 2] + unique.count[2, 2]/2, unique.count[2, 
        2], cex = 1, col = "white")
    dev.new(width = 3.5, height = 5)
    peptide.assignment.ratio <- c(sum(peptide.original.unique.idx)/length(peptide.original.unique.idx), 
        sum(peptide.slice.unique.idx)/length(peptide.slice.unique.idx))
    peptide.assignment.ratio <- round(peptide.assignment.ratio * 
        100)
    barplot(peptide.assignment.ratio, names = c("Original", "Re-Fraction"), 
        ylab = "Percentage of Unique Peptide (\%)", col = c("darkblue", 
            "red"), border = NA)
    text(0.7, peptide.assignment.ratio[1]/2, paste(round(peptide.assignment.ratio[1]), 
        "\%", sep = ""), cex = 1, col = "white")
    text(1.9, peptide.assignment.ratio[2]/2, paste(round(peptide.assignment.ratio[2]), 
        "\%", sep = ""), cex = 1, col = "white")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
