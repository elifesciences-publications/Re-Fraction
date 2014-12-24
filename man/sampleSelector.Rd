\name{sampleSelector}
\alias{sampleSelector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
working on...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleSelector(peptide, slice.idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peptide}{
%%     ~~Describe \code{peptide} here~~
}
  \item{slice.idx}{
%%     ~~Describe \code{slice.idx} here~~
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
function (peptide, slice.idx) 
{
    peptide.filtered <- peptide[grep("REV", invert = T, peptide[, 
        "Proteins"]), ]
    peptide.filtered <- peptide.filtered[grep("CON", invert = T, 
        peptide.filtered[, "Proteins"]), ]
    peptide.assignment <- as.character(peptide.filtered[, "Proteins"])
    peptide.unique.idx <- lapply(strsplit(peptide.assignment, 
        ";"), length) == 1
    peptide.unique.pId <- peptide.assignment[peptide.unique.idx]
    names(peptide.unique.pId) <- rownames(peptide.filtered)[peptide.unique.idx]
    peptide.slices <- peptide.filtered[peptide.unique.idx, slice.idx]
    peptide.range.slices <- peptide.slices[(rowSums(!is.na(peptide.slices)) <= 
        3), ]
    peptide.count.slices <- peptide.range.slices[apply(peptide.range.slices, 
        1, max, na.rm = T) > 1, ]
    selected <- c()
    for (i in 1:nrow(peptide.count.slices)) {
        idx <- which(!is.na(peptide.count.slices[i, ]))
        if (length(idx) == 1) {
            selected <- c(selected, rownames(peptide.count.slices)[i])
        }
        else {
            flag <- 0
            for (j in 1:(length(idx) - 1)) {
                if ((idx[j] + 1) != idx[j + 1]) {
                  flag <- 1
                }
            }
            if (flag == 0) {
                selected <- c(selected, rownames(peptide.count.slices)[i])
            }
        }
    }
    peptide.unique.slices.pId <- peptide.unique.pId[selected]
    peptide.unique.slices <- peptide.slices[selected, ]
    train.psm <- list()
    train.psm$selected <- selected
    train.psm$pId <- peptide.unique.slices.pId
    train.psm$fraction <- peptide.unique.slices
    return(train.psm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
