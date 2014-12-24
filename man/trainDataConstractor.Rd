\name{trainDataConstractor}
\alias{trainDataConstractor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
working on...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trainDataConstractor(peptide.unique.fractions.pId, peptide.unique.fractions, numfraction, protein.weight, protein.length, protein.tryptic, protein.pI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peptide.unique.fractions.pId}{
%%     ~~Describe \code{peptide.unique.fractions.pId} here~~
}
  \item{peptide.unique.fractions}{
%%     ~~Describe \code{peptide.unique.fractions} here~~
}
  \item{numfraction}{
%%     ~~Describe \code{numfraction} here~~
}
  \item{protein.weight}{
%%     ~~Describe \code{protein.weight} here~~
}
  \item{protein.length}{
%%     ~~Describe \code{protein.length} here~~
}
  \item{protein.tryptic}{
%%     ~~Describe \code{protein.tryptic} here~~
}
  \item{protein.pI}{
%%     ~~Describe \code{protein.pI} here~~
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
function (peptide.unique.fractions.pId, peptide.unique.fractions, 
    numfraction, protein.weight, protein.length, protein.tryptic, 
    protein.pI) 
{
    fraction.weight <- list()
    fraction.length <- list()
    fraction.tryptic <- list()
    fraction.pI <- list()
    length(fraction.weight) <- numfraction
    length(fraction.length) <- numfraction
    length(fraction.tryptic) <- numfraction
    length(fraction.pI) <- numfraction
    stp <- 100
    recorder <- 0
    for (i in 1:length(peptide.unique.fractions.pId)) {
        recorder <- recorder + 1
        pId <- peptide.unique.fractions.pId[i]
        sId <- NULL
        opt.sId <- as.numeric(which.max(peptide.unique.fractions[i, 
            ]))
        r <- sum(!is.na(peptide.unique.fractions[i, ]))
        if (r < 3) {
            sId <- opt.sId
        }
        else if (max(peptide.unique.fractions[i, ], na.rm = T) != 
            median(as.numeric(peptide.unique.fractions[i, ]), 
                na.rm = T)) {
            sId <- opt.sId
        }
        else {
            sId <- opt.sId + 1
        }
        tmp <- fraction.weight[[sId]]
        tmp <- c(tmp, protein.weight[pId])
        fraction.weight[[sId]] <- tmp
        tmp <- fraction.length[[sId]]
        tmp <- c(tmp, protein.length[pId])
        fraction.length[[sId]] <- tmp
        tmp <- fraction.tryptic[[sId]]
        tmp <- c(tmp, protein.tryptic[pId])
        fraction.tryptic[[sId]] <- tmp
        tmp <- fraction.pI[[sId]]
        tmp <- c(tmp, protein.pI[pId])
        fraction.pI[[sId]] <- tmp
        if (recorder == stp) {
            print(i)
            recorder <- 0
        }
    }
    fraction <- list()
    for (i in 1:numfraction) {
        fraction[[i]] <- summary(as.factor(names(fraction.weight[[i]])), 
            maxsum = Inf)
    }
    pId <- unique(peptide.unique.fractions.pId)
    for (i in 1:length(pId)) {
        assignment <- c()
        for (j in 1:numfraction) {
            assignment <- c(assignment, fraction[[j]][pId[i]])
        }
        count <- sum(!is.na(assignment))
        if (count > 1) {
            best.sId <- which.max(assignment)
            del.idx <- which(!is.na(assignment))
            for (j in 1:length(del.idx)) {
                idx <- del.idx[j]
                if (idx != best.sId) {
                  del <- which(names(fraction.weight[[idx]]) == 
                    pId[i])
                  fraction.weight[[idx]] <- fraction.weight[[idx]][-c(del)]
                  fraction.length[[idx]] <- fraction.length[[idx]][-c(del)]
                  fraction.tryptic[[idx]] <- fraction.tryptic[[idx]][-c(del)]
                  fraction.pI[[idx]] <- fraction.pI[[idx]][-c(del)]
                }
            }
        }
    }
    train.raw <- list()
    train.raw$fraction.weight <- fraction.weight
    train.raw$fraction.length <- fraction.length
    train.raw$fraction.tryptic <- fraction.tryptic
    train.raw$fraction.pI <- fraction.pI
    return(train.raw)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
