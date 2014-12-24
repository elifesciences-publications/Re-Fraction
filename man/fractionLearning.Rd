\name{fractionLearning}
\alias{fractionLearning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
working on...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fractionLearning(peptide, numSlice, proteinRegressed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peptide}{
%%     ~~Describe \code{peptide} here~~
}
  \item{numSlice}{
%%     ~~Describe \code{numSlice} here~~
}
  \item{proteinRegressed}{
%%     ~~Describe \code{proteinRegressed} here~~
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
function (peptide, numSlice, proteinRegressed) 
{
    peptide.slices <- peptide[, grep("^Slice.\\d", colnames(peptide))]
    peptide.original.assignment <- as.character(peptide[, "Proteins"])
    names(peptide.original.assignment) <- peptide[, "id"]
    peptide.assignment.range1 <- list()
    peptide.unique.range1 <- c()
    peptide.predict <- list()
    count <- 0
    stp <- 100
    recorder <- 0
    for (i in 1:nrow(peptide)) {
        recorder <- recorder + 1
        prots <- as.character(peptide[i, "Proteins"])
        pIds <- unlist(strsplit(prots, ";"))
        selected.pIds.range1 <- c()
        for (j in 1:length(pIds)) {
            if (regexpr("CON__", pIds[j]) > 0) {
                selected.pIds.range1 <- c(selected.pIds.range1, 
                  pIds[j])
            }
            else {
                count <- count + 1
                preds <- proteinRegressed[pIds[j]]
                peptide.predict[[count]] <- paste(i, pIds[j], 
                  preds, sep = ";")
                preds <- round(preds)
                range1 <- preds
                if (range1 < 1) {
                  range1 <- 1
                }
                if (range1 > numSlice) {
                  range1 <- numSlice
                }
                if (sum(peptide.slices[i, range1], na.rm = T) > 
                  0) {
                  selected.pIds.range1 <- c(selected.pIds.range1, 
                    pIds[j])
                }
            }
        }
        if (length(selected.pIds.range1) > 0) {
            peptide.assignment.range1[i] <- paste(selected.pIds.range1, 
                collapse = ";")
            if (length(selected.pIds.range1) == 1) {
                peptide.unique.range1 <- c(peptide.unique.range1, 
                  selected.pIds.range1)
            }
        }
        else {
            peptide.assignment.range1[i] <- paste(pIds, collapse = ";")
        }
        if (recorder == stp) {
            print(i)
            recorder <- 0
        }
    }
    peptide.assignment.range1 <- unlist(peptide.assignment.range1)
    names(peptide.assignment.range1) <- peptide[, "id"]
    result <- list()
    result$predict <- peptide.predict
    result$original <- peptide.original.assignment
    result$range1 <- peptide.assignment.range1
    result$extra.unique.range1 <- peptide.unique.range1
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
