\name{determineProteinID}
\alias{determineProteinID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
working on...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
determineProteinID(peptide, peptide.slice.assignment, peptide.original.assignment, fdr.cutoff = 0.01, fdr.stp = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peptide}{
%%     ~~Describe \code{peptide} here~~
}
  \item{peptide.slice.assignment}{
%%     ~~Describe \code{peptide.slice.assignment} here~~
}
  \item{peptide.original.assignment}{
%%     ~~Describe \code{peptide.original.assignment} here~~
}
  \item{fdr.cutoff}{
%%     ~~Describe \code{fdr.cutoff} here~~
}
  \item{fdr.stp}{
%%     ~~Describe \code{fdr.stp} here~~
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
function (peptide, peptide.slice.assignment, peptide.original.assignment, 
    fdr.cutoff = 0.01, fdr.stp = 20) 
{
    peptide.slice.unique.idx <- lapply(strsplit(peptide.slice.assignment, 
        ";"), length) == 1
    peptide.original.unique.idx <- lapply(strsplit(peptide.original.assignment, 
        ";"), length) == 1
    peptide.slice.unique <- peptide.slice.assignment[peptide.slice.unique.idx]
    peptide.original.unique <- peptide.original.assignment[peptide.original.unique.idx]
    protein.slice.unique <- unique(peptide.slice.unique)
    protein.original.unique <- unique(peptide.original.unique)
    protein.slice.peptideCount <- summary(as.factor(peptide.slice.unique), 
        maxsum = Inf)
    protein.original.peptideCount <- summary(as.factor(peptide.original.unique), 
        maxsum = Inf)
    protein.slice.multiPeptide <- names(protein.slice.peptideCount)[(protein.slice.peptideCount > 
        1)]
    protein.original.multiPeptide <- names(protein.original.peptideCount)[(protein.original.peptideCount > 
        1)]
    peptide.slice.pep <- cbind(peptide[peptide.slice.unique.idx, 
        "PEP"], peptide.slice.unique)
    rownames(peptide.slice.pep) <- names(peptide.slice.unique)
    colnames(peptide.slice.pep) <- c("PEP", "pId")
    peptide.original.pep <- cbind(peptide[peptide.original.unique.idx, 
        "PEP"], peptide.original.unique)
    rownames(peptide.original.pep) <- names(peptide.original.unique)
    colnames(peptide.original.pep) <- c("PEP", "pId")
    peptide.slice.pep.sort <- peptide.slice.pep[order(peptide.slice.pep[, 
        "PEP"]), ]
    peptide.original.pep.sort <- peptide.original.pep[order(peptide.original.pep[, 
        "PEP"]), ]
    protein.slice.PEP <- c()
    for (i in 1:length(protein.slice.multiPeptide)) {
        prot <- peptide.slice.pep.sort[which(peptide.slice.pep.sort[, 
            "pId"] == protein.slice.multiPeptide[i]), "PEP"]
        protein.slice.PEP <- c(protein.slice.PEP, prod(as.numeric(prot)))
    }
    names(protein.slice.PEP) <- protein.slice.multiPeptide
    protein.slice.sort <- protein.slice.PEP[order(protein.slice.PEP)]
    protein.original.PEP <- c()
    for (i in 1:length(protein.original.multiPeptide)) {
        prot <- peptide.original.pep.sort[which(peptide.original.pep.sort[, 
            "pId"] == protein.original.multiPeptide[i]), "PEP"]
        protein.original.PEP <- c(protein.original.PEP, prod(as.numeric(prot)))
    }
    names(protein.original.PEP) <- protein.original.multiPeptide
    protein.original.sort <- protein.original.PEP[order(protein.original.PEP)]
    index <- 0
    count <- 0
    stp <- fdr.stp
    decoy <- 0
    for (i in 1:length(protein.slice.sort)) {
        count <- count + 1
        if (regexpr("REV__", names(protein.slice.sort)[i]) > 
            0) {
            decoy <- decoy + 1
        }
        if (stp == count) {
            count <- 0
            fdr <- 2 * decoy/i
            index <- i
            if (fdr > fdr.cutoff) {
                break()
            }
        }
    }
    protein.slice.filtered <- protein.slice.sort[1:index]
    index <- 0
    count <- 0
    stp <- fdr.stp
    decoy <- 0
    for (i in 1:length(protein.original.sort)) {
        count <- count + 1
        if (regexpr("REV__", names(protein.original.sort)[i]) > 
            0) {
            decoy <- decoy + 1
        }
        if (stp == count) {
            count <- 0
            fdr <- 2 * decoy/i
            index <- i
            if (fdr > fdr.cutoff) {
                break()
            }
        }
    }
    protein.original.filtered <- protein.original.sort[1:index]
    protein.slice.filtered <- protein.slice.filtered[grep("REV", 
        invert = T, names(protein.slice.filtered))]
    protein.slice.filtered <- protein.slice.filtered[grep("CON", 
        invert = T, names(protein.slice.filtered))]
    protein.original.filtered <- protein.original.filtered[grep("REV", 
        invert = T, names(protein.original.filtered))]
    protein.original.filtered <- protein.original.filtered[grep("CON", 
        invert = T, names(protein.original.filtered))]
    result <- list()
    result$fraction <- protein.slice.filtered
    result$original <- protein.original.filtered
    result$fraction.pepSort <- peptide.slice.pep.sort
    result$original.pepSort <- peptide.original.pep.sort
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
