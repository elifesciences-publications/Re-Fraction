\name{proteinClassification}
\alias{proteinClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
working on...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
proteinClassification(peptide, svm_model, protein.weight, protein.length, protein.tryptic, protein.pI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peptide}{
%%     ~~Describe \code{peptide} here~~
}
  \item{svm_model}{
%%     ~~Describe \code{svm_model} here~~
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
function (peptide, svm_model, protein.weight, protein.length, 
    protein.tryptic, protein.pI) 
{
    all.proteins <- unique(unlist(strsplit(as.character(peptide[, 
        "Proteins"]), ";")))
    all.proteins <- all.proteins[-grep("CON__", all.proteins)]
    all.proteins.df <- data.frame(cbind(protein.weight[all.proteins], 
        protein.length[all.proteins], protein.tryptic[all.proteins], 
        protein.pI[all.proteins]))
    rownames(all.proteins.df) <- all.proteins
    protein.regressed <- as.numeric(predict(svm_model, all.proteins.df))
    names(protein.regressed) <- all.proteins
    result <- list()
    result$regressed <- protein.regressed
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
