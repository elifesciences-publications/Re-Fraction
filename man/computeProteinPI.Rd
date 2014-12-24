\name{computeProteinPI}
\alias{computeProteinPI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
protein pI computation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeProteinPI(keyAA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keyAA}{
%%     ~~Describe \code{keyAA} here~~
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
function (keyAA) 
{
    AspNumber <- keyAA[1]
    GluNumber <- keyAA[2]
    CysNumber <- keyAA[3]
    TyrNumber <- keyAA[4]
    HisNumber <- keyAA[5]
    LysNumber <- keyAA[6]
    ArgNumber <- keyAA[7]
    NQ <- 0
    QN1 <- 0
    QN2 <- 0
    QN3 <- 0
    QN4 <- 0
    QN5 <- 0
    QP1 <- 0
    QP2 <- 0
    QP3 <- 0
    QP4 <- 0
    pH <- 6.5
    pHprev <- 0
    pHnext <- 14
    X <- 0
    E <- 0.01
    temp <- 0
    tabl <- 650
    while (TRUE) {
        QN1 <- tab[tabl, 1]
        QN2 <- -AspNumber/tab[tabl, 2]
        QN3 <- -GluNumber/tab[tabl, 3]
        QN4 <- -CysNumber/tab[tabl, 4]
        QN5 <- -TyrNumber/tab[tabl, 5]
        QP1 <- HisNumber/tab[tabl, 6]
        QP2 <- tab[tabl, 7]
        QP3 <- LysNumber/tab[tabl, 8]
        QP4 <- ArgNumber/tab[tabl, 9]
        NQ <- QN1 + QN2 + QN3 + QN4 + QN5 + QP1 + QP2 + QP3 + 
            QP4
        if (NQ < 0) {
            temp <- pH
            pH <- pH - ((pH - pHprev)/2)
            pHnext <- temp
        }
        else {
            temp <- pH
            pH <- pH + ((pHnext - pH)/2)
            pHprev <- temp
        }
        if ((pH - pHprev < E) & (pHnext - pH < E)) {
            break
        }
        n <- pH * 100
        obc <- n
        ulamek <- n - obc
        if (ulamek >= 0.5) {
            obc <- obc + 1
        }
        tabl <- obc
        pH <- obc
        pH <- pH/100
    }
    return(pH)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
