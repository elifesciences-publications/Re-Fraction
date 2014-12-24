\name{trainSVM}
\alias{trainSVM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
working on...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trainSVM(train, numfraction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train}{
%%     ~~Describe \code{train} here~~
}
  \item{numfraction}{
%%     ~~Describe \code{numfraction} here~~
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
function (train, numfraction) 
{
    fraction.weight <- train$fraction.weight
    fraction.length <- train$fraction.length
    fraction.tryptic <- train$fraction.tryptic
    fraction.pI <- train$fraction.pI
    train.list <- list()
    for (i in 1:numfraction) {
        train.list[[i]] <- cbind(fraction.weight[[i]], fraction.length[[i]], 
            fraction.tryptic[[i]], fraction.pI[[i]])
    }
    train.mat <- train.list[[1]]
    for (i in 2:numfraction) {
        train.mat <- rbind(train.mat, train.list[[i]])
    }
    rownames(train.mat) <- NULL
    counts.fraction <- unlist(lapply(train.list, nrow))
    class.numeric <- rep(1:numfraction, counts.fraction)
    library(e1071)
    library(caret)
    set.seed(1)
    fold <- createFolds(class.numeric, 10)
    gm.opt <- 1
    acc.opt <- 0
    slice.acc.opt <- c()
    slice.sen.opt <- c()
    slice.spe.opt <- c()
    for (gm in seq(1, 50, by = 5)) {
        overall.acc <- c()
        slice.acc <- c()
        slice.sen <- c()
        slice.spe <- c()
        for (f in 1:length(fold)) {
            train.data <- train.mat[-c(fold[[f]]), ]
            train.cls <- class.numeric[-c(fold[[f]])]
            test.data <- train.mat[c(fold[[f]]), ]
            test.cls <- class.numeric[c(fold[[f]])]
            model <- svm(train.data, train.cls, gamma = gm)
            preds <- predict(model, test.data)
            count <- 0
            TP <- rep(0, numfraction)
            TN <- rep(0, numfraction)
            names(TP) <- c("TP1", "TP2", "TP3", "TP4", "TP5", 
                "TP6", "TP7", "TP8", "TP9", "TP10")
            FP <- rep(0, numfraction)
            FN <- rep(0, numfraction)
            names(FN) <- c("FN1", "FN2", "FN3", "FN4", "FN5", 
                "FN6", "FN7", "FN8", "FN9", "FN10")
            for (i in 1:length(preds)) {
                if (round(preds[i]) == test.cls[i]) {
                  count <- count + 1
                }
                for (j in 1:numfraction) {
                  if (round(preds[i]) == j) {
                    if (round(preds[i]) == test.cls[i]) {
                      TP[j] <- TP[j] + 1
                    }
                    else {
                      FP[j] <- FP[j] + 1
                    }
                  }
                  else {
                    if (test.cls[i] != j) {
                      TN[j] <- TN[j] + 1
                    }
                    else {
                      FN[j] <- FN[j] + 1
                    }
                  }
                }
            }
            overall.acc <- c(overall.acc, (count/length(preds)))
            for (j in 1:numfraction) {
                acc <- (TP[j] + TN[j])/(FP[j] + FN[j] + TP[j] + 
                  TN[j])
                sen <- TP[j]/(TP[j] + FN[j])
                spe <- 1 - FP[j]/(FP[j] + TN[j])
                slice.acc <- c(slice.acc, acc)
                slice.sen <- c(slice.sen, sen)
                slice.spe <- c(slice.spe, spe)
            }
        }
        dim(slice.acc) <- c(10, numfraction)
        grand.mean <- mean(apply(slice.acc, 1, mean))
        if (acc.opt < grand.mean) {
            acc.opt <- grand.mean
            gm.opt <- gm
            slice.acc.opt <- slice.acc
            slice.sen.opt <- slice.sen
            slice.spe.opt <- slice.spe
            print("current best gamma:")
            print(gm.opt)
        }
    }
    print("peformance by fraction:")
    print("accuracy:")
    dim(slice.acc.opt) <- c(10, numfraction)
    print(apply(slice.acc.opt, 1, mean))
    print("sensitivity:")
    dim(slice.sen.opt) <- c(10, numfraction)
    print(apply(slice.sen.opt, 1, mean))
    print("specificity:")
    dim(slice.spe.opt) <- c(10, numfraction)
    print(apply(slice.spe.opt, 1, mean))
    svm_model <- svm(train.mat, class.numeric, gamma = gm.opt)
    result <- list()
    result$train.mat <- train.mat
    result$class.numeric <- class.numeric
    result$svm_model <- svm_model
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
