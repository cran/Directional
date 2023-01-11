\name{Cross validation for estimating the classification rate}
\alias{dirda.cv}
\title{
Cross validation for estimating the classification rate
}

\description{
Cross validation for estimating the classification rate.
}

\usage{
dirda.cv(x, ina, folds = NULL, nfolds = 10, k = 2:10, stratified = FALSE,
         type = c("vmf", "iag", "esag", "kent", "knn"),
         seed = NULL, B = 1000, parallel = FALSE)
}

\arguments{
\item{x}{
A matrix with the data in Eulcidean coordinates, i.e. unit vectors. The matrix must have three columns, only spherical data are currently supported.
}
\item{ina}{
A variable indicating the groupings.
}
\item{folds}{
Do you already have a list with the folds? If not, leave this NULL.
}
\item{nfolds}{
How many folds to create?
}
\item{k}{
If you choose to use k-NN, what will be the k values?
}
\item{stratified}{
Should the folds be created in a stratified way? i.e. keeping the distribution of the groups 
similar through all folds?
}
\item{seed}{
If seed is TRUE, the results will always be the same.
}
\item{type}{
The type of classifier to use. The avaliable options are "vmf" (von Mises-Fisher distribution),
"esag" (ESAG distribution), "kent" (Kent distribution), "knn" (k-NN algorithm). You can chose 
any of them or all of them. Note that "esag" and "kent" work only with spherical data.
}
\item{B}{
If you used k-NN, should a bootstrap correction of the bias be applied? If yes, 1000 is a good value.
}
\item{parallel}{
If you want the standard -NN algorithm to take place in parallel set this equal to TRUE.
}
}

\details{
Cross-validation for the estimation of the performance of a classifier.

The estimated performance of the best classifier is overestimated. After the cross-valdiation 
procedure, the predicted values produced by all classifiers are colelcted, from all folds, 
in an \eqn{n \times M} matrix, where n is the number of samples and M the number of all 
classifiers used. We sample rows (predictions) with replacement from P and denote them as 
the in-sample values. The non re-sampled rows are denoted as out-of-sample values. The 
performance of each classifier in the insample rows is calculated and the classifier with 
the optimal performance is selected, followed by the calculation of performance in the 
out-of-sample  values. This process is repeated B times and the average performance is returned. 
The only computational overhead is with the repetitive resampling and calculation of the performance, 
i.e. no model or classifier is fitted nor trained. For more information see Tsamardinos et al. (2018). 
This procedure though takes place only for the k-NN algorithm.

The good thing with the function is that you can run any method you want by supplying the 
folds yourselves using the command \code{\link{makefolds}}. Then suppose you want to run  another method. 
By suppying the same folds you will be able to have comparative results for all methods.
}

\value{
A list including:
\item{perf}{
A vector with the estimated performance of each classifier.
}
\item{best}{
The classifier with the optimal performance.
}
\item{boot.perf}{
The bootstrap bias corrected performance.
}
}

\references{
Tsagris M. and Alenazi A. (2019). Comparison of discriminant analysis methods on the sphere. 
Communications in Statistics: Case Studies, Data Analysis and Applications, 5(4), 467--491.

Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Paine P.J. Preston S.P. & Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular 
Gaussian Distribution. Statistics and Computing, 28(3):689--697.

Morris J. E. & Laycock P. J. (1974). Discriminant analysis of directional data. Biometrika, 61(2): 335-341.

Tsamardinos I., Greasidou E. & Borboudakis G. (2018). Machince Learning, 107(12): 1895-1922.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{esag.da}, \link{vmfda.pred}, \link{dirknn}, \link{knn.reg}
}
}

\examples{
x <- rvmf(300, rnorm(3), 10)
ina <- sample.int(4, 300, replace = TRUE)
dirda.cv(x, ina, B = 1000)
}

\keyword{ Discriminant analysis }
\keyword{ Von Mises-Fisher distributions }
