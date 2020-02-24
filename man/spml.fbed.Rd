\name{Forward Backward Early Dropping selection for circular data using the SPML regression}
\alias{spml.fbed}
\title{
Forward Backward Early Dropping selection for circular data using the SPML regression
}

\description{
Forward Backward Early Dropping selection for circular data using the SPML regression.
}

\usage{
spml.fbed(y, x, alpha = 0.05, K = 0, backward = FALSE,
         parallel = FALSE, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The response variable, a numeric vector expressed in rads.
}
\item{x}{
A matrix with continuous independent variables.
}
\item{alpha}{
The significance threshold value for assessing p-values. Default value is 0.05.
}
%\item{ini}{
%If you already have the test statistics and the p-values of the univariate associations
%(the first step of FBED) supply them as a list with the names "stat" and "pvalue" respectively.
%If you have used the EBIc this list contains the eBIC of the univariate associations.
%Note, that the "gam" argument must be the same though.
%}
\item{K}{
How many times should the process be repeated? The default value is 0.
}
\item{backward}{
After the Forward Early Dropping phase, the algorithm proceeds witha the usual Backward Selection phase.
The default value is set to TRUE. It is advised to perform this step as maybe some variables are false positives,
they were wrongly selected. This is rather experimental now and there could be some mistakes in the indices of the
selected variables. \bold{Do not use it for now}.
}
\item{parallel}{
If you want the algorithm to run in parallel set this TRUE.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations Newton-Raphson will perform.
}
}

\details{
The algorithm is a variation of the usual forward selection. At every step, the most
significant variable enters the selected variables set. In addition, only the significant
variables stay and are further examined. The non signifcant ones are dropped. This goes
until no variable can enter the set. The user has the option to re-do this step 1 or more times
(the argument K). In the end, a backward selection is performed to remove falsely selected variables.
Note that you may have specified, for example, K=10, but the maximum value FBED used can be 4 for example.
}

\value{
If K is a single number a list including:
%\item{univ}{
%If you have used the log-likelihood ratio test this list contains the test statistics and the associated p-values
%of the univariate associations tests. If you have used the EBIc this list contains the eBIC of the univariate associations.
Note, that the "gam" argument must be the same though.
%}
\item{res}{
A matrix with the selected variables and their test statistic.
}
\item{info}{
A matrix with the number of variables and the number of tests performed
(or models fitted) at each round (value of K). This refers to the
forward phase only.
}
\item{runtime}{
The runtime required.
}
%\item{back.rem}{
%The variables removed in the backward phase.
%}
%\item{back.n.tests}{
%The number of models fitted in the backward phase.
%}
}

\references{
Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. Journal of Machine Learning Research, 20(8): 1-39.

Tsagis M. (2018). Guide on performing feature selection with the R package MXM.
\url{ https://f1000researchdata.s3.amazonaws.com/manuscripts/22822/f8d7dae0-c9d3-4ce8-afa1-e552c8cf3347_16216_-_michail_tsagris_v2.pdf?doi=10.12688/f1000research.16216.2&numberOfBrowsableCollections=20&numberOfBrowsableInstitutionalCollections=5&numberOfBrowsableGateways=22 }

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data.
Journal of the American Statistical Association, 93(443): 1068-1077.
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
}

\seealso{
\code{ \link{spml.reg}, \link{spml.regs}, \link{spml.mle} }
}

\examples{
x <- matrix( runif(100 * 50, 1, 100), ncol = 50 )
y <- runif(100)
a <- spml.fbed(y, x)
}

