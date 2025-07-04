\name{Directional-package}
\alias{Directional-package}
\docType{package}
\title{
This is an R package that provides methods for the statistical analysis of directional data, including massive (very large scale) directional data.
}
\description{
Circular-linear regression, spherical-spherical regression, spherical regression, discriminant analysis, ANOVA for circular and (hyper-)spherical data, tests for eaquality of conentration parameters, maximum likelihood estimation of the parameters of many distributions, random values generation from various distributions, contour plots and many more functions are included.
}
\details{
\tabular{ll}{
Package: \tab Directional\cr
Type: \tab Package\cr
Version: \tab 7.2\cr
Date: \tab 2025-07-02\cr
License: \tab GPL-2\cr
}
}
\author{
Michail Tsagris \email{mtsagris@uoc.gr}, Giorgos Athineou <gioathineou@gmail.com>, Christos Adam \email{pada4m4@gmail.com}, Zehao Yu <zehaoy@email.sc.edu>, Anamul Sajib <sajibstat@du.ac.bd>, Eli Amson <eli.amson1988@gmail.com>, Micah J. Waldstein <micah@waldste.in> and Panagiotis Papastamoulis <papastamoulis@aueb.gr>.
}

\section{Maintainers}{
Michail Tsagris \email{mtsagris@uoc.gr}.
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley and Sons.
}

\note{
\bold{Acknowledgments:}

Professor Andy Wood and Dr Simon Preston from the university of Nottingham are highly appreciated for being my supervisors during my post-doc in directional data analysis.

Dr Georgios Pappas (former postDoc at the university of Nottingham) helped me construct the contour plots of the von Mises-Fisher and the Kent distribution.

Dr Christopher Fallaize and Dr Theo Kypraios from the university of Nottingham have provided a function for simulating from the Bingham distribution using rejection sampling.
So any questions regarding this function should be addressed to them.

Dr Kwang-Rae Kim (post-doc at the university of Nottingham) answered some of my questions.

Giorgos Borboudakis (PhD student at the university of Crete) pointed out to me a not so clear message in the algorithm of generating random values from the von Mises-Fisher distribution.

Panagiotis (pronounced Panayiotis) Tzirakis (master student at the department of computer science in Heraklion during the 2013-2015 seasons) showed me how to perform parallel computing in R and he is greatly acknowledged and appreciated not only from me but from all the readers of this document. He also helped me with the vectorization of some contour plot functions.

Professor John Kent from the university of Leeds is acknowledged for clarifying one thing with the ovalness parameter in his distribution.

Phillip Paine (postdoc at the university of Nottingham) spotted that the function \code{\link{rfb}} is rather slow and he suggested me to change it. The function has changed now and this is also due to Joshua Davis (from Carleton College, Northfield, MN) who spotted that mistakes could occur, due a vector not being a matrix.

Professor Kurt Hornik from the Vienna university of economics and business is greatly acknowledged for his patience and contast help with this (and not only) R package.

Manos Papadakis is also acknowledged for his programming tips and for his assistance with the "htest" class object.

Dr Mojgan Golzy spotted a mistake in the function \code{\link{desag}} and Michail is very happy for that.

Lisette de Jonge-Hoekstra from the University of Groningen found a wrong sentence in the help file of function \code{\link{spml.reg}} which is now deleted.

Peter Harremoes from the Copenhagen Business College spotted a mistake in the confidence interval of the function  \code{\link{circ.summary}} which has now been corrected.

Dr Gregory Emvalomatis from the University of Crete helped me understand better the EM algorithm for mixture models and I fixed a bug in the function \code{\link{mixvmf.mle}}.

Kinley Russell, PhD student at the Johns Hopkins University School of Medicine, suggested that I include bootstrap ANOVA functions.

Sia Ahmadi found a mistake in the function \code{\link{conc.test}} which has now been corrected.

Rafail Vargiakakis found a bug in the functions \code{\link{circbeta.mle}}, \code{\link{circexp.mle}} and \code{\link{mmvm.mle}}.

If you want more information on many of these algorithms see Chapters 9 and 10 in the following document.
https://www.researchgate.net/publication/324363311_Multivariate_data_analysis_in_R
}

\keyword{ Directional data }
\keyword{ Spherical data }
\keyword{ circular data }
\keyword{ Regression }
\keyword{ Anova }
\keyword{ Discriminant analysis }
\keyword{ Graphs }
\keyword{ Simulation }


