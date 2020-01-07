\name{Directional-package}
\alias{Directional-package}
\docType{package}
\title{
This is an R package that provides methods for the statistical analysis of directional data, including massive (very large scale) directional data.
}
\description{
Circular-linear regression, spherical-spherical regression, spherical regression, discriminant analysis, ANOVA for circular and (hyper-)spherical data, tests for eaquality of conentration parameters, fitting distributions, random values generation, contour plots and many more functions are included in this package.

}
\details{
\tabular{ll}{
Package: \tab Directional\cr
Type: \tab Package\cr
Version: \tab 4.1\cr
Date: \tab 2020-01-07\cr
License: \tab GPL-2\cr
}
}
\author{
Michail Tsagris \email{mtsagris@uoc.gr}, Giorgos Athineou <gioathineou@gmail.com>,
Anamul Sajib <sajibstat@du.ac.bd>, Eli Amson <eli.amson1988@gmail.com> and Micah J. Waldstein <micah@waldste.in>.
}

\section{Maintainers }{
Michail Tsagris \email{mtsagris@uoc.gr}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley and Sons.
}

\note{
\bold{Acknowledgments:}

Professor Andy Wood and Dr Simon Preston from the university of Nottingham are highly appreciated for being my supervisors during my post-doc in directional data analysis.

Dr Georgios Pappas (former postDoc at the university of Nottingham) helped me construct the contour plots of the von Mises-Fisher and the Kent distribution.

Dr Christopher Fallaize and Dr Theo Kypraios from the university of Nottingham have provided a function for simulating from the Bingham distribution using rejection sampling.
So any questions regarding this function should be addressed to them.

Dr Kwang-Rae Kim (post-doc at the university of Nottingham) answered some of my questions.

Giorgos Borboudakis (PhD student at the university of Crete) pointed out to me a not so clear message in the algorithm of generating random values from the von Mises-Fisher distribution.

Panagiotis (pronounced Panayiotis) Tzirakis (master student at the department of computer science in Heraklion during the 2013-2015 seasons) showed me how to perform parallel computing in R and he is greatly acknowledged
and appreciated not only from me but from all the readers of this document. He also helped me with the vectorization of some contour plot functions.

Professor John Kent from the university of Leeds is acknowledged for clarifying one thing with the ovalness parameter in his distribution.

Phillip Paine (postdoc at the university of Nottingham) spotted that the function \code{\link{rfb}} is rather slow and he suggested me to change it. The function has changed now and this is also due to Joshua Davis
(from Carleton College, Northfield, MN) who spotted that mistakes could occur, due a vector not being a matrix.

Professor Kurt Hornik from the Vienna university of economics and business is greatly acknowledged for his patience and contast help with this (and not only) R package.

Manos Papadakis, undergraduate student in the department of computer science at university of Crete, is also acknowledged for his programming tips.

Dr Mojgan Golzy spotted a mistake in the \code{\link{ESAGdensity}} and Michail is very happy for that.

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


