\name{Directional-package}
\alias{Directional-package}
\docType{package}
\title{
This is an R package that provides methods for the statistical analysis of Directional data.
}
\description{
Circular-linear regression, spherical-spherical regression, discriminant analysis, ANOVA for circular and (hyper-)spherical data, tests for eaquality of conentration parameters, fitting distributions, random values generation, contour plots and many more are included in this package.

}
\details{
\tabular{ll}{
Package: \tab Directional\cr
Type: \tab Package\cr
Version: \tab 2.9\cr
Date: \tab 2017-10-17\cr
License: \tab GPL-2\cr
}
}
\author{
Michail Tsagris <mtsagris@yahoo.gr>, Giorgos Athineou <gioathineou@gmail.com> and Anamul Sajib <pmxahsa@nottingham.ac.uk>
}

\section{Maintainers }{
Michail Tsagris <mtsagris@yahoo.gr>
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley and Sons.
}

\note{
Acknowledgments:

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

Dr Eli Amson, postdoctoral fellow at Alexander von Humboldt Foundation found a bug in the \code{\link{spherconc.test}} and he is greatly acknowledged for that.

Dr Mojgan Golzy spotted a mistake in the \code{\link{ESAGdensity}} and Michail is very happy for that.
}

\keyword{ Directional data }
\keyword{ Spherical data }
\keyword{ circular data }
\keyword{ Regression }
\keyword{ Anova }
\keyword{ Discriminant analysis }
\keyword{ Graphs }
\keyword{ Simulation }


