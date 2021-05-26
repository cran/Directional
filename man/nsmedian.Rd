\name{Normalised spatial median for directional data}
\alias{nsmedian}
\title{
Normalised spatial median for directional data
}

\description{
Normalised spatial median for directional data.
}

\usage{
nsmedian(x, tol = 1e-07)
}

\arguments{
\item{x}{
A matrix with Euclidean data, continuous variables.
}
\item{tol}{
A tolerance level to terminate the process.
}
}

\details{
The spatial median, using a fixed point iterative algorithm, for Euclidean data is calculated. It is a robust location estimate. Then it is normalised to become a unit vector. Generally speaking this might be a better alternative than then \code{\link{mediandir}}.
}

\value{
A vector with the spatial median.
}

\references{
Ducharme G. R. and Milasevic P. (1987). Spatial median and directional data. Biometrika, 74(1), 212-215.

Jyrki Mottonen, Klaus Nordhausen and Hannu Oja (2010). Asymptotic theory of the spatial median.
In Nonparametrics and Robustness in Modern Statistical Inference and Time Series Analysis:
A Festschrift in honor of Professor Jana Jureckova.

T. Karkkaminen and S. Ayramo (2005). On computation of spatial median for robust data mining.
Evolutionary and Deterministic Methods for Design, Optimization and Control with Applications to
Industrial and Societal Problems, EUROGEN 2005, R. Schilling, W.Haase, J. Periaux, H. Baier, G. Bugeda (Eds)
FLM, Munich. http://users.jyu.fi/~samiayr/pdf/ayramo_eurogen05.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mediandir}
}
}

\examples{
m <- rnorm(3)
m <- m / sqrt( sum(m^2) )
x <- rvmf(100, m, 10)
nsmedian(x)
mediandir(x)
}

