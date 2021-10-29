\name{Maps of the world and the continents}
\alias{asia}
\alias{africa}
\alias{europe}
\alias{north.america}
\alias{oceania}
\alias{south.america}
\alias{worldmap}
\title{
 maps of the world and the continents
}
\description{
It produces maps of the world and the continents.
}
\usage{
asia(title = "Asia")
africa(title = "Africa")
europe(title = "Europe")
north.america(title = "North America")
oceania(title = "Oceania")
south.america(title = "South America")
worldmap(title = "World map")
}
\arguments{
\item{title}{
A character vector with the title of the map.
}
}
\details{
Maps of the world or the continents are produced. This are experimental functions and plot the countries with specific
colouring at the moment. More functionalities will be added in the future.
}
\value{
A map of the selected continent or the whole world.
}

\author{
Christos Adam.

R implementation and documentation: Christos Adam \email{pada4m4@gmail.com} and Michail Tsagris.
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{sphereplot}
}
}
\examples{
\dontrun{
oceania()
}
}
\keyword{ Median direction }
