\name{Conversion of cosines to azimuth and plunge}
\alias{cosap}
\title{
Conversion of cosines to azimuth and plunge
}

\description{
Conversion of cosines to azimuth and plunge.
}

\usage{
cosap(x,y,z)
}

\arguments{
\item{x}{
x component of cosine.
}
\item{y}{
y component of cosine.
}
\item{z}{
z component of cosine.
}
}

\details{
Orientation: x>0 is 'eastward', y>0 is 'southward', and z>0 is 'downward'.
}

\value{
A list including:
\item{A}{
The azimuth
}

\item{P}{
The plunge}
}

\references{
Amson E, Arnold P, Van Heteren AH, Cannoville A, Nyakatura JA. Trabecular architecture in the forelimb epiphyses of extant xenarthrans (Mammalia). Frontiers in Zoology.
}

\author{
Eli Amson

R implementation and documentation: Eli Amson <eli.amson1988@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{euclid}, \link{euclid.inv}, \link{eul2rot}
}
}

\examples{
cosap(-0.505, 0.510, -0.696)
}

\keyword{ directional data }
