\name{mqReadBrukerFlex}
\alias{mqReadBrukerFlex}

\title{
Reads mass spectrometry data into MALDIquant.
}
\description{
Reads all mass spectrometry data in Bruker Daltonics' XMASS format in
a specified path into \code{\link[MALDIquant]{MALDIquant-package}} 
\code{\link[MALDIquant]{MassSpectrum-class}} objects.
}
\usage{
mqReadBrukerFlex(path, \dots)
}
\arguments{
    \item{path}{
        Path to \emph{directory} or a single \emph{fid} file.
}
    \item{\dots}{
        arguments to be passed to
        \code{\link[readBrukerFlexData]{readBrukerFlexDir}} or
        \code{\link[readBrukerFlexData]{readBrukerFlexFile}} 
}
}
\details{
See \code{\link[readBrukerFlexData]{readBrukerFlexDir}} or
\code{\link[readBrukerFlexData]{readBrukerFlexFile}}.
}
\value{
    Returns a \code{\link[MALDIquant]{MassSpectrum-class}} object
    if \emph{path} is/contains only one fid file. \cr
    Returns a list of \code{\link[MALDIquant]{MassSpectrum-class}} objects
    if \emph{path} contains multiple fid files.
}
\references{                                                                                                                                                                                 
See website: \url{http://strimmerlab.org/software/maldiquant/}                                                                                                                               
}  
\author{
Sebastian Gibb
}
\seealso{
\code{\link[readBrukerFlexData]{readBrukerFlexDir}},
\code{\link[readBrukerFlexData]{readBrukerFlexFile}},
\code{\link[MALDIquant]{MALDIquant-package}},
\code{\link[MALDIquant]{MassSpectrum-class}}
}
\examples{
## load library
library("readBrukerFlexData");

## get examples directory
exampleDirectory <- system.file("Examples", package="readBrukerFlexData");

## read example spectra
spec <- mqReadBrukerFlex(file.path(exampleDirectory, "2010_05_19_Gibb_C8_A1"));

## plot spectra
par(mfrow=c(2, 1));
lapply(spec, function(x)plot(x));
par(mfrow=c(1, 1));
}
\keyword{IO}
