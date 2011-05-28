\name{readBrukerFlexDir}
\alias{readBrukerFlexDir}

\title{
Reads recursively mass spectrometry data in Bruker Daltonics' XMASS format.
}
\description{
Reads recursively all mass spectrometry data in Bruker Daltonics' XMASS format in
a specified directory.
}
\usage{
readBrukerFlexDir(brukerFlexDir, removeCalibrationScans=TRUE, 
    removeMetaData=FALSE, useHpc=TRUE, useSpectraNames=TRUE, verbose=FALSE)
}
\arguments{
    \item{brukerFlexDir}{
        Path to \emph{directory} which should be read recursively.
}
    \item{removeCalibrationScans}{
        logical. If TRUE all scans in directories called [Cc]alibration will be
        ignored.
}
    \item{removeMetaData}{
        logical. To calculate mass data a lot of meta data are needed. To save
        memory they could delete after calculation. 
}
    \item{useHpc}{
        logical. Should Bruker Daltonics' High Precision
        Calibration be used if available? (see also: \code{\link[readBrukerFlexData]{.hpc}})
}
    \item{useSpectraNames}{
        logical. If TRUE all list elements get an unique name from metaData
        otherwise file path is used.
        (If \sQuote{removeMetaData} is TRUE \sQuote{useSpectraNames} has no 
        effect.)
}
    \item{verbose}{
        logical. Print verbose messages?
}
}
\details{
see \code{\link[readBrukerFlexData]{readBrukerFlexFile}}
}
\value{
A list of spectra.

    \item{[[1]]$spectrum$mass}{
        A vector of calculated mass.
}
    \item{[[1]]$spectrum$intensity}{
        A vector of intensity values. 
}
    \item{[[1]]$metaData}{
        A list of metaData depending on read spectrum.
}
}
\author{
Sebastian Gibb
}
\seealso{
    \code{\link[readBrukerFlexData]{readBrukerFlexFile}},
    \code{\link[readBrukerFlexData]{.hpc}}
}
\examples{
## load library
library("readBrukerFlexData");

## get examples directory
exampleDirectory <- system.file("Examples", package="readBrukerFlexData");

## read example spectra
spec <- readBrukerFlexDir(file.path(exampleDirectory, "2010_05_19_Gibb_C8_A1"));

## plot spectra
plot(spec[[1]]$spectrum$mass, spec[[1]]$spectrum$intensity, type="n");

l <- length(spec);
legendStr <- character(l);
for (i in seq(along=spec)) {
    lines(spec[[i]]$spectrum$mass, spec[[i]]$spectrum$intensity, type="l",
    col=rainbow(l)[i]);
    legendStr[i] <- spec[[i]]$metaData$fullName;
}

## draw legend
legend(x="topright", legend=legendStr, col=rainbow(l), lwd=1);
}
\keyword{IO}
