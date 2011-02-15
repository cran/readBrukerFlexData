\name{.double2singlePrecision}
\alias{.double2singlePrecision}

\title{
    Converts double to single precision.
}
\description{
    This method simulate the conversion of floating point numbers from double
    precision (64bit, R: double(), C: double) to single precision (32bit, R: 
    none, C: float). \cr
    It is only needed for comparison between Bruker Daltonics' HPC and 
    \code{\link{.hpc}}. \cr
    It follows IEEE 754 standard: 
    \url{http://754r.ucbtest.org/standards/754.pdf}
}
\usage{
.double2singlePrecision(x)
}
\arguments{
    \item{x}{
        double. Numeric value which should reduce from double precision to 
        single precision.
}
}
\details{
The same could be done in C by using casts: \cr
\preformatted{
double precision32(double value) {
    float x = value;
    return (double)x;
}
}
}
\value{
double (in single precision).
}
\references{
IEEE 754 standard \cr
\url{http://754r.ucbtest.org/standards/754.pdf}
}
\author{
Sebastian Gibb
}
\seealso{
    \code{\link{.hpc}}
}
\examples{
library(readBrukerFlexData);

# more details
options(digits=22);

# a test number
num <- 1/3;

num;
readBrukerFlexData:::.double2singlePrecision(num);
}
\keyword{misc}

