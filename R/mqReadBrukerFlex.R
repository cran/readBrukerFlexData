## Copyright 2010-2011 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of readBrukerFlexData for R and related languages.
##
## readBrukerFlexData is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## readBrukerFlexData is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with readBrukerFlexData. If not, see <http://www.gnu.org/licenses/>

## function mqReadBrukerFlex 
##  imports fid files into MALDIquant MassSpectrum class 
##
##  WARNING: this is a recursive function!
##
## params:
##  path: path to root dir of fid files/or a single fid file e.g. data/
##
## returns:
##  a MALDIquant MassSpectrum 
##
mqReadBrukerFlex <- function(path, ...) {

    if (!file.exists(path)) {
        stop("Path ", sQuote(path), " doesn't exists!");
    }

    if (!require("MALDIquant")) {
        stop("Could not load package ", sQuote("MALDIquant"), ".");
    }

    if (!file.info(path)$isdir) {
        s <- readBrukerFlexFile(fidFile=path, ...);
        return(createMassSpectrum(mass=s$spectrum$mass,
                                  intensity=s$spectrum$intensity,
                                  metaData=s$metaData));
    } else {
        s <- readBrukerFlexDir(brukerFlexDir=path, ...);
        s <- lapply(s, function(x) {
                    return(createMassSpectrum(mass=x$spectrum$mass,
                                              intensity=x$spectrum$intensity,
                                              metaData=x$metaData)); });
        if (length(s) == 1) {
            return(s[[1]]);
        } else {
            return(s);
        }
    }
}

