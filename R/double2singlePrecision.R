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

## This method simulate the conversion of floating point numbers from double
## precision (64bit, R: double(), C: double) to single precision (32bit, R:
## none, C: float)
## It follows IEEE 754 standard. http://754r.ucbtest.org/standards/754.pdf
##
## The same could be done in C by using casts:
## double precision32(double value) {
##   float x = value;
##   return (double)x;
## }
##

## function .double2singlePrecision
##  wrapper function for .changePrecision(x, size=4)
##  see also: .changePrecision
##
## params:
##  x: a vector of double values which should converted
##
## returns:
##  a vector of single precision values (represented by R double type)
##
.double2singlePrecision <- function(x) {
  stopifnot(is.double(x));
  return(.changePrecision(x, size=4));
}

## function .changePrecision
##  converts double values to double values in a given precision 
##  (only correctly working for cut a higher precision to a lower one; e.g.
##  IEEE 754 double precision to IEEE 754 single precision)
##
## params:
##  x: a vector of double values which should converted
##  size: an integer; how many bytes using for target precision
##        e.g. IEEE 754 single precision: size=4
##             IEEE 754 double precision: size=8
##
## returns:
##  a vector of modified values
##
.changePrecision <- function(x, size) {
  ## create a raw object to avoid direct file access
  virtualCon <- raw();
  ## write binary data to raw object and change (mostly cut) precision to size
  ## size==4 # 32bit, single precision
  ## size==8 # 64bit, double precision
  virtualCon <- writeBin(object=x, con=virtualCon, size=size);
  ## re-read data
  x <- readBin(con=virtualCon, what=double(), size=size, n=length(x));
  return(x);
}

