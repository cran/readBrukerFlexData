## $Id: readBrukerFlexData.R 381 2011-02-15 15:58:49Z sgibb $
##
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

## function readBrukerFlexDir
##  reads all Bruker Daltonics' fid files in one directory
##
##  WARNING: this is a recursive function!
##
## params:
##  brukerFlexDir: path to root dir of fid files e.g. data/
##  removeCalibrationScans: default TRUE, don't read spectra from calibration
##                          scans
##  removeMetaData: see .readBrukerFlexFile for details (default: FALSE)
##  useHpc: TRUE/FALSE see .hpc/.readBrukerFlexFile for details 
##          [default: useHPC=TRUE]
##  useSpectraNames: TRUE/FALSE [default: useSpectraNames=TRUE]
##  verbose: TRUE/FALSE [default: verbose=FALSE]
##
##  files:
## ./hc/Pankreas_HB_L_061019_A7/0_a13/1/1SLin/fid
## ./hc/Pankreas_HB_L_061019_A7/0_a14/1/1SLin/fid
## ./hc/Pankreas_HB_L_061019_A7/0_b13/1/1SLin/fid
## ./hc/Pankreas_HB_L_061019_A7/0_b14/1/1SLin/fid
## ./hc/Pankreas_HB_L_061019_B4/0_c8/1/1SLin/fid
##
##  brukerFlexDir="./hc"
##
## returns:
##  a list with metadata and spectra
##
readBrukerFlexDir <- function(brukerFlexDir, removeCalibrationScans=TRUE,
        removeMetaData=FALSE, useHpc=TRUE, useSpectraNames=TRUE,
        verbose=FALSE) {
    if (verbose)
        message("Look for spectra in ", sQuote(brukerFlexDir), " ...");

    if ((!file.exists(brukerFlexDir)) || (!file.info(brukerFlexDir)$isdir)) {
        warning("Directory ", sQuote(brukerFlexDir), " doesn't exists or is no
                directory!");
        return(NA);
    }

    ## look for fid files (alphabetical sort)
    files <- list.files(path=brukerFlexDir, pattern="^fid$", recursive=TRUE);

    ## remove calibrations scans?
    if (removeCalibrationScans) {
        calibrationScans <- grep(pattern="[Cc]alibration", x=files, value=TRUE);
        if (length(calibrationScans) > 0) {
        files <- setdiff(files, calibrationScans);
        }
    }

    ## resort files in correct order, using gtools package
    ## A1, A2, ..., A10, ..., B1, ...
    #files <- .mixedsort(files);

    ## sometimes .mixedsort produces NAs
    #files <- files[!is.na(files)];

    ## generate "path/files"
    files <- sapply(files, function(x) {
            x <- paste(brukerFlexDir, x, sep="/");
            return(x);
        });

    ## read fid files
    brukerFlexData <- lapply(X=files, FUN=function(f) {
            return(readBrukerFlexFile(fidFile=f, removeMetaData=removeMetaData,
                    useHpc=useHpc, verbose=verbose)); });

    if (!removeMetaData & useSpectraNames) {
        ## rewrite names if metadata exists
        if (verbose)
            message("look for spectra names ...");

        names(brukerFlexData) <- sapply(X=brukerFlexData, FUN=function(x) {
                if (!is.null(x$metaData$sampleName)) {
                    if (!is.na(x$metaData$sampleName)) {
                        return(paste("s", x$metaData$fullName, ".",
                                x$metaData$targetIdString, sep=""));
                    }
                }
                return(NA);
            }, USE.NAMES=FALSE);
    }

    return(brukerFlexData);
}

## function readBrukerFlexFile
##  read a single spectrum from Bruker Daltonics' fid and acqu file
##
## params:
##  fidFile: path to fid file e.g. Pankreas_HB_L_061019_A10/0_a19/1/1SLin/fid
##  removeMetaData: if TRUE => don't return metadata to save memory 
##                  [default: removeMetaData=TRUE]
##  useHpc: TRUE/FALSE use HPC if avaiable? [default: useHpc=TRUE]
##  verbose: TRUE/FALSE [default: verbose=FALSE]
##
## returns:
##  a list with intensity, mass and metadata (if removeMetaData == FALSE)
##
##  spectrum$intensity, spectrum$tof, spectrum$mass, metadata
##
readBrukerFlexFile <- function(fidFile, removeMetaData=FALSE, useHpc=TRUE,
    verbose=FALSE) {
    ## try to get absolute file path
    fidFile <- normalizePath(fidFile);

    if (verbose)
        message("Reading spectrum from ", sQuote(fidFile), " ...");
  
    if (!file.exists(fidFile)) {
        warning("File ", sQuote(fidFile), " doesn't exists!");
        return(NA);
    }

    if (file.info(fidFile)$isdir) {
        warning("Not a fid file! ", sQuote(fidFile), " is a directory.");
        return(NA);
    }

    ## read metadata
    metaData <- .readAcquFile(fidFile=fidFile, verbose=verbose);

    ## read peak intensities 
    intensity <- .readFidFile(fidFile, metaData$number);

    ## calculate tof of metadata
    tof <- metaData$timeDelay + ((0:(metaData$number-1)) * metaData$timeDelta);

    ## remove times which have no intensity (intensity == 0), e.g. generated by
    ## REFLECTOR mode
    ## for details see "Release Notes for CompassXport 3.0.3" 
    ## cap. 6 "Filtering of Zero Intensities"
    ## "Bruker Daltonics’ Aqcuisition Software will compress Analysis raw 
    ## data. To save on operation time and to keep export file sizes small,
    ## CompassXport 3.0.3 will filter out zero (0.0) intensities
    ## when exporting to mzXML or mzData ..."
    notNull <- intensity>0;
    intensity <- intensity[notNull];
    tof <- tof[notNull];

    ## calculate mass of TOFs
    mass <- .tof2mass(tof, 
                      metaData$calibrationConstants[[1]],
                      metaData$calibrationConstants[[2]],
                      metaData$calibrationConstants[[3]]);
  
    ## TODO: fix equations in .hpc

    ## was HPC involved?
    ## metaData$hpcUse seems to be always true
    isHPCused <- (useHpc &
                  metaData$hpcUse &
                  metaData$hpcLimits["maxMass"] > 0 &
                  metaData$hpcLimits["minMass"] > 0 & 
                  metaData$hpcOrder > 0);

    if (isHPCused) {
        ## TODO: fix equations in .hpc and remove the following warning 
        warning("The spectrum file ", sQuote(fidFile), " uses HPC.\n",
                "HPC isn't fully supported by readBrukerFlexFile",
                "Please see ", dQuote("?.hpc"), " for details.\n",
                "Original mass are ", sQuote("metaData$backup$mass"), ".");
        metaData$backup$mass <- mass;

        mass <- .hpc(mass=mass,
                     minMass=metaData$hpcLimits["minMass"],
                     maxMass=metaData$hpcLimits["maxMass"],
                     hpcCoefficients=metaData$hpcCoefficients);
    }
  
    if (!removeMetaData) {
        return(list(spectrum=list(tof=tof, mass=mass, intensity=intensity), metaData=metaData));
    } else {
        return(list(spectrum=list(tof=tof, mass=mass, intensity=intensity)));
    }
}

## function .readAcquFile
##  reads acqu file
##
##  we have to import the following data to calculating mass:
##  $TD: total number of measured time periods
##      => metaData$number
##  $DELAY: first measured intensity after $DELAY ns
##      => metaData$timeDelay
##  $DW: ns between measured time periods
##      => metaData$timeDelta
##  $ML1: calibration constant
##      => metaData$calibrationConstants[1]
##  $ML2: calibration constant
##      => metaData$calibrationConstants[2]
##  $ML3: calibration constant
##      => metaData$calibrationConstants[3]
##
##  if we want to use High Precision Calibration (HPC), we need:
##  $HPClBHi: upper mass treshold
##      => metaData$hpcLimits["maxMass"]
##  $HPClBLo: lower mass threshold
##      => metaData$hpcLimits["minMass"]
##  $HPClOrd: polynom order
##      => metaData$hpcOrder
##  $HPClUse: maybe using of HPC? (seems always be "yes" in our test data)
##      => metaData$hpcUse
##  $HPCStr: polynom coefficients in a string
##      => metaData$hpcCoefficients
##
##  we try to import [optinal]:
##  DATATYPE
##  SPECTROMETER/DATASYSTEM
##      => metaData$dataSystem
##  .SPECTROMETER TYPE
##      => metaData$spectrometerType
##  .INLET
##      => metaData$inlet
##  .IONIZATION MODE
##      => metaData$ionizationMode
##  $DATE
##      => metaData$date
##  $ACQMETH
##      => metaData$acquisitionMethod
##  $AQ_DATE
##      => metaData$acquisitionDate
##  $AQ_mod
##      => metaData$acquisitionMode
##  $AQOP_mod
##      => metaData$acquisitionOperatorMode
##      => metaData$tofMode (replaces file path based method)
##  $ATTEN
##      => metaData$laserAttenuation
##  $COM[1:4]
##      => metaData$comments
##  $DEFLON
##      => metaData$deflection
##  $DIGTYP
##      => metaData$digitizerType
##  $DPCAL1
##      => metaData$deflectionPulserCal1
##  $DPMASS
##      => metaData$deflectionPulserMass
##  $FCVer
##      => metaData$flexControlVersion
##  $ID_raw
##      => metaData$id
##  $INSTRUM
##      => metaData$instrument
##  $InstrID
##      => metaData$instrumentId
##  $InstrTyp
##      => metaData$instrumentType
##  $Masserr
##      => metaData$massError
##  $NoSHOTS: number of laser shots
##      => metaData$laserShots
##  $SPType
##      => metaData$spectrumType
##  $PATCHNO: sample postion on target
##      => metaData$patch
##  $PATH: original file path (on bruker flex series controller pc)
##      => metaData$path
##  $REPHZ
##      => metaData$laserRepetition
##  $SPOTNO: same as $PATCHNO (in older files often empty, thats why we use 
##      $PATHNO instead)
##      => metaData$spot
##  $TgIDS: target ids
##      => metaData$targetIdString
##  $TgCount: number of measurement with this target
##      => metaData$targetCount
##  $TgSer: target serial number
##      => metaData$targetSerialNumber
##  $TgTyp: target type number
##      => metaData$targetTypeNumber
##
##  import from file path:
##  full current path to fid file:
##      => metaData$fidFile
##  sample name
##      => metaData$sampleName
##
## params:
##  file: path to corresponding fid file
##           e.g. Pankreas_HB_L_061019_A10/0_a19/1/1SLin/fid
##  verbose: TRUE/FALSE [default: verbose=FALSE]
##
## returns:
##  a list with metadata
##
##  $number
##  $timeDelay
##  $timeDelta
##  $calibrationConstants
##      c1            c2            c3
##
##  $hpc
##  $hpcLimits
##      minMass   maxMass
##  $hpcOrder
##  $hpcUse
##  $hpcCoefficients
##  $hpcCalibrationConstant0
##  $hpcCalibrationConstant2
##
##  $dataType
##  $dataSystem
##  $spectrometerType
##  $inlet
##  $ionizationMode
##  $date
##  $acquisitionMethod
##  $acquisitionDate
##  $acquisitionMode
##  $acquisitionOperatorMode
##  $comments
##  $deflection
##  $digitizerType
##  $deflectionPulserCal1
##  $deflectionPulserMass
##  $file
##  $flexControlVersion
##  $fullName
##  $id
##  $instrument
##  $instrumentId
##  $instrumentType
##  $laserAttention
##  $laserShots
##  $laserRepetition
##  $massError
##  $sampleName
##  $shots
##  $spectrumType
##  $spot
##  $patch
##  $path
##  $targetCount
##  $targetIdString
##  $targetSerialNumber
##  $targetTypeNumber
##  $tofMode
##
.readAcquFile <- function(fidFile, verbose=FALSE) {
    acquFile <- sub(pattern="/fid$", x=fidFile, replacement="/acqu");

    if (verbose)
        message("Reading metadata from ", sQuote(acquFile), " ...");

    if (!file.exists(acquFile)) {
        warning("File ", sQuote(acquFile), " doesn't exists!");
        return(NA);
    }

    con <- file(acquFile, "rt");
    acquLines <- readLines(con, n=-1);
    close(con);

    ## collect data
    metaData <- list();

    ## obligate
    metaData$number <- as.double(.grepAcquValue("##\\$TD=", acquLines));
    metaData$timeDelay <- as.double(.grepAcquValue("##\\$DELAY=", acquLines));
    metaData$timeDelta <- as.double(.grepAcquValue("##\\$DW=", acquLines));
    metaData$calibrationConstants <- c(
        c1=as.double(.grepAcquValue("##\\$ML1=", acquLines)), 
        c2=as.double(.grepAcquValue("##\\$ML2=", acquLines)),
        c3=as.double(.grepAcquValue("##\\$ML3=", acquLines)));

    ## obligate HPC
    metaData$hpcLimits <- c(
        minMass=as.double(.grepAcquValue("##\\$HPClBLo=", acquLines)),
        maxMass=as.double(.grepAcquValue("##\\$HPClBHi=", acquLines)));
    metaData$hpcOrder <- as.double(.grepAcquValue("##\\$HPClOrd=", 
        acquLines));
    metaData$hpcUse <- 
        as.logical(.grepAcquValue("##\\$HPClUse=", acquLines) == "yes");
  
    ## was HPC involved?
    ## metaData$hpcUse seems to be always true
    isHPCused <- (metaData$hpcUse && 
                  metaData$hpcLimits["maxMass"] > 0 &&
                  metaData$hpcLimits["minMass"] && 
                  metaData$hpcOrder > 0);

    if (isHPCused) {
        hpcStr <- .grepAcquValue("##\\$HPCStr=", acquLines);
        hpcConstants <- .extractHPCConstants(hpcStr);
        metaData$hpcCoefficients <- hpcConstants$coefficients;
        metaData$hpcCalibrationConstant0 <- hpcConstants$calibrationConstant0;
        metaData$hpcCalibrationConstant2 <- hpcConstants$calibrationConstant2;
    }

    ## optional
    metaData$dataType <- .grepAcquValue("##DATATYPE=", acquLines);
    metaData$dataSystem <- .grepAcquValue("##SPECTROMETER/DATASYSTEM=", 
        acquLines);
    metaData$spectrometerType <- .grepAcquValue("##.SPECTROMETER TYPE=",
        acquLines);
    metaData$inlet <- .grepAcquValue("##.INLET=", acquLines);
    metaData$ionizationMode <- .grepAcquValue("##.IONIZATION MODE=", acquLines);
    metaData$date <- .grepAcquValue("##\\$DATE=", acquLines);


    metaData$acquisitionMethod <- .grepAcquValue("##\\$ACQMETH=", acquLines);
    metaData$acquisitionDate <- .grepAcquValue("##\\$AQ_DATE=", acquLines);
    metaData$acquisitionMode <- .grepAcquValue("##\\$AQ_mod=", acquLines);

    aqop <- as.double(.grepAcquValue("##\\$AQOP_m=", acquLines));
    switch(aqop,
        "0" = {
            metaData$tofMode <- "LINEAR";
        },
        "1" = {
            metaData$tofMode <- "REFLECTOR";
        },
        {
            metaData$tofMode <- aqop;
        });

    metaData$acquisitionOperatorMode <- metaData$tofMode;

    metaData$laserAttenuation <- 
        as.double(.grepAcquValue("##\\$ATTEN=", acquLines));

    metaData$comments <- c(.grepAcquValue("##\\$CMT1=", acquLines),
                           .grepAcquValue("##\\$CMT2=", acquLines),
                           .grepAcquValue("##\\$CMT3=", acquLines),
                           .grepAcquValue("##\\$CMT4=", acquLines));

    metaData$deflection <- 
        as.logical(.grepAcquValue("##\\$DEFLON=", acquLines) == "yes");
    metaData$digitizerType <-
        as.double(.grepAcquValue("##\\$DIGTYP=", acquLines));
    metaData$deflectionPulserCal1 <- 
        as.double(.grepAcquValue("##\\$DPCAL1=", acquLines));
    metaData$deflectionPulserMass <- 
        as.double(.grepAcquValue("##\\$DPMASS=", acquLines));
    metaData$flexControlVersion <- .grepAcquValue("##\\$FCVer=", acquLines);
    metaData$id <- .grepAcquValue("##\\$ID_raw=", acquLines);

    metaData$instrument <- .grepAcquValue("##\\$INSTRUM=", acquLines);
    metaData$instrumentId <- .grepAcquValue("##\\$InstrID=", acquLines);
    metaData$instrumentType <- .grepAcquValue("##\\$InstTyp=", acquLines);

    metaData$massError <- 
        as.double(.grepAcquValue("##\\$Masserr=", acquLines));

    metaData$laserShots <- 
        as.double(.grepAcquValue("##\\$NoSHOTS=", acquLines));

    if (metaData$laserShots == 0) {
        warning("File ", sQuote(fidFile), " seems to be empty because
            no laser shots applied to this sample.")
    }

    metaData$patch <- .grepAcquValue("##\\$PATCHNO=", acquLines);
    metaData$path <- .grepAcquValue("##\\$PATH=", acquLines);
    metaData$laserRepetition <-
        as.double(.grepAcquValue("##\\$REPHZ=", acquLines));
    metaData$spot <- .grepAcquValue("##\\$SPOTNO=", acquLines);
    
    sptype <- as.double(.grepAcquValue("##\\$SPType=", acquLines));
    switch(sptype,
        "0" = {
            metaData$spectrumType <- "TOF";
        },
        {
            metaData$spectrumType <- sptype;
        });

    metaData$targetCount <- 
        as.double(.grepAcquValue("##\\$TgCount", acquLines));
    metaData$targetIdString <- .grepAcquValue("##\\$TgIDS", acquLines);
    metaData$targetSerialNumber <- .grepAcquValue("##\\$TgSer", acquLines);
    metaData$targetTypeNumber <- .grepAcquValue("##\\$TgTyp", acquLines);

    ## from file path 

    #if (grepl(pattern="/1S?Ref/", x=fidFile, fixed=FALSE)) {
    #    metaData$tofMode <- "REFLECTOR";
    #} else {
    #    metaData$tofMode <- "LINEAR";
    #}

    metaData$file <- fidFile;

    metaData$sampleName <- .sampleName(fidFile);
    metaData$fullName <- paste(metaData$sampleName, metaData$patch, sep=".");

    return(metaData);
}

## function .readFidFile
##  reads binary fid file
##
##  fid files contain intensities for all time periods
##
## params:
##  fidFile: path to fid file e.g. Pankreas_HB_L_061019_A10/0_a19/1/1SLin/fid
##  nIntensities: number of data entries (total count; get from acqu file)
##
## returns:
##  a vector of intensity values 
##
.readFidFile <- function(fidFile, nIntensities) {
    if (!file.exists(fidFile)) {
        warning("File ", sQuote(fidFile), " doesn't exists!");
        return(NA);
    }

    con <- file(fidFile, "rb");
    intensity <- readBin(con, integer(), n=nIntensities, size=4, endian="little");
    close(con);

    return(intensity);
}

###################
## helper functions
###################

## function .grepAcquValue
##  helper function to extract values from acqu file
##
## params:
##  patternStr: pattern to look for
##  srcStr: where to look for patternStr
##
## returns:
##  character vector of the value given in patternStr
##
.grepAcquValue <- function(patternStr, srcStr) {
    tmpLine <- grep(pattern=patternStr, x=srcStr, value=TRUE);

    ## format e.g.
    ##  DATATYPE=  CONTINUOUS MASS SPECTRUM
    ##  .IONIZATION MODE=  LD+
    ##  $INSTRUM= <AUTOFLEX>

    ## remove front pattern
    tmpLine <- sub("^.*= *<?", replacement="", tmpLine);
    ## remove back pattern
    tmpLine <- sub(">? *$", replacement="", tmpLine);

    return(tmpLine);
}

## function .sampleName
##  guess name of current spot from filename
##
##  WARNING: if the 4th upper dir hasn't an unique name you will get
##  equal names for your spot list
##
##  e.g. the following will create a list with 4 elements but
##  only 2 unique spot names (2-100kDa:0_A1 and 2-100kDa:0_B1)
## ./Run1/2-100kDa/0_A1/1/1SLin/fid
## ./Run1/2-100kDa/0_B1/1/1SLin/fid
## ./Run2/2-100kDa/0_A1/1/1SLin/fid
## ./Run2/2-100kDa/0_B1/1/1SLin/fid
##
## params:
##  fidFile: path to fid file e.g. Pankreas_HB_L_061019_A10/0_a19/1/1SLin/fid
##
## returns:
##  sampleName: e.g Pankreas_HB_L_061019_A10
##
.sampleName <- function(fidFile) {
  # remove double slashes created by pasting path="example/path/" and "/file"
  fidFile <- gsub(pattern="//+", replacement="/", x=fidFile);

  # create array of directories (each element == one directory)
  dirs <- strsplit(x=fidFile, split="/", fixed=TRUE)[[1]];

  numDirs <- length(dirs);

  sampleName <- NA;

  if (numDirs > 4) {
    sampleName <- dirs[numDirs-4];

    # -, : or something like that causes errors in names()
    sampleName <- gsub(pattern="[[:punct:]]|[[:space:]]", replacement="_",
            x=sampleName);

  } 

  return(sampleName);
}

## function .extractHPCConstants
##  helper function to extract coefficients and constants values from 
##  metaData$hpcStr
##
## params:
##  hpcStr: metaData$hpcStr, which store coefficents
##      e.g. hpcStr <- " V1.0CHPCData  Order 10 vCoeff V1.0VectorDouble 11 
##                      -0.48579224953906053 0.0009361303203700988 
##                      -6.92711401708155e-008 -1.0992953299897006e-009
##                      1.1718229914003113e-012 -5.392578762547374e-016
##                      9.0176664604755316e-020 1.9704001597871883e-023 
##                      -1.1794161284667635e-026 2.0351573912658823e-030
##                      -1.2617853301428769e-034  
##                      c2 -0.046701600316874939 
##                      c0 237.64781433281422 
##                      minMass 736.50266799999997
##                      maxMass 3698.6377320000001 bUse 1 endCHPCData "
##
## returns:
##  list of doubles;
##  hpcConstants$coefficients: double vector of coefficents
##  hpcConstants$calibrationContant0: c0
##  hpcConstants$calibrationContant2: c2
##
.extractHPCConstants <- function(hpcStr) {
    tmpLine <- strsplit(x=hpcStr, split=" ")[[1]];
    ## remove emtpy elements
    tmpLine <- tmpLine[tmpLine!=""];

    hpcConstants <- list();

    ## extract only coefficients
    hpcConstants$coefficients <- as.double(
        tmpLine[ (which(tmpLine == "V1.0VectorDouble")+2) :
                 (which(tmpLine == "c2")-1) ]);

    hpcConstants$calibrationConstant2 <- 
        as.double(tmpLine[which(tmpLine == "c2")+1]);
    hpcConstants$calibrationConstant0 <- 
        as.double(tmpLine[which(tmpLine == "c0")+1]);

    return(hpcConstants);
}

## function .tof2mass
##  calculate mass from time of flight values
##
##  based on the following article:
##  Mark K Titulaer, Ivar Siccama, Lennard J Dekker, Angelique LCT van Rijswijk,
##  Ron MA Heeren, Peter A Sillevis Smitt, and Theo M Luider
##  "A database application for pre-processing, storage and comparison of mass spectra
##   derived from patients and controls"
##  BMC Bioinformatics. 2006; 7: 403
##  http://www.ncbi.nlm.nih.gov/pubmed/16953879
##
##  params are imported from metadata (acqu-file)
##
## params:
##  tof: vector with times-of-flight
##  c(1:3): metaData$calibrationConstants[1:3]
##
## returns:
##  double vector of mass
##
.tof2mass <- function(tof, c1, c2, c3) {
    ## 0 = A * (sqrt(m/z))^2 + B * sqrt(m/z) + C(times)
    A <- c3;
    B <- sqrt(1e12/c1);
    C <- c2 - tof;

    return( ( (-B + sqrt(B^2 - (4 * A * C)))/(2 * A) )^2  );
}

## function .hpc
##  support of Bruker Daltonics' High Precision Calibration (HPC)
##  Please note that .hpc is not correct! You have been warned.
##
##  maybe HPC based on the following article:
##  Johan Gobom, Martin Mueller, Volker Egelhofer, Dorothea Theiss,
##  Hans Lehrach, and Eckhard Nordhoff
##  "A Calibration Method That Simplifies and Improves Accurate Determination
##   of Peptide Molecular mass by MALDI-TOF MS"
##  Anal Chem. 2002 Aug 1; 74(15): 3915-23
##  http://www.ncbi.nlm.nih.gov/pubmed/12175185
##
##  but using this article for calibration out of acqu files is definitfly
##  the wrong way
##
##  get formula by trying a lot of stupid thinks by hand
##  (trial and error)
##
##  params are imported from metadata (acqu-file)
##
## TODO: 
##  - internal calibration (or something like that)
##  - maybe need the following: 
##      ##$Hpcgc0= 237.647814332814                                                                                                                       
##      ##$Hpcgc2= -0.0467016003168749 
##
## params:
##  mass: vector with alreday calculated mass (used .tof2mass)
##  minMass: metaData$hpc$minMass
##  maxMass: metaData$hpc$maxMass
##  hpcCoeefficents: metaData$hpc$coefficients
##
## returns:
##  double vector of calibrated mass
##
.hpc <- function(mass, minMass, maxMass, hpcCoefficients) {
    ## only defined for mass between minMass and maxMass
    ## QUESTION: Do we exclude values after or before calibration?
    ## ANSWER: I don't know how compassXport do it.
    ##         before: reduce calculation time
    ##         after:  produce better results
    ##         I also don't know whether to use >= or > (<=, <).
    m <- mass[mass>=minMass & mass<=maxMass];

    ## correction = c[0] + c[1]*cal_mass^1 + c[2]*cal_mass^2 + ... + c[n]*cal_mass^n
    ## mass = cal_mass - correction
    l <- length(hpcCoefficients)-1;
    m <- sapply(m, function(x) {
            return(x-sum( hpcCoefficients*x^(0:l) ));
    } );

    mass[mass>=minMass & mass<=maxMass] <- m;

    return(mass);
}

## EOF
