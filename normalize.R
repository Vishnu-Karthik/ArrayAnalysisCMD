    library(methods)
    library(argparser)
    options_ <- arg_parser("")
    options_ <- add_argument(options_, "data_directory", help="Path to directory with .cel files")
    options_ <- add_argument(options_, "output_directory", help="Path to the output directory")
    options_ <- add_argument(options_, "description_file", help="Path to description file")
    options_ <- add_argument(options_, "--normMeth", help="Normlization method (available options: 'RMA','GCRMA','PLIER','MAS5'")
    options_ <- add_argument(options_, "--customCDF", help="Custom CDF (Boolean)", flag=TRUE)
    options_ <- add_argument(options_, "--reOrder", help="Reorder plots (Boolean)", flag=TRUE)
    options_ <- add_argument(options_, "--CDFtype", help="Annotation output type. Default is ENSG.")
    options_ <- add_argument(options_, "--species", help="Species name. If left blank, attempts to deduce species")
    options_ <- add_argument(options_, "--perGroup", help="Normalize per group or per dataset. Default is dataset", flag=TRUE)
    options_ <- add_argument(options_, "--PMAcalls", help="Create a table of PMAcalls (PM-MM only)",flag=TRUE) 
    argv <- parse_args(options_)

    DATA.DIR <-argv$data_directory
    WORK.DIR <-argv$output_directory
    description_file <- argv$description_file
    SCRIPT.DIR <- "~/projects/sols/masters/4th_sem/scripts/"

    normMeth <- argv$normMeth
    if (is.na(normMeth)) {normMeth <- "RMA"}

    normOption1 <- argv$perGroup
    if (is.na(normOption1)) {normOption1<- "dataset"}

    customCDF <- argv$customCDF

    CDFtype <- argv$CDFtype
    if (is.na(CDFtype)) {CDFtype<- "ENSG"}

    species <- argv$species
    if (is.na(species)) {species <- NULL}

    reOrder <- argv$reOrder

    PMAcalls <- argv$PMAcalls

 version_nb <- "1.0.0"
 cat("Script run using R version ",R.Version()$major,".",R.Version()$minor,
   " and affyAnalysisQC version_",version_nb,"\n",sep="")

 if(length(grep("w32",R.Version()$os,fixed=TRUE))>0) memory.size(4095)

   require("affy", quietly = TRUE)
   require("affycomp", quietly = TRUE)
   require("affyPLM", quietly = TRUE)
   require("bioDist", quietly = TRUE)
   require("simpleaffy", quietly = TRUE)
   require("affyQCReport", quietly = TRUE)
   require("plier", quietly = TRUE)
   if(exists("samplePrep")) require("yaqcaffy", quietly = TRUE)
   require("gdata", quietly = TRUE) #trim function
   print("Libraries have been loaded")

   reload <- function() {
     source(paste(SCRIPT.DIR,"functions_processingQC.R",sep=""))
     source(paste(SCRIPT.DIR,"functions_imagesQC.R",sep=""))
     print ("Functions have been loaded")
   }
   reload();

   setwd(DATA.DIR)
   rawData <- ReadAffy()
   print("Raw data have been loaded in R")

   setwd(WORK.DIR)

   # Make sure that the CDF environment works
   rawData <- addStandardCDFenv(rawData)   # if already works, won't be changed

   # Verify the array type (PMMM or PMonly)
   aType <- getArrayType(rawData)

   # When refName does not exist, use the empty string
   if(!exists("refName")) refName <- ""

 WIDTH <- 1000
 HEIGHT <- 1414
 POINTSIZE <- 24
 if(!exists("maxArray")) maxArray <- 41

 if(description_file !=""){
   # Information is available: groups will be created
   # 1- read the description file and trim spaces
   # 2- define the array names and classes (experimentFactor)
   if(!exists("DESC.DIR")) DESC.DIR <- ""

   descfile <- paste(DESC.DIR, description_file, sep="")
   extension<-strsplit(descfile,"\\.")
   extension<-paste(".",extension[[1]][length(extension[[1]])],sep="")
   description = NULL;
   switch(extension,
          ".txt" = description<-trim(read.delim(descfile, fill = FALSE, as.is=TRUE)),
          ".csv" = description<-trim(read.csv(descfile, fill = FALSE, as.is=TRUE)),
          ".xls" = {library(gdata); description<-trim(read.xls(descfile, as.is=TRUE))},
          ".xlsx" = {library(gdata); description<-trim(read.xls(descfile, as.is=TRUE))}
	 )
   if(is.null(description)) stop(paste("extension",extension,"not recognised"))

  # description <- trim(read.table(paste(DESC.DIR, description_file , sep=""),
  # 	  header = TRUE, as.is = TRUE, sep="\t"))

   if(length(grep(".CEL",toupper(colnames(description)[1]),
     ignore.case = TRUE))>0) {
     stop(paste("The description file may not contain a header, as the first",
     	 "column header seems to be a CEL file name"))
   }
   file_order <- match(description[,1],sampleNames(rawData))
   if(sum(is.na(file_order)) > 0) stop("file names in data directory and file names in description file do not match")
   if(length(unique(file_order)) < length(file_order)) stop("file names in description file are not unique")
   rawData <- rawData[,file_order]

   sampleNames(rawData)<- as.character(description[,2])
   experimentFactor <- factor(description[,3])

   # if required reorder the arrays according to group levels in order to keep
   # groups together in all plots
   if(reOrder) {
     rawData <- rawData[,order(experimentFactor)]
     experimentFactor <- experimentFactor[order(experimentFactor)]
   }
 } else {
   # No information: arrays will be computed/colored independently
   sampleNames(rawData) <- as.character(sampleNames(rawData))
   experimentFactor <- factor(rep(1, length(sampleNames(rawData))))
   description <- cbind(sampleNames(rawData),sampleNames(rawData),
     experimentFactor)
   colnames(description) <- c("ArrayDataFile","SourceName","FactorValue")
 }

 if (aType == "PMonly") {
   if (normMeth == "MAS5") {
     warning("MAS5 cannot be applied to PMonly arrays. Changed MAS5 to PLIER")
     normMeth <- "PLIER"
   }
   if (normMeth == "GCRMA") {
     warning("GCRMA cannot be applied to PMonly arrays. Changed GCRMA to RMA")
     normMeth <- "RMA"
   }
 }

 if(!is.null(normMeth)) {
   if(customCDF) {
     if(species=="") {
       warning("Species has not been set and custom cdf requested, attempting to deduce species for chip type")
       species <- deduceSpecies(rawData@annotation)
     }
	 if(species!=""){
		 normData <- normalizeData(rawData,normMeth,perGroup=(normOption1=="group"),
		   experimentFactor, aType=aType, customCDF, species, CDFtype,WIDTH=WIDTH,
		   HEIGHT=HEIGHT)
	 }else{
		 warning("Could not define species; the CDF will not be changed")
		 normData <- normalizeData(rawData,normMeth,perGroup=(normOption1=="group"),
		   experimentFactor, aType=aType, customCDF,WIDTH=WIDTH,HEIGHT=HEIGHT)
	 }

   } else {
     normData <- normalizeData(rawData,normMeth,perGroup=(normOption1=="group"),
	   experimentFactor, aType=aType, customCDF,WIDTH=WIDTH,HEIGHT=HEIGHT)
     saveRDS(normData, file = paste0(DATA.DIR,"my_data.rds"))  
     print("Saving normalized data table")

     normDataTable <- createNormDataTable(normData,
                                          customCDF=(sum(featureNames(normData)!=featureNames(rawData)[1:length(featureNames(normData))])>0),
                                          species,
                                          CDFtype)

   #output normalised expression data to file
   refName <- sub("(_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}_\\d{2})", "", refName)
   normFileName <- paste(normMeth,"NormData_",refName,".txt",sep="")
   print(paste("Normalized data table:", normFileName))
   write.table(normDataTable, normFileName, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
   }
 }

 if(PMAcalls) {
   if(customCDF) {
     if(species=="") {
       warning("Species has not been set and custom cdf requested, attempting to deduce species for chip type")
       species <- deduceSpecies(rawData@annotation)
     }
 	 if(species!=""){
		 PMAtable <- computePMAtable(rawData,customCDF,species,CDFtype)
	 }else{
		 warning("Could not define species; the CDF will not be changed")
		 PMAtable <- computePMAtable(rawData,customCDF)
	 }
   } else {
     PMAtable <- computePMAtable(rawData,customCDF)
   }
   if(!is.null(PMAtable)) {
     write.table(PMAtable, "PMAtable.txt", sep="\t", row.names=FALSE,
	   col.names=TRUE, quote=FALSE)
   }
 }
print("DONE")
