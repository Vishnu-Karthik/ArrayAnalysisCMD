    library(methods)
    library(argparser)
    options_ <- arg_parser("")
    options_ <- add_argument(options_, "data_directory", help="Full path to directory with .cel files")
    options_ <- add_argument(options_, "output_directory", help="Full path to the output directory")
    options_ <- add_argument(options_, "description_file", help="Full path to description file")
    options_ <- add_argument(options_, "--all", help="Generate all possible plots", flag=TRUE)
    # Sample quality
    options_ <- add_argument(options_, "--samplePrep", help="Boolean for sample prep plot", flag=TRUE)
    options_ <- add_argument(options_, "--ratio", help="Boolean for 3'/5' ratio plots", flag=TRUE)
    options_ <- add_argument(options_, "--rnadeg", help="Boolean for RNA degradation plot", flag=TRUE)
    # Hybridization and signal quality
    options_ <- add_argument(options_, "--hybrid", help="Boolean for spike-in controls", flag=TRUE)
    options_ <- add_argument(options_, "--bgPlot", help="Boolean for background intensities", flag=TRUE)
    options_ <- add_argument(options_, "--percPres", help="Boolean for 3'/5' for b-actin and GAPDH", flag=TRUE)
    options_ <- add_argument(options_, "--posnegDist", help="Boolean for positive and negative controls distribution ", flag=TRUE)
    options_ <- add_argument(options_, "--control", help="Boolean for plots of the AFFX controls on the arrays", flag=TRUE)
    # Signal comparability and bias diagnostic
    options_ <- add_argument(options_, "--scaleFact", help="Boolean for scale factor", flag=TRUE)
    options_ <- add_argument(options_, "--box", help="Boolean for boxplot of raw intensities", flag=TRUE)
    options_ <- add_argument(options_, "--density", help="Boolean for density histogram of raw intensities", flag=TRUE)
    options_ <- add_argument(options_, "--MARaw", help="Boolean for MA-plot", flag=TRUE)
    options_ <- add_argument(options_, "--layout", help="Boolean for plot of array layout", flag=TRUE)
    options_ <- add_argument(options_, "--posnegCOI", help="Boolean for positive and negative control plots", flag=TRUE)
    options_ <- add_argument(options_, "--PLMimage", help="Boolean for 2D PLM plots", flag=TRUE)
    options_ <- add_argument(options_, "--spatialImage", help="Boolean for 2D images", flag=TRUE)
    options_ <- add_argument(options_, "--Nuse", help="Boolean for Nuse plot", flag=TRUE)
    options_ <- add_argument(options_, "--Rle", help="Boolean for RLE plot", flag=TRUE)
    # Array correlation
    options_ <- add_argument(options_, "--correlRaw", help="Boolean for correlation plot", flag=TRUE)
    options_ <- add_argument(options_, "--PCARaw", help="Boolean for PCA analysis plot", flag=TRUE)
    options_ <- add_argument(options_, "--clusterRaw", help="Boolean for hierarchical clustering", flag=TRUE)
    # Other options
    options_ <- add_argument(options_, "--Reorder", help="Reorder arrays by group", flag=TRUE)
    options_ <- add_argument(options_, "--perGroup", help="Generate MAplot per group. Default is for dataset", flag=TRUE)
    argv <- parse_args(options_)

    DATA.DIR <-argv$data_directory
    WORK.DIR <-argv$output_directory
    description_file <- argv$description_file
    SCRIPT.DIR <- "~/projects/sols/masters/4th_sem/scripts/"

    if (isTRUE(argv$all)){
    # Sample quality
      samplePrep <- TRUE
      ratio <- TRUE
      degPlot <-TRUE
    # Hybridization and signal quality
      hybrid <- TRUE
      bgPlot <- TRUE
      percPres <-TRUE
      posnegDistrib <- TRUE
      controlPlot <- TRUE
    # Signal comparability and bias diagnostic
      scaleFact <- TRUE
      boxplotRaw <- TRUE
      densityRaw <- TRUE
      MARaw <- TRUE
      layoutPlot <- TRUE
      posnegCOI <- TRUE
      PLMimage <- TRUE
      spatialImage <- TRUE
      Nuse <- TRUE
      Rle<- TRUE
    # Array correlation
      correlRaw <- TRUE
      PCARaw <- TRUE
      clusterRaw <- TRUE
    } else {
    # Sample quality
      samplePrep <- argv$samplePrep 
      ratio <- argv$ratio
      degPlot <- argv$rnadeg
    # Hybridization and signal quality
      hybrid <- argv$hybrid
      bgPlot <- argv$bgPlot
      percPres <- argv$percPres
      posnegDistrib <- argv$posnegDist
      controlPlot <- argv$control
    # Signal comparability and bias diagnostic
      scaleFact <- argv$scaleFact
      boxplotRaw <- argv$box
      densityRaw <- argv$density
      MARaw <- argv$MARaw
      layoutPlot <- argv$layout
      posnegCOI <- argv$posnegCOI
      PLMimage<- argv$PLMimage
      spatialImage <- argv$spatialImage
      Nuse <- argv$Nuse
      Rle<- argv$Rle
    # Array correlation
      correlRaw <- argv$correl
      PCARaw <- argv$PCA
      clusterRaw <- argv$cluster
    }

    # Other options
    reOrder <- argv$Reorder
    perGroup <- argv$perGroup

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
   require("gplots", quietly = TRUE) #heatmap.2 functions
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

 colList <- colorsByFactor(experimentFactor)
 plotColors <- colList$plotColors
 legendColors <- colList$legendColors
 rm(colList)

 plotSymbols <- 18-as.numeric(experimentFactor)
 legendSymbols <- sort(plotSymbols, decreasing=TRUE)

 WIDTH <- 1000
 HEIGHT <- 1414
 POINTSIZE <- 24
 if(!exists("maxArray")) maxArray <- 41

 #create a cover sheet for the report to be created later
 #and create a page indicating the naming and grouping used
 coverAndKeyPlot(description, refName,WIDTH=WIDTH,HEIGHT=HEIGHT)

 #create a table with several QC indicators
 if(samplePrep || ratio || hybrid || percPres || bgPlot || scaleFact) {

   # The indicators are calculated only for PM-MM arrays as the calculation
   # based on MAS5 does not work for PM-only arrays

   quality <- NULL
   try(quality <- qc(rawData),TRUE) # calculate Affymetrix quality data for PMMM
   if(is.null(quality)) {
     warning("Plots based on the simpleaffy qc function cannot be created for this chip type")
   }

   if(samplePrep) {
     # find the data
     try(yack <- yaqc(rawData),TRUE)
     if(exists("yack")) {
       spnames<-rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", # only 3'
       rownames(yack@morespikes), ignore.case = TRUE),])
       sprep<-t(yack@morespikes[spnames,])
     } else {
       sprep <- NULL
       warning("Plots based on the yaqc function cannot be created for this chip type")
     }

     try({calls<-detection.p.val(rawData)$call
     lys<-calls[rownames(calls)[grep("lys.*3",rownames(calls),ignore.case=TRUE)],]
     rm(calls)},TRUE)
     if(!exists("lys")) {
       lys <- NULL
       warning("Plots based on the detection.p.val function cannot be created for this chip type")
     }else{
		 if(length(lys) > length(sampleNames(rawData))) { lys<-lys[1,] }
     }
   }

   QCtablePlot(rawData,quality,sprep,lys,samplePrep=samplePrep,ratio=ratio,
       hybrid=hybrid,percPres=percPres,bgPlot=bgPlot,scaleFact=scaleFact,
	   WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
 }

 print("Graphs ready to be computed")

 if(samplePrep && !is.null(sprep) && !is.null(lys)) {
   print ("   plot sample prep controls"  )
   samplePrepPlot(rawData,sprep,lys,plotColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(ratio && !is.null(quality)) {
   print ("   plot beta-actin & GAPDH 3'/5' ratio")
   ratioPlot(rawData,quality=quality,experimentFactor,plotColors,legendColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(degPlot) {
   print ("   plot degradation plot"  )
   RNAdegPlot(rawData,plotColors=plotColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(hybrid && !is.null(quality)) {
   print ("   plot spike-in hybridization controls"  )
   hybridPlot(rawData,quality=quality,plotColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(bgPlot && !is.null(quality)) {
   print ("   plot background intensities"  )
   backgroundPlot(rawData,quality=quality,experimentFactor,plotColors,legendColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(percPres && !is.null(quality)) {
   print ("   plot percent present"  )
   percPresPlot(rawData,quality=quality,experimentFactor,plotColors,legendColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(posnegDistrib) {
   print ("   plot pos & neg control distribution"  )
   PNdistrPlot(rawData,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
 }

 if(controlPlot) {
   print ("   plot control profiles and/or boxplots")
   controlPlots(rawData,plotColors,experimentFactor,legendColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(scaleFact && !is.null(quality)) {
   print ("   plot scale factors")
   scaleFactPlot(rawData,quality=quality,experimentFactor,plotColors,
      legendColors,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
	  MAXARRAY=maxArray)
 }

 if(boxplotRaw){
   print ("   plot boxplot for raw intensities")
   boxplotFun(Data=rawData, experimentFactor, plotColors, legendColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(densityRaw){
   print ("   plot density histogram for raw intensities")
   densityFun(Data=rawData, plotColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
   #densityFunUnsmoothed(Data=rawData, plotColors,
   #  WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(MARaw){
   print ("   MA-plots for raw intensities")
   maFun(Data=rawData, experimentFactor, perGroup,
      aType=aType,WIDTH=WIDTH,HEIGHT=HEIGHT,MAXARRAY=maxArray)
 }

 if(layoutPlot) {
   print ("   plot array reference layout")
   plotArrayLayout(rawData,aType,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
 }

 if(posnegCOI){
   print ("   Pos/Neg COI")
   PNposPlot(rawData,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
 }

 # fit a probe level model on the raw data, used by nuse and rle plot as well
   rawData.pset <- NULL
   if(spatialImage || PLMimage || Nuse || Rle) {
   print ("   Fit a probe level model (PLM) on the raw data")
     rawData.pset <- fitPLM(rawData)
   }

 if(spatialImage) {
   print ("   2D virtual images")
   valtry<-try(spatialImages(rawData, Data.pset=rawData.pset, TRUE,FALSE,FALSE,FALSE,
	           WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE),
		       silent=TRUE)
   if(class(valtry)=="try-error") {
	 print("      Use array.image instead of spatialImages function")
	 if(length(sampleNames(rawData))>6){
		 # Usage of a median array is interesting when there are enough arrays
		 array.image(rawData,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
	 }else{
		 # Usage when few arrays in dataset (one page for 3 arrays -> max: 2 pages)
		 array.image(rawData,relative=FALSE,col.mod=4,symm=TRUE,WIDTH=WIDTH,
		   HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
	 }
   }
 }

 if(PLMimage) {
   print ("   Complete set of 2D PLM images")
   valtry<-try(spatialImages(rawData, Data.pset=rawData.pset, TRUE, TRUE, TRUE, TRUE,
	             WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray),
				 silent=TRUE)
   if(class(valtry)=="try-error") {
	 print("      Could not create the PLM images.")
   }
 }

 if(Nuse){
   print ("   NUSE boxplot")
   nuseFun(rawData, Data.pset=rawData.pset, experimentFactor, plotColors,
      legendColors,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
	  MAXARRAY=maxArray)
 }

 if(Rle){
   print ("   RLE boxplot")
   rleFun(rawData, Data.pset=rawData.pset, experimentFactor, plotColors,
      legendColors,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
	  MAXARRAY=maxArray)
 }

 if(correlRaw){
   print ("   Correlation plot of raw data")
   correlFun(Data=rawData, experimentFactor=experimentFactor, legendColors=legendColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }

 if(PCARaw){
   print("   PCA analysis of raw data")
   pcaFun(Data=rawData, experimentFactor=experimentFactor,
	 plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
	 legendSymbols=legendSymbols, namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&&
	 (length(sampleNames(rawData))<=(maxArray/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,
	 POINTSIZE=POINTSIZE)
 }

 if(clusterRaw){
   print ("   Hierarchical clustering of raw data")
   clusterFun(Data=rawData, experimentFactor=experimentFactor,
    clusterOption1=clusterOption1, clusterOption2=clusterOption2,
    plotColors=plotColors, legendColors=legendColors,
    plotSymbols=plotSymbols, legendSymbols=legendSymbols,
    WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
 }
print("DONE")
