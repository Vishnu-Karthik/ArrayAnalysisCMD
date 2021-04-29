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
   normData <- readRDS(paste0(DATA.DIR,"my_data.rds"))

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

   if(boxplotNorm){
     print ("   plot boxplot for normalized intensities")
     boxplotFun(Data=normData, experimentFactor, plotColors, legendColors,
	   normMeth=normMeth,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
	   MAXARRAY=maxArray)
   }

   if(densityNorm){
     print ("   plot density histogram for normalized intensities")
     densityFun(Data=normData, plotColors, normMeth=normMeth,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
     #densityFunUnsmoothed(Data=normData, plotColors, normMeth=normMeth,
     #  WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
   }

   if(MANorm){
     print ("   MA-plots for normalized intensities")
     maFun(Data=normData, experimentFactor, perGroup=(MAOption1=="group"),
	  normMeth=normMeth,WIDTH=WIDTH,HEIGHT=HEIGHT,MAXARRAY=maxArray)
   }

   if(correlNorm){
     print ("   Correlation plot of normalized data")
     correlFun(Data=normData, normMeth=normMeth, experimentFactor=experimentFactor, legendColors=legendColors,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
   }

   if(PCANorm){
     print("   PCA graph for normalized data")
     pcaFun(Data=normData, experimentFactor=experimentFactor,normMeth=normMeth,
	   plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
	   legendSymbols=legendSymbols, namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&&
	   (length(sampleNames(rawData))<=(maxArray/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,
	   POINTSIZE=POINTSIZE)
   }

   if(clusterNorm){
     print ("   Hierarchical clustering of normalized data")
     clusterFun(Data=normData, experimentFactor=experimentFactor,
     clusterOption1=clusterOption1, clusterOption2=clusterOption2,
     normMeth=normMeth, plotColors = plotColors, legendColors = legendColors,
     plotSymbols=plotSymbols, legendSymbols=legendSymbols,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
   }
 }
