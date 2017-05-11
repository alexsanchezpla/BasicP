## ----echo=FALSE, include=FALSE-------------------------------------------
require(devtools)
install("/home/mguerrero/Dropbox/Practiques Grau (Biotec-UVic)-Mercedes_Guerrero/BasicP")
require(affy)
require(oligo) #ExonStudy

## ----message=FALSE-------------------------------------------------------
readCELS <- TRUE 
my.targets <- "./celfiles/targets.txt"
targets<- read.table("./celfiles/targets.txt", head=TRUE, sep="\t", row.names = 1)
my.fileNames <-paste("./celfiles/",rownames(targets),sep="")
rawDataFileName <- "rawData.Rda"
my.outputDir <- "."
isExonStudy <- FALSE
orgPackage <- "org.Hs.eg" # Necesita este parametro para funciones internas
rawData <- BasicP::readOrLoad.RawData(readCELS = readCELS, phenoDat = my.targets,
                      fileNames = my.fileNames, dataFName =rawDataFileName,
                      outputDir = my.outputDir, exonSt = isExonStudy)

## ----message=FALSE-------------------------------------------------------
loadAnnotations <- FALSE
chipPackAvailable <- TRUE
platformDesignPackAvailable <- FALSE
chipPackage <- "hgu133a2"
platformDesignPackage <- NULL
outputDir <- "./ResultsDir"
annotationsFileName <- "Annotations"
entrezTableFileName <-"Entrezs.Rda"
symbolsTableFileName <-"Symbols.Rda"
controlsTableFileName <- "controls.Rda"

anotacions <- BasicP::createOrLoadAnnotations (loadAnnotations= loadAnnotations, chipPackAvailable = chipPackAvailable, platformDesignPackAvailable = platformDesignPackAvailable,chipPackage = chipPackage, platformDesignPackage = platformDesignPackage, outputDir = outputDir,annotationsFileName = annotationsFileName,entrezTableFileName = entrezTableFileName, symbolsTableFileName = symbolsTableFileName, controlsTableFileName = controlsTableFileName)


## ----message=FALSE , warning=FALSE---------------------------------------
load("rawData.Rda")
rawData <- my.raw
normMethod <- "RMA"
my.targets <- read.AnnotatedDataFrame("./celfiles/targets.txt", header = TRUE, row.names = 1)
celFilesDir <-"./celfiles"
loadFile <- FALSE 
normalized.eset.FileName <-  "normalizedData.Rda"   
outputDir <- "./ResultsDir"
exonStudy <- FALSE

eset_norm <- BasicP::normalization(my.data = rawData, method = normMethod, targetsinfo = my.targets, inputDir = celFilesDir, loadFile = loadFile , normalizedFName = normalized.eset.FileName, outputDir = outputDir, exonSt = exonStudy)

## ----message=FALSE,warning=FALSE-----------------------------------------
load("./ResultsDir/normalizedData.Rda")
repes <- duplicated(exprs(my.norm), MARGIN=1)
exprs(my.norm) <- exprs(my.norm)[!repes,]
eset_norm <- my.norm
my.colors <- rainbow(length(sampleNames(eset_norm)))
my.names <- pData(eset_norm)$ShortName
myCex<- 0.8
dim3 <- FALSE
fileName <- "NormalizedPlots.pdf"
outputDir <- "./ResultsDir"
PCAPlots <- TRUE
csv <- "csv2"

BasicP::normplots2File(my.data = eset_norm, sampleNames = my.names, my.colors = my.colors, my.groups = pData(eset_norm)$Group, my.method = "average",my.cex = myCex ,posText = 2, dim3 = FALSE,fileName = fileName, outputDir = outputDir,PCAPlots = TRUE, csv = fileType)


## ------------------------------------------------------------------------
load("./ResultsDir/normalizedData.Rda")
repes <- duplicated(exprs(my.norm), MARGIN=1)
exprs(my.norm) <- exprs(my.norm)[!repes,]
eset_norm <- my.norm
load("./ResultsDir/controls.Rda")
removeNAs <- TRUE
load("./ResultsDir/Entrezs.Rda")
entrezs <- entrezTable
SignalFilter <- TRUE
signalThreshold <- 50
signalFilter.Function <- "filt.by.Signal"
signalThreshold.as.percentage <- TRUE
VarFilter <- TRUE 
variabilityThreshold <- 50
variability.Function <- "sdf"
variabilityThreshold.as.percentage <- TRUE
pairing.Function <- NULL
my.targets <-read.AnnotatedDataFrame("./celfiles/targets.txt", header = TRUE, row.names = 1)
doReport <- TRUE
outputDir <- "./ResultsDir"
FilteringReportFileName <- "FilteringReport.txt"

exprs.filtered <- BasicP::filterData(expres = exprs(eset_norm),controls = names(controlsTable),removeNAs = TRUE, entrezs = entrezs ,bySignal = SignalFilter,signalThr = signalThreshold, grups = pData(eset_norm)$grupo, sigFun.Name = signalFilter.Function, sigThr.as.perc = signalThreshold.as.percentage, byVar = VarFilter, variabilityThr = variabilityThreshold, varFun.Name = variability.Function,
varThr.as.perc = variabilityThreshold.as.percentage, pairingFun.Name = pairing.Function, targets = my.targets, doReport = doReport, outputDir = outputDir, filteringReportFName = FilteringReportFileName)

## ------------------------------------------------------------------------
load("./ResultsDir/normalizedData.Rda")
repes <- duplicated(exprs(my.norm), MARGIN=1)
exprs(my.norm) <- exprs(my.norm)[!repes,]
eset_norm <- my.norm
normalized.all.FileName <- "normalized.all"
fileType <-"csv2"
symbolsTable <- load("./ResultsDir/Symbols.Rda")
expres.all.FileName <- "expres.Rda"
linksFileName <- "Links.txt"
outputDir <- "./ResultsDir"

BasicP::saveData(expres = exprs(eset_norm), expres.csv.FileName = normalized.all.FileName, csvType=fileType, description = "Normalized values for all genes", anotPackage = NULL, symbolsVector = symbolsTable, SYMBOL = "SYMBOL", expres.bin.FileName = expres.all.FileName, linksFile = linksFileName, outputDir = outputDir)


