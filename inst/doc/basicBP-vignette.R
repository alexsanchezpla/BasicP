## ----echo=FALSE, include=FALSE-------------------------------------------
require(devtools)
install("/home/mguerrero/Dropbox/Practiques-Grau_Biotec-UVic_Mercedes_Guerrero/BasicP")
require(affy)
require(oligo) #ExonStudy
require("limma")

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


## ------------------------------------------------------------------------
targets <- read.table(file ="./celfiles/targets.txt" , header = TRUE, sep = "\t")
column2design<- 5
lev <- targets[,column2design] 
design <- model.matrix( ~ 0 + lev)        
colnames(design) <- levels(lev)
rownames(design) <- targets$ShortName
numParameters <- ncol(design)
print(design)

cont.matrix<-makeContrasts(AvsB= (A -B), AvsL= (A-L), BvsL=(B-L),levels=design)
print(cont.matrix)


## ------------------------------------------------------------------------
load("./ResultsDir/exprs.filtered.Rda")
contrasts2test <- 1:ncol(cont.matrix)
anotPackage = NULL
comparison =  "Estudi"
outputDir = "./ResultsDir"
ENTREZIDs = "entrezTable"
SYMBOLIDs = "symbolsTable"
linksFile = "Links.txt"
fitFileName = "fit.Rda"
csvType= "csv"
rows2HTML= NULL
anotFileName <- "Annotations"
runMulticore = 0 
toTIFF= FALSE

fitMain <- 
  lmAnalysis(exprs.filtered = exprs.filtered, design = design, cont.matrix = cont.matrix, contrasts2test = contrasts2test, anotPackage = anotPackage, outputDir = outputDir, comparison = comparison, Expressions_And_Top = TRUE , showParams = FALSE , use.dupCorr = FALSE, block = NULL, nDups = 1 , ENTREZIDs = ENTREZIDs, SYMBOLIDs = SYMBOLIDs, linksFile = linksFile,fitFileName = fitFileName , csvType=csvType, rows2HTML = NULL, anotFileName = anotFileName)



## ----echo=FALSE, include=FALSE-------------------------------------------
#Function to make a list for make the dolmAnalysis.
add2parsList <- function(oneList,object)
{
  pos <- length(oneList) + 1
  oneList[[pos]] <- object
  names(oneList)[pos] <- oneList[[pos]]$comparisonName
  return(oneList)
}


## ------------------------------------------------------------------------
lmParsList <- list()
Estudi <- list(dades = NULL,
               expresFileName = "exprs.filtered.Rda",
               targets = targets,
               designMat = design,
               contMat = cont.matrix,
               whichContrasts = 1:ncol(cont.matrix),
               anotPack = NULL,
               outputDir = outputDir,
               ExpressionsAndTop = TRUE,
               showLmParams = FALSE, 
               use.dupCorr = FALSE,
               block = NULL,
               nDups = 1,
               comparisonName = comparison,  
               ENTREZIDs = "entrezTable",
               SYMBOLIDs = "symbolsTable",
               fileOfLinks = linksFile,
               fitFileName = fitFileName,
               csvType=csvType,
               rows2HTML = NULL,
               anotFilename = anotFileName
               )

lmParsList <- add2parsList(lmParsList, Estudi)
               
for(ix in 1:length(lmParsList))
{
  fit.Main   <- doLmAnalysis(lmParsList[ix])
}

