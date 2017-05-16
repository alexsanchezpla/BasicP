###################################################
loadFromFile <-function (fileName, pos=1){
  tempEnv =new("environment")
  load (fileName, envir = tempEnv)
  varNames <-ls(tempEnv)
  myVarName <- varNames[pos]
  load (fileName)
  myVar <- eval(parse(text=myVarName))
  return(myVar)
}
###################################################
#' doMultCompAnalysis
#'
#' Function for make the comparative analysis.
#' @param mcPar List object that contains the parameters.
#' @return genelist
#' @export


doMultCompAnalysis <- function(mcPar)
{

  p <- mcPar[[1]]
  if (!is.null(p$fitFileName)){
    fitMain <- loadFromFile(file.path(p$outputDir,p$fitFileName))
  }else{
    if (!is.null(p$fitMain)) {
      fitMain <- eval(parse(text = p$fitMain))
      # Posar-hi un tryCatch per poder sortir si dona error!!!
    }else{
      stop("Error, Cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
    }
  }

  geneList <-  multipleComp(fitMain = fitMain,
                            whichContrasts = p$whichContrasts,
                            comparisonName = p$comparisonName,
                            titleText = p$titleText,
                            outputDir = p$outputDir,
                            anotPackage = p$anotPackage,
                            my.symbols = symbolsTable,
                            linksFile = p$fileOfLinks,
                            multCompMethod = p$multCompMethod,
                            adjustMethod = p$adjustMethod,
                            selectionType = p$selectionType,
                            P.Value.cutoff = p$P.Value.cutoff,
                            plotVenn = p$plotVenn,
                            colsVenn = p$colsVenn,
                            vennColors = p$vennColors,
                            cexVenn = p$cexVenn,
                            geneListFName=p$geneListFName,
                            csvType=p$csvType,
                            minLFC=p$minLogFC)

  return (geneList)
}
