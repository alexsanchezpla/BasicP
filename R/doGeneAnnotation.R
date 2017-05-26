######################################################
loadFromFile <-function (fileName, pos=1){
  tempEnv =new("environment")
  load (fileName, envir = tempEnv)
  varNames <-ls(tempEnv)
  myVarName <- varNames[pos]
  load (fileName)
  myVar <- eval(parse(text=myVarName))
  return(myVar)
}
######################################################


#'doGeneAnnotation
#'
#' Function that annotates a data set.
#'
#' @param anotList List object that contains the parameters
#' @return genesAnnotated (??)
#' @export
#' @exemples doGeneAnnotation(anotList)
#'
#'
doGeneAnnotation <- function(AnotList)
{

  p <- AnotList
  
  if(!is.null(p$my.IDs))
  {
    EntrezIDs <-  eval(parse(text = p$my.IDs))
  }

  if (!is.null(p$fitFileName))
  {
    fitMain <- loadFromFile(file.path(p$outputDir, p$fitFileName))
  }else{
    if (!is.null(p$fitMain))
    {
      fitMain <- eval(parse(text = p$fitMain)) # Posar-hi un tryCatch per poder sortir si dona error!!!
    }else{
      stop("Error, cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
    }
  }

  genes2annotate <- EntrezIDs[unique(rownames(fitMain$p.value))]

  genesAnnotated <- GeneAnnotation(egIDs = genes2annotate,
                                   anotPackage = orgPackage,
                                   toHTML = p$toHTML,
                                   outputDir = p$outputDir,
                                   filename = p$anotFilename,
                                   myTitle = p$titleAnotations,
                                   specie = organisme,
                                   info2show = p$info2show,
                                   linksFile = p$linksFile,
                                   maxGenes = p$numGenesPerPage)
  return(genesAnnotated)
}
