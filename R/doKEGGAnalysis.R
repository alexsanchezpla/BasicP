#####################################################
loadFromFile <-function (fileName, pos=1){
  tempEnv =new("environment")
  load (fileName, envir = tempEnv)
  varNames <-ls(tempEnv)
  myVarName <- varNames[pos]
  load (fileName)
  myVar <- eval(parse(text=myVarName))
  return(myVar)
}
#####################################################
#'doKEGGAnalysis
#'
#' Function to interate with KEGGAnalysis
#' @param KEGGPar List object that contains the parameters.
#' @return KEGGResult 
#' @examples  for(i in 1:length(KEGGParsList))
#'  {
#'  KEGGList <- doKEGGAnalysis(KEGGParsList[i])
#'  }
#' @export
#' 
doKEGGAnalysis <- function(KEGGPar)
{
  
  p <- KEGGPar[[1]]
  if(!is.null(p$my.IDs)){
    EntrezIDs <-  eval(parse(text = p$my.IDs)) #  Seran el EntrezTable
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
  
  KEGGResult <- KEGGAnalysis(fitMain = fitMain,
                             whichContrasts = p$whichContrasts, 
                             comparison.Name = p$comparisonName, 
                             outputDir = p$outputDir, 
                             anotPackage = orgPackage,
                             organisme = organisme,
                             my.IDs = EntrezIDs, # era p$my.IDs
                             addGeneNames = p$addGeneNames,
                             fileOfLinks = p$fileOfLinks,   
                             cutoffMethod = p$cutoffMethod, 
                             P.Value.cutoff = p$P.Value.cutoff,
                             pval = p$pvalKEGGterms,
                             thrLogFC = p$minLogFC,
                             minNumGens = p$minNumGens)
  
  return(KEGGResult)
}
