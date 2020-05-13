suppressMessages({
  library(RDAVIDWebService)  
})

setup_david <- function(email, url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService"){
  david <- DAVIDWebService$new(email=email, url=url)
  return(david)
}

run_david <- function(out_prefix, top_table, david_db, genelist_prefix, gene_col='ensembl',
                      symbol_col='hugo', logFcColStr="logFC", 
                      trend_col='regulation', trend_classes=c('UPREG', 'DOWNREG'),
                      num_clusters=5, num_attempts=10) {
  
  cat("Beginning a DAVID analysis:\n")
  cat("Attempting to connect to DAVID....")
  attempts = 0
  result1 = NULL
  result2 = NULL
  while (attempts < num_attempts) {
    try({
      if (is.null(result1)){
        result1 <- addList(david_db, 
                           top_table[top_table[[logFcColStr]]>0, gene_col],
                           idType="ENSEMBL_GENE_ID",
                           listName=paste(genelist_prefix, trend_classes[1], trend_col, sep='_'), 
                           listType="Gene")        
      }
      if (is.null(result2)){
        result2 <- addList(david_db, 
                          top_table[top_table[[logFcColStr]]<0, gene_col],
                          idType="ENSEMBL_GENE_ID",
                          listName=paste(genelist_prefix, trend_classes[2], trend_col, sep='_'), 
                          listType="Gene")
      }
    }, silent = TRUE)
    
    if (!is.null(result1) & !is.null(result2)){
      break
    } else {
      attempts <- attempts + 1
    }
  }
  
  if (is.null(result1) | is.null(result2)){
    stop(paste0('ERROR: Could not reach DAVID after ', num_attempts, ' attempts'))
  } else {
    cat("Done.\n")
  }
  
  cat("Adding annotations \n")
  result <- setAnnotationCategories(david_db, c("GOTERM_BP_ALL",  # Biological Process
                                                "GOTERM_MF_ALL",  # Molecular Function
                                                 "GOTERM_CC_ALL",  # Cellular Component
                                                 "UP_KEYWORDS"  # Uniprot Keywords
                                                ))

  geneLists = getGeneListNames(david_db)
  expectedGeneLists = c(paste(genelist_prefix, trend_classes[1], trend_col, sep="_"),
                        paste(genelist_prefix, trend_classes[2], trend_col, sep="_"))
  for (geneListPosition in c(1:length(geneLists))) {
    if (!geneLists[geneListPosition] %in% expectedGeneLists){
      next
    }
    cat(paste0("Writing/Plotting results for ",
               geneLists[geneListPosition], "\n"))
    result <- setCurrentGeneListPosition(david_db, geneListPosition)
    termCluster <- getClusterReport(david_db, type="Term")
    if (length(termCluster@cluster)==0){
      cat(paste0('WARNING: ', geneLists[geneListPosition], ' returned no hits.\n'))
      next
    }
    result <- getClusterReportFile(david_db, type="Term", 
                                   fileName=paste0(out_prefix, 
                                                   "_DAVID_",
                                                   geneLists[geneListPosition],
                                                   "_termClusterReport.tab"))
    
    for (clustNumber in c(1:num_clusters)){
        png(paste0(out_prefix, 
                   "_DAVID_",
                   trend_classes[geneListPosition], 
                   trend_col,
                   "_cluster",
                   clustNumber,
                   ".png"), width=500, height=500, units='px')
        print(plot2D(termCluster, clustNumber))
        dev.off()
    }
  }
}
