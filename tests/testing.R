sourceEntireFolder <- function(folderName, verbose=FALSE, showWarnings=TRUE) {
  files <- list.files(folderName, full.names=TRUE)
  # Grab only R files
  files <- files[ grepl("\\.[rR]$", files) ]

  if (!length(files) && showWarnings)
    warning("No R files in ", folderName)

  for (f in files) {
    if (verbose)
      cat("sourcing: ", f, "\n")
    ## TODO:  add caught whether error or not and return that
    try(source(f, local=FALSE, echo=FALSE), silent=!verbose)
  }
  return(invisible(NULL))
}

sourceEntireFolder('R\\')

model <- CreateModel(n_segments=3, n=1000, p = 100, RandomNetwork)
x <- SimulateFromModel(model)

delete_values <- function(x,p){
  tmp <- matrix(runif(length(x)), nrow=nrow(x))
  x[tmp <= p] <- NA
  x
}

res <- list()
x_missing <- list()

pvalues <- c(0.1,0.2,0.3,0.4,0.5)

for (p in c(0,pvalues)){
  x_missing[[as.character(p)]] <- delete_values(x,p)
}

method <- c('averaging', 'pairwise_covariance', 'interpolation')

for (p in pvalues){
  for (mth in method){
    res[[method]][[as.character(p)]] <- hdcd(x_missing[[as.character(0)]],lambda=0.1, gamma=0.1, optimizer='section_search', parallel=F,
                                             FUN=NULL, method='glasso'
                                               #NA_SegmentLoss(mth)
                                             , n_folds=5, delta=0.1)
  }
}
