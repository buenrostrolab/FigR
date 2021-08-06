# Misc helper functions

### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University

#' Split and fetch string subparts
#'
#' Function to split and fetch string components based on a specific delimiter
#'
#' @param vec character vector to apply splitting to
#' @param delim character specifying the delimitor to separate string by
#' @param part integer specifying the part to fetch after splitting
#' @return substring based on splitting criteria
#' @examples
#' my_string <- c("This_is_random_string.111","Yes_it_is.222")
#' # Get first part splitting on _
#' splitAndFetch(vec=my_string,delim="_",part=1)
#' # Get first two parts splitting on _
#' splitAndFetch(vec=my_string,delim="_",part=1:2)
#' # Get second part splitting on .
#' splitAndFetch(vec=my_string,delim=".",part=2)
#'@export
#'@author Vinay Kartha
splitAndFetch <- function(vec,
                          delim,
                          part){
  if(length(part)==1){
    sapply(strsplit(as.character(vec),delim,fixed=TRUE),"[[",part) } else {
      sapply(strsplit(as.character(vec),delim,fixed = TRUE),function(x) paste(x[part],collapse = delim))
    }
}

#' Center single-cell ATAC-seq count data
#'
#'Function to center raw read counts for scATAC-seq data based on mean reads in features/rows (e.g. peaks) per cell
#'@param obj Either a \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} or \code{\link[Matrix]{dgeMatrix-class}} object of the single-cell ATAC data (peaks x cells) to center raw data for
#'@param doInChunks boolean value whether or not to score cells in chunks (useful for large scATAC datasets of > 5,000 cells). If TRUE, cells are centered sequentially in chunks of size=1000 cells at a time to avoid memory limitations, and eventually merged (column-wise). Overriden and set to TRUE if > 10,000 cells in dataset
#'@param chunkSize numeric specifying the number of cells to perform centering for at once, if running in chunks (to save memory). Default is 1000
#'@return Either a \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} or \code{\link[Matrix]{dgeMatrix-class}} object with counts centered by mean reads in peaks per cell
#'@import SummarizedExperiment Matrix
#'@export
#'@author Vinay Kartha
centerCounts <- function(obj,
                         doInChunks=TRUE,
                         chunkSize=1000){
  if(!class(obj) %in% c("SummarizedExperiment","RangedSummarizedExperiment","dgCMatrix","dgeMatrix","Matrix"))
    stop("Supplied object must be either of class SummarizedExperiment or sparse Matrix ..\n")

  if(ncol(obj) > 10000)
    doInChunks <- TRUE

  if(doInChunks){
    cat("Centering counts for cells sequentially in groups of size ",
        chunkSize, " ..\n\n")
    starts <- seq(1,ncol(obj),chunkSize)
  } else{
    starts <- 1
  }

  counts.l <- list()

  for(i in 1:length(starts)){
    beginning <- starts[i]
    if(i==length(starts)) # If it's the last one
    {
      ending <- ncol(obj)
    } else {
      ending <- starts[i]+chunkSize-1
    }

    cat("Computing centered counts for cells: ",beginning," to ", ending,"..\n")

    if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
      m <- SummarizedExperiment::assay(obj[, beginning:ending])} else {
        m <- obj[,beginning:ending] # Assumes Matrix format
      }
    cellMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    # Center cell counts based on its mean RIP count
    cCounts <- Matrix::t(Matrix::t(m)/cellMeans)

    counts.l[[i]] <- cCounts

    gc()
  }

  cat("Merging results..\n")
  centered.counts <- do.call("cbind",counts.l)
  cat("Done!\n")

  if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  } else {
    return(centered.counts)
  }
}


#' Smooth single cell matrix using cell kNNs
#'
#' Method to smooth sparse scores (e.g. RNA expression, ATAC gene scores, or DORC scores) per cell per feature using cell K nearest neighbors (NNs). Useful for visualizing scores in single cells for tSNE and UMAP plots.
#'
#'@param NNmat matrix of K Nearest neighbor indices for each cell (row), where each column is one of K nearest neighbors per cell. See the \code{\link[FNN]{knn.index}} function for obtaining these K nearest neighbor maps. Rownames must be valid cell IDs.
#'@param mat matrix of single cell scores. Rownames must be valid gene symbols, and column names must be valid cell IDs matching the rownames in the cell NN matrix.
#'@param geneList vector of valid gene names to compute smoothed TSS accessibility scores for. If geneList is NULL, all genes (rows) in the TSS matrix are used
#'@param barcodesList vector of valid barcode IDs to subset data using and compute smoothed TSS accessibility scores for. If barcodesList is NULL, all cell barcodes (columns) in the TSS matrix are used
#'@param nCores integer specifying the number of cores to use, if running in parallel. Default is 1 (i.e. no parallelization)
#'@importFrom parallel mclapply
#'@import Matrix doParallel
#'@return a matrix of scores for each gene TSS and each cell, smoothed over its K nearest neighbors
#'@export
#'@author Vinay Kartha
smoothScoresNN <- function(NNmat,
                           mat,
                           geneList = NULL,
                           barcodesList=NULL,
                           nCores = 1)
{
  if (is.null(rownames(NNmat)))
    stop("NN matrix has to have matching cell IDs as rownames\n")
  if (!all.equal(rownames(NNmat), colnames(mat)))
    stop("Nearest-neighbor matrix and cell data matrix don't have matching cells barcodes ..\n")
  cat("Number of cells in supplied matrix: ", ncol(mat),
      "\n")
  cat("Number of genes in supplied matrix: ", nrow(mat),
      "\n")
  cat("Number of nearest neighbors being used per cell for smoothing: ",
      ncol(NNmat), "\n")
  if (!is.null(geneList)) {
    if (!(all(geneList %in% rownames(mat)))) {
      cat("One or more of the gene names supplied is not present in the matrix provided: \n")
      cat(geneList[!geneList %in% rownames(mat)], sep = ", ")
      cat("\n")
      stop()
    }
    cat("Running smoothing only on genes:", geneList,
        sep = "\n")
    cat("........\n")
    mat <- mat[rownames(mat) %in% geneList, ]
  }
  else {
    if(nrow(mat) > 10000){
      cat("Running smoothing for all genes in matrix! (n = ",
          nrow(mat), ") This is bound to take more time than querying specific markers ..\n",
          sep = "")
    }
  }
  opts <- list()
  pb <- txtProgressBar(min = 0, max = ncol(mat), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()


  cl <- parallel::makeCluster(nCores)
  doSNOW::registerDoSNOW(cl)

  if(!is.null(barcodesList)){
    cat("Subsetting to ",length(barcodesList)," barcodes in dataset..\n")
    NNmat <- NNmat[barcodesList,]
  }
  cat("Running in parallel using ", nCores, "cores ..\n")
  matL <- foreach::foreach(x=1:nrow(NNmat),.options.snow = opts,.packages = c("Matrix","data.table","dplyr")) %dopar% {
    smoothedScore <- data.table(Matrix::rowMeans(mat[, NNmat[x,]]))
    rownames(smoothedScore) <- rownames(mat)
    colnames(smoothedScore) <- rownames(NNmat)[x]
    smoothedScore
  }

  parallel::stopCluster(cl)

  close(pb)
  cat("Merging results ..\n")
  smoothedMat <- dplyr::bind_cols(matL) %>% data.matrix() %>% Matrix(sparse=TRUE)
  rownames(smoothedMat) <- rownames(mat)

  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed),
            "\n"))

  return(smoothedMat)
}

#' Peak-based TF motif enrichment Z-test
#'
#' Wrapper function to assess a specific set of peaks for enrichment with respect to TF binding motifs using a Z test after adjusting for GC bias and Tn5 accessibility
#'
#' @param peakSet vector of peak indices you wish to test for motif enrichment
#' @param bgPeaks background peaks matrix returned by \code{\link[chromVAR]{getBackgroundPeaks}} function in chromVAR
#' @param tfMat binary match \code{\link[Matrix]{sparseMatrix}} object of (reference) peaks by TF motifs returned by \code{\link[motifmatchr]{matchMotifs}} function
#' @return data frame of motif enrichment results rank sorted by significance and frequency of motifs in the peakset, respectively
#' @import chromVAR SummarizedExperiment dplyr
#' @export
#' @author Vinay Kartha, Samuel Rose, and Jason D. Buenrostro
motifPeakZtest <- function(peakSet,
                           bgPeaks,
                           tfMat
) {

  if(nrow(tfMat)!=nrow(bgPeaks))
    stop("Reference peak set used for TF and background peaks matrix must match..\n")

  if(!all(peakSet %in% 1:nrow(bgPeaks)))
    stop("One or more of the provided peak indices are out of the background peak set range ..\n")


  # Filter out TFs whose motifs that overlap no peaks (if they exist)
  tfMat <- tfMat[,!Matrix::colSums(tfMat)!=0]

  # get selected peak motif frequencies
  cat("Getting selected peak motif frequencies ..\n")

  # get frequency of motifs in test set (observed)
  p.tab <- Matrix::colMeans(tfMat[peakSet, ])

  # get the background frequency in peak sets of the same size
  cat("Getting background peak motif frequencies ..\n")
  # extract relevant rows (i.e. peakset being tested) from background peak matrix
  bg.f <- as.matrix(bgPeaks[peakSet, ])

  # calculate (background) motif frequencies in each iteration of background peaks corresponding to peakset
  bg.tab <- apply(bg.f[, c(1:ncol(bgPeaks))], 2, function(bg_iter) {

    b.i <- Matrix::colMeans(tfMat[bg_iter, ])
    return(b.i)

  })

  cat("Calculating empirical p values and z score p values ..\n")

  # loop over each motif and generate enrichment statistics compared to background
  m.p <- dplyr::bind_rows(lapply(names(p.tab), function(motif) {

    # calculate sd and mean frequencies for bg and selected peaks
    s <- sd(bg.tab[motif, ])
    bg_freq <- mean(bg.tab[motif, ])

    z_score <- (p.tab[motif] - bg_freq) / s

    if(is.nan(z_score))
      z_score <- 0

    # generate data.frame object of relevant statistics
    d <- data.frame(
      motifID = motif,
      gene = extractTFNames(motif),
      motif_obs_freq = p.tab[motif],
      motif_bg_freq = mean(bg.tab[motif, ]),
      motif_counts = p.tab[motif] * length(peakSet),
      emp_pval = 1 - (sum(bg.tab[motif, ] < p.tab[motif]) / ncol(bg.tab)),
      z_test = z_score,
      pval.z = 2 * pnorm(-abs(z_score)),
      signed.log10p = -log10(2 * pnorm(-abs(z_score))) * sign(z_score)
    )
    return(d)
  }))
  # sort by enrichment pval, motif observed frequency
  m.p <- dplyr::arrange(m.p,pval.z, motif_obs_freq)
  # return df of enrichment scores
  return(m.p)
}

#' Fetch TF names
#'
#' This function fetches TF names from full motif IDs (see \url{https://github.com/GreenleafLab/chromVARmotifs})
#'
#'@param motifIDs vector of full motif IDs that are part of the original human or mouse pwms matrix (see \url{https://github.com/GreenleafLab/chromVARmotifs})
#'@return a vector of TF names extracted from the full motif ID (if older format with '_'), or returned as is otherwise
#'@export
#'@author Vinay Kartha
extractTFNames <- function(motifIDs){
  if(all(grepl("_",motifIDs,fixed = TRUE))){
    sapply(strsplit(sapply(strsplit(motifIDs,"_LINE.",fixed=FALSE),"[[",2),"_",fixed=FALSE),"[[",2)
  } else {
    message("One or more provided motif IDs do not contain any '_' characters .. returning IDs as is")
    motifIDs
  }
}

#'@export
#'@author Vinay Kartha
all.unique <- function(x){
  length(x)==length(unique(x))
}
