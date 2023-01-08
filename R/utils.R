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
  # Check if SE or Matrix
  # Added to avoid error https://github.com/buenrostrolab/FigR/issues/16
  if(any(sapply(c("SummarizedExperiment","RangedSummarizedExperiment"),function(x){ inherits(obj,x)}))){
    cat("SummarizedExperiment object input detected .. Centering counts under assay")
    isSE <- TRUE
  } else {
    if(any(sapply(c("dgCMatrix","dgeMatrix","Matrix"),function(x){ inherits(obj,x)}))){
      cat("Matrix object input detected")
      isSE <- FALSE
    } else {
      stop("Supplied object must be either of class SummarizedExperiment or Matrix ..\n")
    }
  }

  #if(!class(obj) %in% c("SummarizedExperiment","RangedSummarizedExperiment","dgCMatrix","dgeMatrix","Matrix"))


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

    #if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
    if(isSE){
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

  #if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
  if(isSE){
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
  stopifnot(all.equal(nrow(NNmat),ncol(mat)))

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

  if(Sys.info()['sysname'] %in% "Windows"){
    message("Windows OS detected .. Cannot support parallilzation using mclapply for mc.cores > 1")
    message("Using 1 core instead ..\n")
    nCores <- 1
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
  smoothedMat <- dplyr::bind_cols(matL) %>% data.matrix() %>% Matrix::Matrix(sparse=TRUE)
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
  tfMat <- tfMat[,Matrix::colSums(tfMat)!=0]

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
  # Note: this returned an error saying pval.z not found in the arrange, look into it
  #m.p <- dplyr::arrange(m.p,pval.z, motif_obs_freq)

  # return df of enrichment scores
  return(m.p)
}

#' Plot marker scores on single cell scatter plots
#'
#' Function to plot variable of interest (e.g. scATAC-seq gene score, motif score, scRNA-seq expression etc.) onto existing tSNE/UMAP coordinates of single cell clusters
#'
#'@param df data.frame object of 2-D coordinates, where each row is a cell and the first two columns corresponds to the first (x-axis) and second (y-axis) dimensions to use for the plot (e.g. UMAP coordinates). Can have additional cell metadata columns used for splitting (see splitBy parameter)
#'@param markerMat matrix of scores containing markers to visualize. Can be a sparse matrix. Must have rownames
#'@param markers character vector of markers to visualize (can list multiple). Must be one of the rownames in markerMat
#'@param pointSize integer specifying point size to use for ggplot (passed to the geom_point() layer in ggplot). Default is 0.5
#'@param plotClean boolean indicating whether or not to return a 'clean' plot without any axes
#'@param rasteRize boolean indicating whether to use raster point rendering (to avoid giant file sizes and too many point layers in PDFs). Requires ggrastr package if set to TRUE (Default)
#'@param splitBy character indicating a single variable to split cells by using facets when visualizing. Must be a valid column name in df
#'@param minCutoff cut off to use for capping  minimum values prior to visualization. Can either be a raw numeric value, a percentile value cut-off (0-1) or the number of standard deviations below the mean, used to set the minimum value in the dataset prior to plotting. For percentiles, must be a character with 'q' before the percentile cut-off value (e.g. 'q0.1' for 10\%ile minimum). For standard deviation-based capping, must be a character with 'sd' before the cut-off value (e.g. 'sd2' for using 2 standard deviations below the mean as the minimum value displayed)
#'@param maxCutoff cut off to use for capping  maximum values prior to visualization. Can either be a raw numeric value, a percentile value cut-off (0-1) or the number of standard deviations above the mean, used to set the minimum value in the dataset prior to plotting. For percentiles, must be a character with 'q' before the percentile cut-off value (e.g. 'q0.9' for 90\%ile maximum). For standard deviation-based capping, must be a character with 'sd' before the cut-off value (e.g. 'sd2' for using 2 standard deviations above the mean as the maximum value displayed)
#'@param colorPalette color palette name specification to use for coloring points. Default is "solar_extra". See \href{https://github.com/caleblareau/BuenColors}{BuenColors} package developed by Caleb Lareau for more palette options
#'@param legend.position character specifying where to place the legend. Valid options are "top", "bottom", "left", "right", or "none" (to remove the legend completely). Default is bottom
#'@param showDataRange boolean indicating whether or not to show the data numeric range in the legends (only applies if legends are being displayed). Default is TRUE. If set to FALSE, just shows as min and max
#'@param combine boolean indicating whether or not to merge plots and display as grid.Default is TRUE
#'@param ... extra parameters passed to cowplot's \code{\link[cowplot]{plot_grid}} function to control organization of plot grids. Useful parameters include nrow and ncol that control the number of rows and columns in the layout, respectively
#'@return If a single marker is queried, or combine is set to TRUE for multiple markers, a ggplot object will be printed to graphics window unless assigned to variable. If combine is set to FALSE, then a list object of individual plots per marker is returned
#'@import ggplot2 ggrastr BuenColors cowplot
#'@export
#'@author Vinay Kartha
plotMarker2D <- function(df, # data frame of tSNE or UMAP coordinates for single cells
                         markers, # name associated with markerScore
                         markerMat, # Matrix of scores
                         pointSize=0.5,
                         plotClean=TRUE,
                         rasteRize=TRUE,
                         splitBy=NULL,
                         minCutoff=NULL,
                         maxCutoff=NULL,
                         colorPalette="solar_extra",
                         legend.position="bottom",
                         showDataRange=TRUE,
                         combine=TRUE,
                         ... # Extra to cowplot
){

  # CHECK ON DATA FRAME
  stopifnot(class(df)=="data.frame")

  if(!is.null(splitBy))
    if(!all(splitBy %in% colnames(df)))
      stop("One or more of the splitBy variables is not a column in the provided dataframe\n")
  # Assume first two columns are UMAP/tSNE/PCA coords
  dimnames <- colnames(df)[1:2]


  # CHECK FEATURES

  if(is.null(markers)) {
    stop("Must provide at least 1 marker to plot")  } else {

      if (!(all(markers %in% rownames(markerMat)))) {
        message("One or more of the features supplied is not present in the matrix provided: \n")
        message(markers[!markers %in% rownames(markerMat)], sep = ", ")
        cat("\n")
        stop()
      }
    }

  markerMat <- markerMat[markers,,drop=FALSE]

  if(length(markers) > 20 & combine)
    stop("Plotting this many features at once can lead to squished graphics ..\n")


  # Score transformations
  transformScores <- function(markerScore,minCutoff,maxCutoff){

    if(!is.null(minCutoff)){
      # Quantile cut-off
      if(grepl(pattern = '^q.*',minCutoff)){
        quantile.cutMin <- as.numeric(strsplit(minCutoff,"q",fixed = TRUE)[[1]][2])
        if(!quantile.cutMin > 0 & quantile.cutMin < 1)
          stop("Must provide a fractional value to use for quantile cutoff. See help page for details\n")
        minCutoff <- quantile(markerScore,quantile.cutMin)}

      # SD cut-off
      if(grepl(pattern = '^sd.*',minCutoff)){
        sd.cutMin <- as.numeric(strsplit(minCutoff,"sd",fixed = TRUE)[[1]][2])
        if(sd.cutMin > 0 & sd.cutMin < 1)
          stop("SD cut-off must be an integer number representing the number of standard deviations away from the mean to use. See help page for details\n")
        minCutoff <- mean(markerScore) - (sd.cutMin * sd(markerScore))
      }
    } else {
      minCutoff <- min(markerScore)
    }

    if(!is.null(maxCutoff)) {
      # Quantile cut-off
      if(grepl(pattern = '^q.*',maxCutoff)){
        quantile.cutMax <- as.numeric(strsplit(maxCutoff,"q",fixed = TRUE)[[1]][2])
        if(!quantile.cutMax > 0 & quantile.cutMax < 1)
          stop("Must provide a fractional value to use for quantile cutoff. See help page for details\n")
        maxCutoff <- quantile(markerScore,quantile.cutMax)}

      # SD cut-off
      if(grepl(pattern = '^sd.*',maxCutoff)){
        sd.cutMax <- as.numeric(strsplit(maxCutoff,"sd",fixed = TRUE)[[1]][2])
        if(sd.cutMax > 0 & sd.cutMax < 1)
          stop("SD cut-off must be an integer number representing the number of standard deviations away from the mean to use. See help page for details\n")
        maxCutoff <- mean(markerScore) + (sd.cutMax * sd(markerScore))
      }
    } else {
      maxCutoff <- max(markerScore)
    }

    # Transform
    markerScore[markerScore <= minCutoff] <- minCutoff
    markerScore[markerScore >= maxCutoff] <- maxCutoff
    markerScore
  }

  gglist <- lapply(markers,function(i) {

    cat("Plotting ",i,"\n")

    # Transform
    mScore <- transformScores(markerScore = markerMat[i,],
                              minCutoff=minCutoff,
                              maxCutoff=maxCutoff)

    if (sum(mScore) == 0)
      warning("No counts detected for marker: ", i, "\n")
    if (anyNA(mScore))
      warning("NAs detected for marker: ", i, "\n")

    i <- gsub(pattern = "-",replacement = "",x = i)

    df[,i] <- mScore


    if(legend.position!="bottom"){
      if(!legend.position %in% c("top","right","left","none"))
        stop("Must specify a valid legend position speciciation .. See function options for details\n")
    }

    if(rasteRize){
      require(ggrastr)
      pointfun <- geom_point_rast
    } else {
      pointfun <- geom_point
    }

    if(!showDataRange){
      myColfun <- scale_color_gradientn(colours = BuenColors::jdb_palette(colorPalette),na.value = "orange",breaks=range(df[,i]),labels=c("min","max"))
    } else {
      myColfun <- scale_color_gradientn(colours = BuenColors::jdb_palette(colorPalette),na.value = "orange",breaks=scales::breaks_pretty(n=3))
    }

    g <- ggplot(df,aes_string(dimnames[1],dimnames[2],color=i)) +
      pointfun(size=pointSize) +
      theme_classic() +
      labs(title= i) +
      theme(plot.title = element_text(hjust=0.5),legend.position=legend.position,
            legend.key.size = unit(0.3, "cm"),legend.text = element_text(size=5),legend.title=element_blank())+
      #scale_color_gradientn(colours = BuenColors::jdb_palette(colorPalette),na.value = "orange",breaks=myLegBreaks)
      myColfun

    # Split?
    if(!is.null(splitBy))
      g <- g + facet_wrap(as.formula(paste("~", splitBy))) + theme(strip.background = element_blank())

    if(plotClean)
      g <- g + theme(axis.line = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank())

    g
  })


  if(length(gglist)==1){
    return(gglist[[1]])
  } else if(combine){
    cat("Merging plots ..\n")
    return(cowplot::plot_grid(plotlist = gglist,align="hv",...))
  } else {
    cat("Returning list of plots\n")
    names(gglist) <- markers
    return(gglist)
  }

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
    #message("One or more provided motif IDs do not contain any '_' characters .. returning IDs as is")
    motifIDs
  }
}

#'@export
#'@author Vinay Kartha
all.unique <- function(x){
  length(x)==length(unique(x))
}
