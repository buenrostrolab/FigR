# Script housing main DORC functions

### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University

chunkCore <- function(chunk,
                      A, # ATAC matrix
                      R, # RNA matrix
                      O, # Gene-Peak overlap pairing data.frame
                      met # Correlation method ("spearman" or "pearson")
){

  # Get indices of genes and peaks from overlap object for chunk
  # Assumes query hits are genes and subject hits are peaks in the overlap object
  geneIndices <- O$Gene[chunk[1]:chunk[2]]
  peakIndices <- O$Peak[chunk[1]:chunk[2]]

  pairnames <- cbind(rownames(A)[peakIndices],rownames(R)[geneIndices])

  uniquegenes <- unique(geneIndices)
  uniquepeaks <- unique(peakIndices)


  M1 <- as.matrix(Matrix::t(A[uniquepeaks,,drop=FALSE])) # In case only 1 match, keep matrix structure
  M2 <- as.matrix(Matrix::t(R[uniquegenes,,drop=FALSE])) # In case only 1 match, keep matrix structure

  # Peak x Gene correlation matrix, subset by peak-gene pair names to get corresponding correlation vector
  # NOTE: This correlation call fails if you have maps with just 1 gene / peak. This is unlikely for large chunk sizes
  cor(x = M1,y = M2,method = met)[pairnames]

}

PeakGeneCor <- function(ATAC, # Normalized reads in peaks counts (rownames should  be "Peak1","Peak2" etc.)
                        RNA, # Normalized gene expression counts
                        OV, # Gene TSS - Peak overlap pairs object (Genes: query, Peaks: subject)
                        ncores=4,
                        chunkSize=200,
                        metric="spearman",
                        bg=NULL){

  stopifnot(ncol(ATAC)==ncol(RNA))

  if(chunkSize > 1000)
    stop("Do not specify very large chunk sizes. Please use chunkSize < 1000")


  # Number of total gene-peak pairs to chunk up for parallelization
  n <- length(OV)
  starts <- seq(1, n, chunkSize)
  ends <- starts + chunkSize -1
  ends[length(ends)] <- n

  OVd <- OV %>% as.data.frame() %>% dplyr::rename("Gene"="queryHits","Peak"="subjectHits")

  chunkList <- mapply(c, starts, ends, SIMPLIFY = FALSE)

  time_elapsed <- Sys.time()

  cat("Running in parallel using ", ncores, "cores ..\n")

  cat("Computing observed correlations ..\n")

  corList <- pbmcapply::pbmclapply(X=chunkList,
                                   FUN=function(x) {FigR::chunkCore(chunk=x,A=ATAC,R=RNA,O=OVd,met=metric)},mc.cores = ncores)


  if(any(unlist(sapply(corList,is.null)))){
    message("One or more of the chunk processes failed unexpectedly (returned NULL) ..")
    message("Please check to see you have enough cores/memory allocated")
    message("Also make sure you have filtered down to non-zero peaks/genes")
  }

  OVd$rObs <- unlist(corList)


  cat("Finished!\n")

  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")

  if(!is.null(bg)){
    n_iter <- ncol(bg)
    cat("Computing background correlations ..\n")

    time_elapsed <- Sys.time()

    bgCor <- foreach(i=1:n_iter,.combine = 'cbind',
                     .export = c("chunkCore","t"),.packages = c("pbmcapply","FigR","Matrix")) %do% {
                       OVdBg <- OVd[,1:2] # Initialize gene-peak pairing to observed
                       OVdBg$Peak <- bg[OVdBg$Peak,i] # Swap actual peaks with bg peaks for given iteration in pairing
                       bgCorList <- pbmcapply::pbmclapply(X=chunkList,
                                                          FUN=function(x) {chunkCore(chunk=x,A=ATAC,R=RNA,O=OVdBg,met=metric)},mc.cores = ncores)
                       unlist(bgCorList) # Vector of background permuted correlation values for that set of background peaks
                     }


    if(sum(is.null(bgCor))!=0 | sum(is.na(bgCor))!=0)
      stop("One or more of the chunk processes failed unexpectedly (returned NULL) .. Please check to see you have enough cores/m
           emory allocated")

    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")

    colnames(bgCor) <- paste0("rBg",1:ncol(bgCor))
    OVd <- cbind(OVd,bgCor)
  }

  return(OVd)
}

#' Gene-Peak correlations used for calling domains of regulatory chromatin (DORCs)
#'
#'Function to compute correlation between RNA expression and peak accessibility for peaks falling within a window around each gene, across single cells, using background peak correlations for significance testing
#'@param ATAC.se SummarizedExperiment object of the scATAC-seq reads in peak counts.
#'@param RNAmat Matrix object of the normalized scRNA-seq counts. For example, this will be the normalized count data housed within the `@assays$RNA@data` field if processed using Seurat. Must have same number of cells (i.e. matched) with the scATAC-seq data
#'@param genome character specifying the reference genome build to use. Must be one of "hg19", "hg38" or "mm10", with no default
#'@param TSSwindow numeric specifying the window size (in base pairs) to pad around either side of each TSS, to fetch overlapping peaks. Default is 50 kb (so 100 kb window is drawn arund each TSS)
#'@param normalizeATACmat boolean indicating whether or not to normalize the counts present in the ATAC.se object prior to computing correlations. Default is TRUE (i.e. assumes peak counts in ATAC.se are raw)
#'@param nCores numeric indicating the number of cores to use if parallelizing tasks
#'@param keepPosCorOnly boolean indicating whether to only filter for positively correlated gene-peak associations (and hence assumes performing a right-tailed Z-test with permuted correlations). Default is TRUE (only keep positively correlated associations)
#'@param keepMultiMappingPeaks boolean indicating whether or not to keep peaks mapping to more than 1 gene. Default is FALSE (i.e. force 1-1 mapping for peak-gene, keeping peak with higher correlation if mapping to more than one gene)
#'@param n_bg numeric indicating the number of background correlations to compute per gene-peak pair. Default is 100 (increasing this can significantly increase run time)
#'@param p.cut numeric indicating p-value cut-off to apply to gene-peak correlations. Default is NULL (i.e. don't apply any filtering of results). Good option is to set to 0.05
#'@return a data.frame of correlations for each gene-peak overlap pair, with the observed correlation level and significance p-value based on background correlations per pair
#'@import Matrix dplyr pbmcapply SummarizedExperiment GenomicRanges
#'@export
#'@author Vinay Kartha
runGenePeakcorr <- function(ATAC.se, # SummarizedExperiment object of scATAC data
                            RNAmat, # Paired normalized scRNA-seq data, with gene names as rownames
                            genome, # Must be one of "hg19", "mm10", or "hg38"
                            geneList=NULL, # 2 or more valid gene symbols (if only running on subset of genes)
                            windowPadSize=50000, # base pairs padded on either side of gene TSS
                            normalizeATACmat=TRUE, # Whether or not to normalize scATAC counts (default is yes, assumes raw counts)
                            nCores=4, # Number of cores if parallelization support
                            keepPosCorOnly=TRUE,
                            keepMultiMappingPeaks=FALSE,
                            n_bg=100, # Number of background peaks to use
                            p.cut=NULL # Optional, if specified, will only return sig peak-gene hits
) {

  stopifnot(inherits(ATAC.se,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNAmat,c("Matrix","matrix")))

  if(!all.equal(ncol(ATAC.se),ncol(RNAmat)))
    stop("Input ATAC and RNA objects must have same number of cells")

  message("Assuming paired scATAC/scRNA-seq data ..")

  peakRanges.OG <- granges(ATAC.se) # Peak ranges in reference input SE (pre-filtering)

  # Function needs rownames for both matrices or gives error
  rownames(ATAC.se) <- paste0("Peak",1:nrow(ATAC.se))
  ATACmat <- assay(ATAC.se) # Rownames preserved

  # Normalize peak counts
  if(normalizeATACmat)
    ATACmat <- centerCounts(ATACmat) # Rownames preserved

  if(is.null(rownames(RNAmat)))
    stop("RNA matrix must have gene names as rownames")

  # Check for peaks/genes with 0 accessibility/expression

  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATAC.se))!=0
    ATAC.se <- ATAC.se[peaksToKeep,] # Subset ranges
    ATACmat <- ATACmat[peaksToKeep,]
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }


  peakRanges <- granges(ATAC.se) # Peak ranges

  if(any(Matrix::rowSums(RNAmat)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAmat)!=0
    RNAmat <- RNAmat[genesToKeep,]
  }

  cat("Number of peaks in ATAC data:",nrow(ATACmat),"\n")
  cat("Number of genes in RNA data:",nrow(RNAmat),"\n")


  if (!genome %in% c("hg19", "hg38", "mm10"))
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  switch(genome, hg19 = {
    TSSg <- FigR::hg19TSSRanges
  }, hg38 = {
    TSSg <- FigR::hg38TSSRanges
  }, mm10 = {
    TSSg <- FigR::mm10TSSRanges
  })

  # Keep genes that have annotation and are in RNA matrix
  names(TSSg) <- as.character(TSSg$gene_name)

  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")

    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }

    TSSg <- TSSg[geneList]
  }

  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg),rownames(RNAmat))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")

  # Match gene order in RNA matrix and TSS ranges
  RNAmat <- RNAmat[genesToKeep,]
  TSSg <- TSSg[genesToKeep]

  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg,
                                   width = windowPadSize,
                                   both = TRUE)

  # Get peak summit
  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges,width = 1,fix = "center")

  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlaps(query = TSSflank,subject = peakSummits)
  numPairs <- length(genePeakOv)

  cat("Found ",numPairs,"total gene-peak pairs for given TSS window ..\n")

  cat("Number of peak summits that overlap any gene TSS window: ",length(unique(subjectHits(genePeakOv))),"\n")
  cat("Number of gene TSS windows that overlap any peak summit: ",length(unique(queryHits(genePeakOv))),"\n\n")

  # For each gene, determine observed correlation of each overlapping peak to its associated gene (gene expression)

  # For each of those genes, also determine correlation based on background peaks (run in parallel) and save
  # Do this in parallel, and get data frame of gene-peak-pearson values
  # Fetch background peaks for each peak tested (i.e. that has overlap in window with gene)
  set.seed(123)
  cat("Determining background peaks ..\n")

  if(is.null(rowData(ATAC.se)$bias)){
    if(genome %in% "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if(genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if(genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    ATAC.se <- chromVAR::addGCBias(ATAC.se,genome=myGenome) }

  cat("Using ",n_bg," iterations ..\n\n")

  set.seed(123)
  bg <- chromVAR::getBackgroundPeaks(ATAC.se,niterations=n_bg)

  cat("Computing gene-peak correlations ..\n")

  pairsPerChunk <- 500

  # This defines the outer (larger chunks)
  largeChunkSize <- 5000

  startingPoint <- 1 # If for any reason the workers fail, resume from where it failed by specifying the starting point here
  chunkStarts <- seq(startingPoint, numPairs, largeChunkSize)
  chunkEnds <- chunkStarts + largeChunkSize -1
  chunkEnds[length(chunkEnds)] <- numPairs

  library(doParallel)

  dorcList <- list()
  for(i in 1:length(chunkStarts)){
    cat("Running pairs: ",chunkStarts[i], "to",chunkEnds[i],"\n")
    # This fill further chunk this up and run in parallel, saving the merged output ObsCor
    ObsCor <- FigR::PeakGeneCor(ATAC = ATACmat,
                                RNA = RNAmat,
                                OV = genePeakOv[chunkStarts[i]:chunkEnds[i]],
                                chunkSize = pairsPerChunk,
                                ncores = nCores,
                                bg = bg)
    gc()

    dorcList[[i]] <- ObsCor
  }

  cat("\nMerging results ..\n")
  dorcTab <- bind_rows(dorcList)

  cat("Performing Z-test for correlation significance ..\n")
  permCols <- 4:(ncol(bg)+3)


  if (keepPosCorOnly){
    # Filter to positive correlations
    cat("Only considering positive correlations ..\n")
    dorcTab <- dorcTab %>% dplyr::filter(rObs > 0)
  }

  if(!keepMultiMappingPeaks){
  # Remove multi-mapping peaks (force 1-1 mapping)
  cat("Keeping max correlation for multi-mapping peaks ..\n")
  dorcTab <- dorcTab %>% dplyr::group_by(Peak) %>% dplyr::filter(rObs==max(rObs))
  }

  # Swap gene number for gene symbol from TSS annotation lookup
  dorcTab$Gene <- as.character(TSSg$gene_name)[dorcTab$Gene]

  # Swap peak numbers to match reference input peak numbers
  # This only changes if some peaks had zero accessibility and were filtered out internally
  # Use rownames from reference matching
  dorcTab$Peak <- as.numeric(splitAndFetch(rownames(ATACmat)[dorcTab$Peak],"Peak",2))

  # # Z test pval
  dorcTab$rBgSD <- matrixStats::rowSds(as.matrix(dorcTab[,permCols]))
  dorcTab$rBgMean <- rowMeans(dorcTab[,permCols])
  dorcTab$pvalZ <- 1-stats::pnorm(q = dorcTab$rObs,mean = dorcTab$rBgMean,sd = dorcTab$rBgSD)


  cat("\nFinished!\n")

  if(!is.null(p.cut)){
    cat("Using significance cut-off of ",p.cut," to subset to resulting associations\n")
    dorcTab <- dorcTab[dorcTab$pvalZ <= p.cut,] # Subset to significant correlations only
  }

  # Add peak ranges to final data frame output
  dorcTab$PeakRanges <- paste(as.character(seqnames(peakRanges.OG[dorcTab$Peak])),paste(start(peakRanges.OG[dorcTab$Peak]),end(peakRanges.OG[dorcTab$Peak]),sep="-"),sep=":")

  return(as.data.frame(dorcTab[,c("Peak","PeakRanges","Gene","rObs","pvalZ")],stringsAsFactors=FALSE))
}

#' DORC scATAC-seq score summarization for single cells
#'
#'Function to compute single cell DORC scores for each gene identified as a DORC (see \code{\link[FigR]{runGenePeakcorr}} for determining DORCs through peak-gene correlations)
#'@param ATAC.se SummarizedExperiment object of the scATAC-seq reads in peak counts.
#'@param dorcTab data.frame object containing significant peak-gene pairs using which DORC scores will be computed. Must be a filtered set returned from \code{\link[FigR]{runGenePeakcorr}}. IMPORTANT: Make sure the exact same scATAC SE (peak set) was used when determining DORCs that is used here to get corresponding DORC peak counts
#'@param normalizeATACmat boolean indicating whether or not to normalize the counts present in the ATAC.se object prior to summing peak counts per gene. Default is TRUE (i.e. assumes peak counts are raw).
#'@param geneList character vector specifying a subset of genes to compute scores for. Useful if you already have a set list of DORCs and you only want to compute scores for those genes/peaks.
#'@param nCores numeric indicating the number of cores to use if parallelizing tasks
#'@return a Matrix object of scores pertaining to the sum of normalized scATAC peak counts per gene, as determined for significant gene-peak pairs when calling DORCs
#'@import Matrix dplyr pbmcapply SummarizedExperiment GenomicRanges
#'@export
#'@author Vinay Kartha
getDORCScores <- function(ATAC.se,
                          dorcTab,
                          normalizeATACmat=TRUE,
                          geneList=NULL,
                          nCores=4){


  if(!all(c("Peak","Gene") %in% colnames(dorcTab)))
    stop("The provided gene-peak table must have columns named Peak and Gene ..")

  if(any(dorcTab$Peak > nrow(ATAC.se)))
    stop("One or more peak indices in the gene-peak table are larger than the total number of peaks in the provided ATAC SE object ..\n Make sure the exact same SummarizedExperiment object is provided here as was for running the runGenePeakcorr function ..\n")

  if(!is.null(geneList)){
    if(!(all(geneList %in% as.character(dorcTab$Gene))))
      stop("One or more of the gene names supplied is not present in the gene-peak table provided..\n")

    if(length(geneList) > 50){
	    message("Running DORC scoring for ",length(geneList)," genes: ",paste(geneList[1:20],collapse=", "),", ... , ... , ... (truncated display)")
    } else {
      message("Running DORC scoring for ",length(geneList)," genes: ",paste(geneList,collapse = "\n"))
    }

    cat("........\n")

    dorcTab <- dorcTab[dorcTab$Gene %in% geneList,] # Filter only to these genes
    dorcGenes <- sort(as.character(unique(dorcTab$Gene)))
  } else {
    dorcGenes <- sort(as.character(unique(dorcTab$Gene)))
    cat("Running DORC scoring for all genes in annotation! (n = ",length(dorcGenes),")\n",sep="")
  }


  if(normalizeATACmat){
    # Normalize
    cat("Normalizing scATAC counts ..\n")
    ATAC.mat <- assay(centerCounts(ATAC.se,chunkSize = 5000))
    gc()
  } else {
    cat("Assuming provided scATAC counts are normalized ..\n")
    ATAC.mat <- assay(ATAC.se)
  }

  time_elapsed <- Sys.time()

  if(Sys.info()['sysname'] %in% "Windows"){
    message("Windows OS detected .. Cannot support parallilzation using mclapply for mc.cores > 1")
    message("Using 1 core instead ..\n")
    nCores <- 1
  }

  cat("Computing DORC scores ..\n")
  cat("Running in parallel using ", nCores, "cores ..\n")

  dorcMatL <- pbmcapply::pbmclapply(X=dorcGenes,
                                    FUN=function(x) {

                                      dorcPeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% x])

                                      if(length(dorcPeaks) > 1) {
                                        dorcCounts <- Matrix::colSums(ATAC.mat[dorcPeaks,])
                                      } else if(length(dorcPeaks==1)) {
                                        dorcCounts <- ATAC.mat[dorcPeaks,]
                                      }
                                    },mc.cores = nCores)

  dorcMat <- Matrix::Matrix(do.call('rbind',dorcMatL),sparse=TRUE)

  rownames(dorcMat) <- dorcGenes

  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)),"\n\n")

  return(dorcMat)

}


#' DORC J scatter plot
#'
#'Function to visualize genes based on number of significant gene-peak correlations (per gene), used to define DORCs
#'@param dorcTab data.frame object containing significant peak-gene pairs using which DORC scores will be computed. Must be a filtered set returned from \code{\link[FigR]{runGenePeakcorr}}. IMPORTANT: Make sure the exact same scATAC SE (peak set) was used when determining DORCs that is used here to get corresponding DORC peak counts
#'@param cutoff numeric indicating cut-off to use for number of significant peaks per gene to call it a DORC. Default is 7 significant peaks per gene.
#'@param labelTop numeric indicating the number of top genes to add labels for. Setting this to a very high number can lead to plot rendering difficulties (since labels get crowded).
#'@param returnGeneList boolean indicating whether to also return the DORCs passing the provided filter. Default is FALSE (will just plot and nothing else)
#'@param cleanLabels boolean indicating whether or not to try to have cleaner DORC gene labels. If TRUE, uses \code{\link[ggrepel]{geom_text_repel}}, if FALSE, will use regular \code{\link[ggplot2]{geom_text}}
#'@param labelSize numeric indicating what to set label size to. Default is 4.
#'@param ... additional parameters passed to either \code{\link[ggrepel]{geom_text_repel}} (if cleanLabels is set to TRUE) or \code{\link[ggplot2]{geom_text}}. For example, set size to control font size, or segment.length or segment.color to control call-out segment length and color etc.
#'@return a ggplot object of the DORC J plot. If returnGeneList is set to TRUE, then also returns the vector of DORC genes passing the specified threshold
#'@import dplyr ggplot2 ggrepel scales
#'@export
#'@author Vinay Kartha
dorcJPlot <- function(dorcTab,
                      cutoff=7, # Eventually we should use a knee caller for this
                      labelTop=25,
                      returnGeneList=FALSE,
                      cleanLabels=TRUE,
                      labelSize=4,
                      ...){

  stopifnot(all(c("Peak","Gene","pvalZ") %in% colnames(dorcTab)))

  # Count the number of significant peak associations for each gene (without pre-filtering genes)
  numDorcs <- dorcTab  %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
  numDorcs$Index <- 1:nrow(numDorcs) # Add order index
  numDorcs %>% as.data.frame(stringsAsFactors=FALSE) -> numDorcs
  rownames(numDorcs) <- numDorcs$Gene

  dorcGenes <- numDorcs$Gene[numDorcs$n >= cutoff]

  numDorcs <- numDorcs %>%
    mutate(isDORC=ifelse(Gene %in% dorcGenes,"Yes","No")) %>%
    mutate(Label=ifelse(Gene %in% dorcGenes[1:labelTop],Gene,""))

  # Plot
  dorcG <- ggplot(numDorcs,aes(x=Index,y=n,color=isDORC,label=Label)) +
    geom_hline(linetype="dotted",yintercept = cutoff)+
    geom_vline(linetype="dotted",xintercept = max(numDorcs[numDorcs$Gene %in% dorcGenes,"Index"]))+
    geom_point(size=0.8) +
    geom_line()+
    scale_color_manual(values=c("gray65","firebrick"))+
    scale_y_continuous(breaks = scales::pretty_breaks())+
    theme_classic() +
    labs(y="Number of correlated peaks",x="Ranked genes",title=paste0("# DORCs: ( n >= ",cutoff,") = ",length(dorcGenes)))+
    theme(axis.text = element_text(color = "black"),legend.position = "none",plot.title=element_text(hjust=0.5)) +
    scale_x_reverse() # flip so we can add labels later, if needed, with more space

  if(cleanLabels){
    dorcG <- dorcG + ggrepel::geom_label_repel(size=labelSize,max.iter = 100,max.overlaps = Inf,fontface="italic",...)
  } else {
    dorcG <- dorcG + ggplot2::geom_text(size=labelSize,fontface="italic",...)
  }

  print(dorcG)

  if(returnGeneList)
    return(dorcGenes)

}
