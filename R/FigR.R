# Script housing main FigR wrapper functions

### Author: Vinay Kartha
### Contact: <vinay_kartha@g.harvard.edu>
### Affiliation: Buenrostro Lab, Department of Stem Cell and Regenerative Biology, Harvard University

#' Infer FigR TF-DORC associations
#'
#'Function to run TF motif-to-gene associations using reference DORC peak-gene mappings and TF RNA expression levels
#'@param ATAC.se SummarizedExperiment object of peak x cell scATAC-seq data, the same as used to compute DORCs using \code{\link[FigR]{runGenePeakcorr}}
#'@param dorcK numeric specifying the number of dorc nearest-neighbors to pool peaks from for the motif enrichment per DORC. Default is 30, i.e. set to ~3 percent of total DORCs determined
#'@param dorcTab data.frame object containing significant peak-gene pairs using which DORC scores will be computed. Must be a filtered set returned from \code{\link[FigR]{runGenePeakcorr}}. IMPORTANT: Make sure the exact same scATAC SE peak set was used when determining DORCs that is used here to get corresponding DORC peak counts
#'@param n_bg number of background peaks to use for
#'@param genome character specifying a valid genome assembly to use for peak GC content estimation and background peak determination. Must be one of "hg19","hg38", or "mm10", and requires the corresponding genomes package e.g. \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}} for hg19
#'@param dorcMat Matrix object of smoothed single-cell DORC accessibility scores
#'@param rnaMat Matrix object of smoothed single-cell RNA expression values
#'@param dorcGenes character vector specifying the subset of DORCs to test, if not running on everything. Note: We still use the entire list of DORCs found in dorcMat to determine dorc KNNs from, but will only test and include results for these specified genes (also must exist in the provided RNA matrix rnaMat as rownames)
#'@param nCores numeric specifying the number of cores to run DORCs in parallel. Default is 1, i.e. don't use parallel backend
#'@return a data.frame with all TF-DORC motif enrichment and correlation associations, and the corresponding FigR regulation score for each association
#'@import dplyr Matrix SummarizedExperiment chromVAR
#'@export
#'@author Vinay Kartha
#'
runFigRGRN <- function(ATAC.se, # SE of scATAC peak counts. Needed for chromVAR bg peaks etc.
                       dorcK=30, # How many dorc kNNs are we using to pool peaks
                       dorcTab, # peak x DORC connections (should contain indices relative to peaks in ATAC.se)
                       n_bg=50, # No. of background peaks to use for motif enrichment Z test
                       genome, # One of mm10, hg19, hg38, with no default
                       dorcMat, # Expect smoothed
                       rnaMat, # Expect smoothed
                       dorcGenes=NULL, # If only running on a subset of genes
                       nCores=1
){
  # Must be matched data
  stopifnot(all.equal(ncol(dorcMat),ncol(rnaMat)))

  # Expects "Gene" / "Peak" in dorcTab
  if(!all(c("Peak","Gene") %in% colnames(dorcTab)))
    stop("Expecting fields Peak and Gene in dorcTab data.frame .. see runGenePeakcorr function in BuenRTools")

  if(all(grepl("chr",dorcTab$Peak,ignore.case = TRUE))) {
    usePeakNames <- TRUE
    message("Detected peak region names in Peak field")

    if(!(all(grepl("chr",rownames(ATAC.se),ignore.case = TRUE))))
      stop("Peak regions provided in dorcTab data.frame but not found as rownames in input SE")

    if(!all(dorcTab$Peak %in% rownames(ATAC.se)))
      stop("Found DORC peak region not present in input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  } else{
    usePeakNames <- FALSE
    message("Assuming peak indices in Peak field")
  # If using index, make sure no indices are outside range of SE
    if(max(dorcTab$Peak) > nrow(ATAC.se))
      stop("Found DORC peak index outside range of input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  }


  if(is.null(dorcGenes)) {
    dorcGenes <- rownames(dorcMat)
  } else {
    cat("Using specified list of dorc genes ..\n")
    if (!(all(dorcGenes %in% rownames(dorcMat)))) {
      cat("One or more of the gene names supplied is not present in the DORC matrix provided: \n")
      cat(dorcGenes[!dorcGenes %in% rownames(dorcMat)], sep = ", ")
      cat("\n")
      stop()
    }
  }

  DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))),k = dorcK)$nn.index # Scaled
  rownames(DORC.knn) <- rownames(dorcMat)

  if (is.null(SummarizedExperiment::rowData(ATAC.se)$bias)) {
    if (genome %in% "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if (genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
  }

  # Set data subfolder path
  packagePath <- find.package("FigR", lib.loc=NULL, quiet = TRUE)

  if(grepl("hg",genome)){
    pwm <- readRDS(paste0(packagePath,"/data/cisBP_human_pfms_2021.rds"))
  } else {
    pwm <- readRDS(paste0(packagePath,"/data/cisBP_mouse_pfms_2021.rds"))
  }

  # Old motif naming convention
  if(all(grepl("_",names(pwm),fixed = TRUE)))
     names(pwm) <- FigR::extractTFNames(names(pwm))

  message("Removing genes with 0 expression across cells ..\n")
  rnaMat <- rnaMat[Matrix::rowSums(rnaMat)!=0,]
  myGeneNames <- gsub(x = rownames(rnaMat),pattern = "-",replacement = "") # NKX2-1 to NKX21 (e.g.)
  rownames(rnaMat) <- myGeneNames

  # Only non-zero expression TFs (also found in rnaMat)
  motifsToKeep <- intersect(names(pwm),myGeneNames)

  # This has to be done on the full SE (same peakset used as input to dorc calling)
  cat("Getting peak x motif matches ..\n")
  motif_ix <- motifmatchr::matchMotifs(subject = ATAC.se,pwms = pwm[motifsToKeep],genome=genome)

  # Keep TFs with some peak x motif match
  motif_ix <- motif_ix[,Matrix::colSums(assay(motif_ix))!=0]

  cat("Determining background peaks ..\n")
  cat("Using ", n_bg, " iterations ..\n\n")
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    ATAC.mat <- assay(ATAC.se)
    ATAC.mat <- cbind(ATAC.mat,1)
    ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=ATAC.mat),rowRanges = granges(ATAC.se))
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
  } else {
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  }

  # For each DORC, do motif enrichment among dorc sig Peaks, and correlation of DORC accessibility (smoothed) to TF RNA levels

  cat("Testing ",length(motifsToKeep)," TFs\n")
  cat("Testing ",nrow(dorcMat)," DORCs\n")
  library(doParallel)
  if(nCores > 1)
    message("Running FigR using ",nCores," cores ..\n")
  opts <- list()
  pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(nCores)
  clusterEvalQ(cl, .libPaths())
  doSNOW::registerDoSNOW(cl)
  mZtest.list <- foreach(g=dorcGenes,
                         .options.snow = opts,
                         .packages = c("FigR", "dplyr","Matrix","Rmpfr")) %dopar%   {
                           # Take peaks associated with gene and its k neighbors
                           # Pool and use union for motif enrichment
                           DORCNNpeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% c(g,rownames(dorcMat)[DORC.knn[g,]])])

                           if(usePeakNames)
                             DORCNNpeaks <- which(rownames(ATAC.se) %in% DORCNNpeaks) # Convert to index relative to input

                           mZ <- FigR::motifPeakZtest(peakSet = DORCNNpeaks,
                                                bgPeaks = bg,
                                                tfMat = assay(motif_ix))

                           mZ <- mZ[,c("gene","z_test")]
                           colnames(mZ)[1] <- "Motif"
                           colnames(mZ)[2] <- "Enrichment.Z"
                           mZ$Enrichment.P <- 2*pnorm(abs(mZ$Enrichment.Z),lower.tail = FALSE) # One-tailed
                           mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
                           mZ <- cbind("DORC"=g,mZ)
                           # Correlate smoothed dorc with smoothed expression, with spearman
                           corr.r <- cor(dorcMat[g,],t(as.matrix(rnaMat[mZ$Motif,])),method = "spearman")
                           stopifnot(all.equal(colnames(corr.r),mZ$Motif))

                           mZ$Corr <- corr.r[1,] # Correlation coefficient
                           mZ$Corr.Z <- scale(mZ$Corr,center = TRUE,scale = TRUE)[,1] # Z-score among all TF correlations
                           mZ$Corr.P <- 2*pnorm(abs(mZ$Corr.Z),lower.tail = FALSE) # One-tailed
                           mZ$Corr.log10P <- sign(mZ$Corr.Z)*-log10(mZ$Corr.P)
                           return(mZ)
                         }
  cat("Finished!\n")
  cat("Merging results ..\n")
  # Merge and save table for downstream filtering and plotting (network)
  TFenrich.d <- do.call('rbind',mZtest.list)
  dim(TFenrich.d)
  rownames(TFenrich.d) <- NULL

  # Make combined score based on multiplication
  # Here, we only sign by corr
  # Since sometimes we lose digit precision (1 - v small number is 1, instead of 0.9999999..)
  # Use Rmpfr, increase precision limits above default (100 here)
  TFenrich.d <- TFenrich.d %>% dplyr::mutate("Score"=sign(Corr)*as.numeric(-log10(1-(1-Rmpfr::mpfr(Enrichment.P,100))*(1-Rmpfr::mpfr(Corr.P,100)))))
  TFenrich.d$Score[TFenrich.d$Enrichment.Z < 0] <- 0
  TFenrich.d
}


#' Rank TF drivers
#'
#' Ranked plot of TF activators and repressors based on their inferred regulation score
#'@param figR.d data.frame of results returned by \code{\link[FigR]{runFigRGRN}})
#'@param rankBy character specifying one of "meanScore" or "nTargets" to either rank TFs by the mean regulation score across all genes, or by the total number of inferred activated or repressed targets passing a specified (absolute) regulation score, respectively
#'@param myLabels character vector specifying the subset of TFs to highlight on the plot, if rankBy is set to "meanScore". Useful if you want to see where your TFs of interest lie. If NULL (Default), we label the top and bottom 95 percentile TFs
#'@param score.cut numeric specifying the absolute regulation score to threshold TF-DORC connections on, only if rankBy is set to "nTargets". Default is 1 if "nTargets" and no custom cut-off is specified
#'@param interactive boolean indicating whether or not to allow interactive hover-over utility for more label information (useful if visualizing too many TFs and labels are hard to distinguish). Default is FALSE
#'@return a ggplot2 object of the resulting plot
#'@import dplyr ggplot2 ggrepel plotly
#'@export
#'@author Vinay Kartha
rankDrivers <- function(figR.d,
                        rankBy=c("meanScore","nTargets"),
                        myLabels=NULL,
                        score.cut=NULL,
                        interactive=FALSE){

  if(!rankBy %in% c("meanScore","nTargets"))
    stop("rankBy parameter has to be one of meanScore or nTargets to rank drivers using ..\n")

  if(rankBy %in% "meanScore"){
    message("Ranking TFs by mean regulation score across all DORCs ..\n")

    # Prep summary stats
    figR.summ <- figR.d %>%  group_by(Motif) %>%
      dplyr::summarise(Score=mean(Score)) %>%
      arrange(desc(Score)) %>% # Activating to Repressing
      mutate(Motif=factor(Motif,levels=as.character(Motif)))

    # Top and bottom %ile labels
    figR.summ$TF <- as.character(figR.summ$Motif)

    if(is.null(myLabels)){
      # Use quantiles to define what labels are drawn
      figR.summ$TF[figR.summ$Score >= quantile(figR.summ$Score,0.05) & figR.summ$Score <= quantile(figR.summ$Score,0.95)] <- ""
    } else {
      # Only highlight user-specified
      figR.summ$TF[!figR.summ$TF %in% myLabels] <- ""
    }

    library(ggrepel)

    gAll <- ggplot(figR.summ,aes(x=Motif,y=Score,label=TF)) +
      geom_bar(size=0.1,stat="identity",fill="darkorange",color=NA) +
      theme_classic() + theme(axis.text.x = element_blank(),axis.text=element_text(color="black")) +
      ggrepel::geom_text_repel(size=3,min.segment.length = 0.1,segment.size = 0.2,max.overlaps = 20) +
      geom_hline(yintercept = 0) + labs(x="TF Motifs",y="Regulation Score")

  } else {
    message("Ranking TFs by total number of associated DORCs ..\n")

    if(is.null(score.cut)){
      message("Regulation score cut-off not specified ..\n")
      score.cut <- 1
    }

    message("Using absolute score cut-off of: ",score.cut," ..\n")

    # Prep summary stats
    figR.summ  <- figR.d %>% dplyr::filter(abs(Score) >= score.cut) %>%
      group_by(Motif) %>% dplyr::select(-DORC) %>%
      dplyr::summarize(numActivated = sum(Score > 0), numRepressed = sum(Score < 0)) %>%
      dplyr::mutate(diff = numActivated -numRepressed) %>% # log2FC works too
      mutate(numActivatedY=ifelse(diff > 0,numActivated,-numActivated),numRepressedY=ifelse(diff > 0,numRepressed,-numRepressed)) %>%
      dplyr::arrange(desc(diff)) %>%
      mutate(Motif=factor(Motif,levels=as.character(Motif))) %>%
      dplyr::select(-diff) %>% reshape2::melt(id.vars=c("Motif","numActivated","numRepressed"))

    # New
    # Make ggplot
    gAll <- figR.summ %>%
      ggplot(aes(x=Motif,y=value,fill=variable,numActivated=numActivated,numRepressed=numRepressed)) +
      geom_bar(stat="identity",color="lightgray",size=0.1) + theme_classic() + geom_hline(yintercept = 0) +
      scale_fill_manual(values=c("firebrick3","steelblue4"),labels=c("Activated","Repressed")) +
      theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=6),axis.text = element_text(color="black")) +
      labs(x = "Ranked TF Motifs", y = paste0("# Associated genes \nabs(Score) >= ", score.cut),fill="Class") +
      scale_y_continuous(labels=abs)
}

 if(!interactive) {
  gAll
 } else {
   if(rankBy %in% "meanScore"){
     plotly::ggplotly(gAll)
   }    else{
   plotly::ggplotly(gAll + theme(legend.position = "none",
                                 axis.text.x = element_blank()),
                    tooltip=c("Motif","numActivated","numRepressed"))
  }

  }
}



#' Plot FigR scatter profile
#'
#' Scatter plot visualization of filtered TF-DORC associations based on the enrichment of each motif among the queried DORC's peaks and the correlation of the TF RNA to the DORC accessibility score
#'@param figR.d data.frame of results returned by \code{\link[FigR]{runFigRGRN}})
#'@param marker character specifying a valid DORC gene to restrict TF drivers for
#'@param score.cut numeric specifying the absolute regulation score to threshold TF-DORC connections on. Default is 1
#'@param label boolean indicating whether or not to add text labels for TF drivers passing score filter
#'@return a ggplot2 object of the scatter plot of predicted TF drivers for the specified DORC
#'@import dplyr ggplot2 scales
#'@export
#'@author Vinay Kartha
plotDrivers <- function(figR.d,
                        marker,
                        score.cut=1,
                        label=TRUE){

  if(!marker %in% figR.d$DORC)
    stop("Marker specified is not a valid DORC symbol found in the data.frame")

  d <- figR.d %>% dplyr::filter(DORC %in% marker) %>% mutate(isSig=ifelse(abs(Score) >= score.cut,"Yes","No"))
  if(label){
  d$Label <- d$Motif
  d$Label[d$isSig %in% "No"] <- ""
  } else {
    d$Label <- ""
  }

  gScatter <- d %>%  ggplot(aes(x=Corr.log10P,y=Enrichment.log10P,color=isSig,label=Label)) +
    geom_hline(yintercept = 0,color="gray60",linetype="dashed") +
    geom_vline(xintercept = 0,color="gray60",linetype="dashed") +
    geom_point(size=0.8) + theme_classic() +
    scale_color_manual(values=c("gray66","firebrick3"))+
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    scale_y_continuous(breaks=scales::pretty_breaks()) +
    labs(y="Enrichment log10 P",x="Correlation log10 P",title=marker) +
    ylim(-ceiling(max(abs(d$Enrichment.log10P))), ceiling(max(abs(d$Enrichment.log10P))))+
    xlim(-ceiling(max(abs(d$Corr.log10P))), ceiling(max(abs(d$Corr.log10P))))+
    theme(legend.position = "none",axis.text = element_text(color="black"),
          plot.title = element_text(hjust=0.5,face="italic"),panel.background = element_rect(fill=NA)) +
    geom_text(hjust=1.1,fontface="italic",color="black",size=3)

  gScatter

}

#' Plot FigR heatmap
#'
#' Heatmap visualization of TF-DORC associations based on the regulation scores inferred by FigR
#'@param figR.d data.frame of results returned by \code{\link[FigR]{runFigRGRN}}).
#'@param score.cut numeric specifying the absolute regulation score to threshold TF-DORC connections on. Default is 1
#'@param DORCs character specifying valid DORC gene symbols to subset heatmap to. Default is NULL (no subsetting)
#'@param TFs character specifying valid TF gene symbols to subset heatmap to. Default is NULL (no subsetting)
#'@param ... additional parameters passed to the \code{\link[ComplexHeatmap]{Heatmap}})
#'@return a TF-DORC filtered Heatmap generatd using \code{\link[ComplexHeatmap]{Heatmap}})
#'@import dplyr ComplexHeatmap BuenColors scales reshape2 tibble circlize
#'@export
#'@author Vinay Kartha
plotfigRHeatmap <- function(figR.d,
                            score.cut=1,
                            DORCs=NULL,
                            TFs=NULL,
                            ... # Additional params passed to ComplexHeatmap
){


  message("Using absolute score cut-off of: ",score.cut," ..\n")

  DORCsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut) %>% pull(DORC) %>% unique()
  TFsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut) %>% pull(Motif) %>% unique()


  if(!is.null(DORCs)){
    if(!all(DORCs %in% figR.d$DORC))
      stop("One or more DORCs specified is not a valid DORC symbol found in the data.frame")
    DORCsToKeep <- intersect(DORCsToKeep,DORCs)
    TFsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut & DORC %in% DORCsToKeep) %>% pull(Motif) %>% unique()
  }


  if(!is.null(TFs)){
    if(!all(TFs %in% figR.d$Motif))
      stop("One or more TFs specified is not a valid TF symbol found in the data.frame")
    TFsToKeep <- intersect(TFsToKeep,TFs)
    DORCsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut & Motif %in% TFsToKeep) %>% pull(DORC) %>% unique()
  }


  net.d <- figR.d %>% dplyr::filter(DORC %in% DORCsToKeep & Motif %in% TFsToKeep) %>%
    reshape2::dcast(DORC ~ Motif) %>%
    tibble::column_to_rownames("DORC") %>% as.matrix()

  message("Plotting ",nrow(net.d)," DORCs x ",ncol(net.d), "TFs\n")

  # Heatmap view

  myCols <- circlize::colorRamp2(seq(-2,2,length.out = 9),colors = BuenColors::jdb_palette("solar_flare"))
  myHeat <- ComplexHeatmap::Heatmap(net.d,
                                    col=myCols,
                                    clustering_distance_rows = "pearson",
                                    clustering_distance_columns = "pearson",
                                    name="Score",border = TRUE,
                                    row_names_gp = gpar(fontsize=5,fontface="italic"),...)

  myHeat

}

#' Plot FigR Network
#'
#' Network visualization of TF-DORC associations based on the regulation scores inferred by FigR
#'@param figR.d data.frame of results returned by \code{\link[FigR]{runFigRGRN}})
#'@param score.cut numeric specifying the absolute regulation score to threshold TF-DORC connections on. Default is 1
#'@param DORCs character specifying valid DORC gene symbols to subset heatmap to. Default is NULL (no subsetting)
#'@param TFs character specifying valid TF gene symbols to subset heatmap to. Default is NULL (no subsetting)
#'@param weight.edges boolean specifying whether or not to weight edges by FigR regulation score. Default is FALSE
#'@param TFnodecol character specifying valid color name to use for TF nodes. Default is Tomato
#'@param DORCnodecol character specifying valid color name to use for DORC nodes. Default is Sky Blue
#'@param posEdgecol character specifying valid color name to use for Activating edges between TFs and DORCs. Default is Forest Green
#'@param negEdgecol character specifying valid color name to use for Repressive edges between TFs and DORCs. Default is Purple
#'@param labelSize numeric specifying font size to use for all labels. Default is 13
#'@param showLegend boolean indicating whether to show color legend. Default is TRUE
#'@return a network plot of the resulting filtered TF-DORC associations
#'@import dplyr networkD3 BuenColors scales reshape2 tibble circlize
#'@export
#'@author Vinay Kartha
plotfigRNetwork <- function(figR.d,
                        score.cut=1,
                        DORCs=NULL,
                        TFs=NULL,
                        weight.edges=FALSE,
                        TFnodecol='Tomato',
                        DORCnodecol='Sky Blue',
                        posEdgecol='Forest Green',
                        negEdgecol='Purple',
                        labelSize=13,
                        showLegend=TRUE){
# Network view

# Filter
net.dat <-  figR.d %>% dplyr::filter(abs(Score) >= score.cut)

if(!is.null(DORCs))
  net.dat <- net.dat %>% dplyr::filter(DORC %in% DORCs)

if(!is.null(TFs))
  net.dat <- net.dat %>% dplyr::filter(Motif %in% TFs)

net.dat$Motif <- paste0(net.dat$Motif, ".")
net.dat$DORC <- paste0(net.dat$DORC)


dorcs <- data.frame(name = unique(net.dat$DORC), group = "DORC", size = 8)
tfs <- data.frame(name = unique(net.dat$Motif), group = "TF", size = 3)
nodes <- rbind(dorcs,tfs)

edges <- as.data.frame(net.dat)

# Make edges into links (subtract 1 for 0 indexing)
links <- data.frame(source=unlist(lapply(edges$Motif, function(x) {which(nodes$name==x)-1})),
                    target=unlist(lapply(edges$DORC, function(x) {which(nodes$name==x)-1})),
                    corr=edges$Corr,
                    enrichment=edges$Enrichment.P)

links$Value <- scales::rescale(edges$Score)*20

# Set of colors you can choose from for TF/DORC nodes
#colors <- c("Red", "Orange", "Yellow", "Green", "Blue", "Purple", "Tomato", "Forest Green", "Sky Blue","Gray","Steelblue3","Firebrick2","Brown")
colors=c(TFnodecol,DORCnodecol,posEdgecol,negEdgecol)
nodeColorMap <- data.frame(color = colors, hex = gplots::col2hex(colors))

getColors <- function(tfColor, dorcColor = NULL) {
  temp <- c(as.character(nodeColorMap[nodeColorMap$color==tfColor,]$hex),
            as.character(nodeColorMap[nodeColorMap$color==dorcColor,]$hex))
  if (is.null(dorcColor)) {
    temp <- temp[1]
  }
  colors <- paste(temp, collapse = '", "')
  colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')
  colorJS
}

networkD3::forceNetwork(Links = links,
             Nodes = nodes,
             Source = "target",
             Target = "source",
             NodeID = "name",
             #NodeID="myLabel",
             Group = "group",
             Value = "Value",
             Nodesize = "size",
             #linkWidth = 1.5, # Fixed link weight
             radiusCalculation = "Math.sqrt(d.nodesize)*2",
             arrows = FALSE,
             opacityNoHover = 0.6,
             opacity = 1,
             zoom = TRUE,
             bounded = TRUE,
             charge = -15,
             fontSize = labelSize,
             legend = showLegend,
             colourScale = getColors(tfColor = TFnodecol,dorcColor =  DORCnodecol), # TF then DORC
             linkColour = ifelse(links$corr > 0, as.character(nodeColorMap[nodeColorMap$color==posEdgecol,]$hex),
                                 as.character(nodeColorMap[nodeColorMap$color==negEdgecol,]$hex)))

}
