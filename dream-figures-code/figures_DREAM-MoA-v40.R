


####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ########      ########     ########    #########    ####     ##  ##########                  #########        ##      ##########       ##                                    
##     ##     ##     ##     ##     ##          ##           ##           ## ##    ##      ##                      ##      ##      ####         ##          ####                              
##     ##     ##     ##    ##      ##          ##           ##           ##  ##   ##      ##                      ##       ##    ##  ##        ##         ##  ##                             
##     ########      #######       ########     #######     ########     ##   ##  ##      ##                      ##       ##   ##    ##       ##        ##    ##                                
##     ##            ##    ##      ##                 ##    ##           ##    ## ##      ##                      ##       ##  ##########      ##       ##########                                                                                             
##     ##            ##     ##     ##                 ##    ##           ##     ####      ##                      ##      ##   ##      ##      ##       ##      ##                           
##     ##            ##      ##    ########     #######     #########    ##      ###      ##                      #########    ##      ##      ##       ##      ##                                      
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#DATE: 4-1-2020
#START MAKING FIGURES FOR PAPER:
#  (1)  INTRO FIGURES TO DATA-SET
#  (2)  ANALYSIS of results


#setwd("/Volumes/ifs/home/c2b2/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA")
setwd("/Volumes/ac_lab/efd2115/02 Research/02 Coding-Data/2019-10-10 - DREAM-MoA/")



load("DREAM_kinaseInhibitor_FINAL-DATA.RData")
rm(FINAL_DREAM_Dose_Responses)
rm(FINAL_DREAM_RNAseq)

load("HALLMARKS_drugAVERAGE.RData")
rm(PANGEA_ges_HALLMARKS_aREA)
rm(PANGEA_ges_HALLMARKS_targetAVG)
rm(PANGEA_vpr_HALLMARKS_drugAVG)
rm(PANGEA_vpr_HALLMARKS_target)
PANACEA_ges_drugAvg_Hallmarks<-PANGEA_ges_HALLMARKS_drugAVG$nes
rm(PANGEA_ges_HALLMARKS_drugAVG)


colnames(PANACEA_ges_drugAvg_Hallmarks)<-gsub("-","",gsub(" ","",gsub("_","",toupper(colnames(PANACEA_ges_drugAvg_Hallmarks)))))

OVERLAP_ges_drugAvg_Hallmarks<-PANACEA_ges_drugAvg_Hallmarks[,rownames(FINAL_DREAM_kinome_Kds)]
##################################################################################################################################
#PLOTTING FUNCTIONS:
##################################################################################################################################

library("gplots")
library("RColorBrewer")
library("devtools")

heatmap.3 <- function(x,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                      breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                      trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}



hcluster_corr_n1screen_DrugTargets_hallmarks<-function(matrix,database,target_vector,target_colors,main_title,num_var=30){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)
  master_indx_filter<-c()
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #2 extract class-indices in rlab
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(rownames(matrix),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 make master index to filter matrix
    master_indx_filter<-c(master_indx_filter,rlab_class_idx)
  }
  master_indx_filter_unique<-unique(master_indx_filter)
  
  matrix<-matrix[master_indx_filter_unique,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix_topvar,method="spearman")
  cor_rows <- cor(t(matrix_topvar),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(20,12),
              RowSideColors=rlab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(20,12),
              symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }
  
  return(rownames(matrix_topvar[hclust_rows$order,]))
}



hcluster_euclid_n1screen_DrugTargets_kinasesDREAM<-function(matrix,database,target_vector,target_colors,main_title,num_var=30){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)
  master_indx_filter<-c()
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #2 extract class-indices in rlab
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(rownames(matrix),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 make master index to filter matrix
    master_indx_filter<-c(master_indx_filter,rlab_class_idx)
  }
  master_indx_filter_unique<-unique(master_indx_filter)
  
  matrix<-matrix[master_indx_filter_unique,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  
  dist_columns <- dist(t(matrix_topvar))
  dist_rows <- dist(matrix_topvar)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(12,12),
              RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=as.dendrogram(hclust_rows),
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="both", margins=c(12,12),
              symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }
  
  return(rownames(matrix_topvar[hclust_rows$order,]))
}



hcluster_corr_n1screen_DrugTargets_hallmarks_kinomeCLUSTER<-function(matrix,database,target_vector,target_colors,main_title,num_var=30,
                                                                     rownames_order){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)

  
  matrix<-matrix[kinome_row_order,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix_topvar,method="spearman")
  cor_rows <- cor(t(matrix_topvar),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=NA,
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="column", margins=c(20,12),
              RowSideColors=rlab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=NA,
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="column", margins=c(20,12),
              symbreaks=TRUE, key=TRUE, symkey=TRUE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1,col=colorRampPalette(c("blue","blue", "white", "red", "red"))(n = 20),
              KeyValueName="NES",keysize = 0.8,xlab="Hallmarks",ylab="Drug-Averages (across cell-lines)")
  }
  
  
}

hcluster_euclid_n1screen_DrugTargets_kinasesDREAM_panaceaCLUSTER<-function(matrix,database,target_vector,target_colors,main_title,num_var=30,
                                                                           rownames_order){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix_idx<-which(rownames(matrix) %in% rownames(database))
  
  matrix<-matrix[matrix_idx,]
  
  
  
  #FILTER MATRIX BY ONLY DRUGS ANNOTATED: (reduce noise of drugs that do nothing...)

  
  matrix<-matrix[rownames_order,]
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(target_vector),ncol=length(rownames(matrix)))
  rownames(rlab)<-target_vector
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(target_vector)){
    target<-target_vector[row]
    
    #1 extract drugs-vector from n1screen_drugDB
    targets_list<-strsplit(database$target,",")
    db_class_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    class_names_tmp<-rownames(database[db_class_idx,])
    
    
    rlab_class_idx<-which(sapply(strsplit(colnames(rlab),"\\."),function(x) x[[1]][1]) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-target_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  top_var<-names(sort(apply(matrix,2,var),decreasing=T))[1:num_var]
  matrix_topvar<-matrix[,top_var]
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  
  dist_columns <- dist(t(matrix_topvar))
  dist_rows <- dist(matrix_topvar)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(length(target_vector)>1){
    heatmap.3(matrix_topvar, 
              Rowv=NA,
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="column", margins=c(12,12),
              RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              RowSideColorsSize=5, KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }else{
    heatmap.3(matrix_topvar, 
              Rowv=NA,
              Colv=as.dendrogram(hclust_columns), 
              na.rm = TRUE, scale="none", dendrogram="column", margins=c(12,12),
              symbreaks=FALSE, key=TRUE, symkey=FALSE,
              density.info="none", trace="none", main=main_title, cexRow=1,cexCol=1.2,col=colorRampPalette(c("blue","white", "red"))(n = 20),
              KeyValueName="-log10(Kd)",keysize = 0.8,xlab="Kinases",ylab=paste0(as.character(nrow(matrix_topvar))," Drugs"))
  }
  
  return(rownames(matrix_topvar[hclust_rows$order,]))
}



###################################################################################################################################################
#1 INTRODUCTION FIGURES
###################################################################################################################################################

#MIX OF INHIBITORS:  some dissagreement with Drug-Banks
kinome_row_order<-hcluster_euclid_n1screen_DrugTargets_kinasesDREAM(matrix=FINAL_DREAM_kinome_Kds,
                                                  database=FINAL_DREAM_drugLibrary_annotation,
                                                  target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                                  target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                                  main_title = "KINOME-BINDING DATABASE")



panacea_row_order<-hcluster_corr_n1screen_DrugTargets_hallmarks(matrix=t(OVERLAP_ges_drugAvg_Hallmarks),
                                                  database=FINAL_DREAM_drugLibrary_annotation,
                                                  target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                                  target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                                  main_title = "Differential-mRNA Database")


#############################
#RE-ORDER CLUSTERS TO MATCH OTHER-DATA-SET DRUG-order
##############################
hcluster_euclid_n1screen_DrugTargets_kinasesDREAM_panaceaCLUSTER(matrix=FINAL_DREAM_kinome_Kds,
                                                                    database=FINAL_DREAM_drugLibrary_annotation,
                                                                    target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                                                    target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                                                    main_title = "KINOME-BINDING DATABASE",
                                                                    rownames_order=panacea_row_order[length(panacea_row_order):1])

hcluster_corr_n1screen_DrugTargets_hallmarks_kinomeCLUSTER(matrix=t(OVERLAP_ges_drugAvg_Hallmarks),
                                             database=FINAL_DREAM_drugLibrary_annotation,
                                             target_vector = c("EGFR","FLT3","ABL1","RET","MET","AKT2"),
                                             target_colors = c("darkred","darkorange","steelblue","blue","darkgreen","red"),
                                             main_title = "Differential-mRNA Database",
                                             rownames_order=kinome_row_order)




###################################################################################################################################################
#2 PLOTTING LEADER-BOARD SCORES:
###################################################################################################################################################
LB_SC1<-read.csv("raw_scores/LB_sc1.csv",header=F)$"V1"
LB_SC2<-read.csv("raw_scores/LB_sc2.csv",header=F)$"V1"


hist(LB_SC2[1:25],breaks=20,col="yellowgreen",xlim=c(0,60),ylim=c(0,35),cex.main=1.5,
     main="Subchallenge 2 (SC2) Leaderboard Results",xlab="-log10(p-value)",cex=1.5,cex.axis=1.5,cex.lab=1.5)
hist(LB_SC2[26:53],breaks=6,col="gold2",add=T)
hist(LB_SC2[54:86],breaks=1,col="firebrick",add=T)





hist(LB_SC1[1:14],breaks=6,col="yellowgreen",xlim=c(0,17),ylim=c(0,10),cex.main=1.5,
     main="Subchallenge 1 (SC1) Leaderboard Results",xlab="-log10(p-value)",cex=1.5,cex.axis=1.5,cex.lab=1.5)
hist(LB_SC1[14:39],breaks=20,col="gold2",add=T)
hist(LB_SC1[40:86],breaks=3,col="firebrick",add=T)



































































####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########   ##         ##     ##     ########   ##########   #########   ########     ########    ###     ##      #########                                                      
##     ##          ##         ##     ##    ##              ##       ##          ##     ##       ##       ####    ##     ## 
##    ##           ##         ##     ##    ##              ##       ##          ##     ##       ##       ## ##   ##    ##                 
##    ##           ##         ##     ##     #######        ##       #######     ########        ##       ##  ##  ##    ##    ####
##    ##           ##         ##     ##           ##       ##       ##          ##    ##        ##       ##   ## ##    ##       ##
##     ##          ##         ##     ##           ##       ##       ##          ##     ##       ##       ##    ####     ##      ##       
##      ########   ########    #######     ########        ##       #########   ##      ##   ########    ##     ###      ########     
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################


load("finalScores_drugs.RData")

load("finalScores_drugs_FISHER.RData")


library("gplots", lib.loc="~/Library/R/3.3/library")
library("RColorBrewer", lib.loc="~/Library/R/3.3/library")


hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  dist_columns <- dist(t(matrix))
  dist_rows <- dist(matrix)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "NES",key.xlab = "NES")
}

hcluster_corr<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(10,12), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "NES",key.xlab = "NES")
}


#GSEA- all predictions

hcluster_corr(t(finalScores_drugs_ceiling[,-12]),title="Drug-wise Scores per Team",col_size = 0.8,row_size=1.5)

hcluster_corr(t(apply(finalScores_drugs_ceiling[,-12],2,rank)),title="Drug-wise Scores per Team",col_size = 1.3,row_size=1.5)


#FISHER - top 10 predictions

hcluster_corr(t(finalScores_drugs_top10_ceiling[,-12]),title="Drug-wise Scores per Team",col_size = 0.8,row_size=1.5)

hcluster_corr(t(apply(finalScores_drugs_top10_ceiling[,-12],2,rank)),title="Drug-wise Scores per Team",col_size = 1.3,row_size=1.5)



####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ########     ########    ##     ##     ########              ########       ##       ###     ##   ##     ##                                       
##    ##     ##    ##     ##   ##     ##    ##                     ##     ##     ####      ####    ##   ##    ##                     
##    ##      ##   ##     ##   ##     ##   ##                      ##     ##    ##  ##     ## ##   ##   ##   ##                         
##    ##      ##   ########    ##     ##   ##    ####    ######    ########    ##    ##    ##  ##  ##   ######                    
##    ##      ##   ##    ##    ##     ##   ##       ##             ##     ##  ##########   ##   ## ##   ##   ##                      
##    ##     ##    ##     ##   ##     ##    ##      ##             ##     ##  ##      ##   ##    ####   ##    ##          
##    ########     ##      ##   #######      ########              ########   ##      ##   ##     ###   ##     ##                                                 
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##      ########     #########    ##########        ##              
##     ##           ##            ##               ####                                        
##    ##            ##            ##              ##  ##                              
##    ##    ####     ########     ########       ##    ##                                                                                                                                  
##    ##       ##           ##    ##            ##########                                                         
##     ##      ##           ##    ##           ##        ##                                                 
##      ########     ########     ##########  ##          ##                                                                                
####################################################################################################################################################################################################################################################################
##    ########     ########    ##     ##     ########             ##############     ##        ########      ########    ########   ##########                  ##       ##      ##     #########                                                                  
##    ##     ##    ##     ##   ##     ##    ##                          ##          ####       ##     ##    ##           ##             ##                     ####      ##      ##    ##                                              
##    ##      ##   ##     ##   ##     ##   ##                           ##         ##  ##      ##     ##   ##            ##             ##                    ##  ##     ##      ##   ##                                               
##    ##      ##   ########    ##     ##   ##    ####    ######         ##        ##    ##     ########    ##     ####   ########       ##      ######       ##    ##     ##    ##    ##     #####                                 
##    ##      ##   ##    ##    ##     ##   ##       ##                  ##       ##########    ##    ##    ##        ##  ##             ##                  ##########     ##  ##     ##         ##                                      
##    ##     ##    ##     ##   ##     ##    ##      ##                  ##      ##        ##   ##     ##    ##      ##   ##             ##                 ##        ##     ####       ##       ##                                        
##    ########     ##      ##   #######      ########                   ##     ##          ##  ##      ##    ########    ########       ##                ##          ##     ##         #########                                                                         
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################



load("n1druglibrary_annotation_CTEP.RData")
n1druglibrary_annotation_CAPS<-n1druglibrary_annotation_CTEP
rownames(n1druglibrary_annotation_CAPS)<-make.names(gsub("_","",gsub("-","",gsub(" ","",toupper(rownames(n1druglibrary_annotation_CAPS))))),unique=T)

load("finalScores_drugs.RData")
load("finalScores_drugs_FISHER.RData")



########################
#SC1 & SC2 tightly correlated
########################
plot(rowMeans(apply(finalScores_drugs_ceiling[,],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),],2,rank)),col="white")
text(rowMeans(apply(finalScores_drugs_ceiling,2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),],2,rank)),labels=rownames(finalScores_drugs_ceiling))

########################
#Netphar vs Atom
########################

plot(rowMeans(apply(finalScores_drugs_ceiling[,c("netphar","netphar")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("Atom","Atom")],2,rank)),col="white")
text(rowMeans(apply(finalScores_drugs_ceiling[,c("netphar","netphar")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("Atom","Atom")],2,rank)),labels=rownames(finalScores_drugs_ceiling))

plot(rowMeans(apply(finalScores_drugs_ceiling[,c("netphar","netphar")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("SBNB","SBNB")],2,rank)),col="white")
text(rowMeans(apply(finalScores_drugs_ceiling[,c("netphar","netphar")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("SBNB","SBNB")],2,rank)),labels=rownames(finalScores_drugs_ceiling))

plot(rowMeans(apply(finalScores_drugs_ceiling[,c("Atom","Atom")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("SBNB","SBNB")],2,rank)),col="white")
text(rowMeans(apply(finalScores_drugs_ceiling[,c("Atom","Atom")],2,rank)),rowMeans(apply(finalScores_drugs_top10_ceiling[rownames(finalScores_drugs_ceiling),c("SBNB","SBNB")],2,rank)),labels=rownames(finalScores_drugs_ceiling))


library("viper")
library("gplots")
library("RColorBrewer", lib.loc="~/Library/R/3.3/library")
library("devtools")
heatmap.3_sidebarLabel <- function(x,drug_label_size=1,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                                   breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                                   trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                                   labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                                   densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE,cex.axis=drug_label_size)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

final_targets<-unique(unlist(strsplit(n1druglibrary_annotation_CAPS$target,",")))


DREAM_aREA_drugTargets<-function(matrix,database,DB_targets,top_drug_num=25,stat_sig=c("all_sig","inhibit_sig","activate_sig","top_num")){

  trim25_pangea_matrix<-matrix
  
  #######################################
  #RANKED DRUGS & LINES PLOTTING:
  #######################################
  RANKED_DRUGS_all<-names(sort(rowMeans(trim25_pangea_matrix,na.rm=T),decreasing=FALSE))
  RANKED_LINES_all<-names(sort(colMeans(trim25_pangea_matrix[RANKED_DRUGS_all[1:20],],na.rm=T),decreasing=FALSE))
  
  if (stat_sig=="activate_sig"){
    RANKED_LINES_all<-names(sort(colMeans(trim25_pangea_matrix[names(sort(rowMeans(trim25_pangea_matrix,na.rm=T),decreasing=T))[1:20],],na.rm=T),decreasing=T))
    
  }
  
  matrix<-trim25_pangea_matrix[RANKED_DRUGS_all,RANKED_LINES_all]
  
  
  
  avg_drug_rank<-sort(rowMeans(trim25_pangea_matrix,na.rm=T),decreasing=T)
  
  
  
  
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #0 RENAME MATRIX to names compatible with database 
  rownames(database)<-toupper(rownames(database))
  rownames(matrix)<-toupper(rownames(matrix))
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(DB_targets),ncol=length(rownames(matrix)))
  rownames(rlab)<-DB_targets
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  #target="MAP2K1"
  for (target in DB_targets){
    
    #1 extract drug-set given target
    targets_list<-strsplit(database$target,",")
    target_drug_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    target_drug_names<-rownames(database[target_drug_idx,])
    rlab_class_idx<-which(colnames(rlab) %in% target_drug_names)
    
    #3 replace "white" w/ color assigned to each class:
    rlab[target,rlab_class_idx]<-"black"
  }
  
  
  
  
  #######################################
  #RANKED DRUG-CLASSES & LINES PLOTTING:
  #######################################
  #extract average NES-vector to run aREA on:
  avg_nes_gene_vector<-as.matrix(sort(rowMeans(matrix,na.rm=T),decreasing=FALSE))
  
  drug_class_regulon<-list()
  
  #MAKE drug-class regulon:
  #target="MAP2K1"
  for (target in DB_targets){
    
    
    targets_list<-strsplit(database$target,",")
    target_drug_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
    target_drug_names<-rownames(database[target_drug_idx,])
    target_drug_names_filtered<-target_drug_names[which(target_drug_names %in% rownames(matrix))]
    num_drugs<-length(target_drug_names_filtered)
    
    #3 define regulon:
    drug_class_regulon[[target]]<-list()
    drug_class_regulon[[target]]$tfmode<-rep(1,num_drugs)
    drug_class_regulon[[target]]$likelihood<-rep(1,num_drugs)
    names(drug_class_regulon[[target]]$tfmode)<-target_drug_names_filtered
    
  }
  
  
  ################
  #USE: STAT-SIG or RANKING
  ################
  #RUN aREA to get ranking:
  aREA_object<-aREA(avg_nes_gene_vector,drug_class_regulon,minsize = 1)
  ranked_drug_classes_values<-sort(aREA_object$nes[,1],decreasing=F)
  ranked_drug_classes<-names(ranked_drug_classes_values)
  
  #PLOT ALL STATISTICALLY SIGNIFICANT DRUGS:
  if (stat_sig=="all_sig"){
    rlab_sorted<-rlab[ranked_drug_classes[which(abs(ranked_drug_classes_values)>1.96)],]
  }
  #PLOT TOP # OF DRUGS
  if(stat_sig=="top_num"){
    rlab_sorted<-rlab[ranked_drug_classes[1:top_drug_num],]
  }
  
  #PLOT ONLY INHIBITING STAT-SIG DRUGS:
  if (stat_sig=="inhibit_sig"){
    rlab_sorted<-rlab[ranked_drug_classes[which(ranked_drug_classes_values<(-1.96))],]
  }
  
  #PLOT ONLY ACTIVATING STAT-SIG DRUGS:
  if (stat_sig=="activate_sig"){
    rlab_sorted<-rlab[ranked_drug_classes[which(ranked_drug_classes_values>-10)],]
  }
  
  
  label_size=min(max(6/nrow(rlab_sorted),1),2.75)
  
  
  
  
  
  
  
  #########################################
  #COLUMN TOP-BAR: cell-lines:
  #########################################
  
  
  
  
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  if(stat_sig=="activate_sig"){
    heatmap.3_sidebarLabel(matrix[(nrow(matrix):1),], drug_label_size=label_size,
                           Rowv=NA,
                           Colv=NA, 
                           na.rm = TRUE, scale="none", dendrogram="none", margins=c(13,8),
                           RowSideColors=rlab_sorted[nrow(rlab_sorted):1,ncol(rlab_sorted):1], symbreaks=FALSE, key=TRUE, symkey=FALSE,
                           density.info="none", trace="none", main=NA, 
                           cexRow=1,cexCol=1.5,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
                           RowSideColorsSize=45, KeyValueName="NES",keysize = 0.8,xlab=NA,ylab=NA)
  }else{
    heatmap.3_sidebarLabel(matrix, drug_label_size=label_size,
                           Rowv=NA,
                           Colv=NA, 
                           na.rm = TRUE, scale="none", dendrogram="none", margins=c(10,6),
                           RowSideColors=rlab_sorted, symbreaks=TRUE, key=TRUE, symkey=TRUE,
                           density.info="none", trace="none", main=NA, 
                           cexRow=0.2,cexCol=1.5,col=colorRampPalette(c("blue", "blue","white", "red", "red"))(n = 20),
                           RowSideColorsSize=15, KeyValueName="NES",keysize = 0.8,xlab="cell-lines",ylab="Drugs")
  }
  
  
  
  return(avg_drug_rank)
}

####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################

##################################
#GSEA metric = SC2
##################################
#final targets include representative Dr

#RAW SCORE:
nes_avg<-DREAM_aREA_drugTargets(matrix=finalScores_drugs_ceiling[,c("Atom","netphar","SBNB")],
                       database=n1druglibrary_annotation_CTEP,
                       DB_targets = final_targets,
                       stat_sig = "activate_sig")

#RANK-Transformed Score:
rank_avg<-DREAM_aREA_drugTargets(matrix=apply(finalScores_drugs_ceiling[,c("Atom","netphar","SBNB")],2,rank),
                                 database=n1druglibrary_annotation_CTEP,
                                 DB_targets = final_targets,
                                 stat_sig = "activate_sig")


n1druglibrary_annotation_CAPS[names(rank_avg),"DrugBank_targets"]


rank_transform<-apply(finalScores_drugs_ceiling[,-12],2,rank)

write.csv(rank_transform,file="rank-tansformed_SC2.csv")
##################################
#FISHER metric = SC1
##################################


#RAW SCORE:
nes_avg<-DREAM_aREA_drugTargets(matrix=finalScores_drugs_top10_ceiling[,c("Atom","netphar","SBNB")],
                                database=n1druglibrary_annotation_CTEP,
                                DB_targets = final_targets,
                                stat_sig = "activate_sig")

#RANK-Transformed Score:
rank_avg<-DREAM_aREA_drugTargets(matrix=apply(finalScores_drugs_top10_ceiling[,c("Atom","netphar","SBNB")],2,rank),
                                 database=n1druglibrary_annotation_CTEP,
                                 DB_targets = final_targets,
                                 stat_sig = "activate_sig")


n1druglibrary_annotation_CAPS[names(rank_avg),"DrugBank_targets"]



####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##      ########    #######     ########      #######    ##              ##      ########### ########    ########   ####    ##                                                                            
##     ##          ##     ##    ##     ##     ##    ##   ##             ####         ##         ##      ##      ##  ## ##   ##                
##    ##           ##     ##    ##     ##     ##    ##   ##            ##  ##        ##         ##      ##      ##  ##  ##  ##                         
##    ##           ##     ##    ########      #######    ##           ##    ##       ##         ##      ##      ##  ##   ## ##                 
##    ##           ##     ##    ##    ##      ##   ##    ##          ##########      ##         ##      ##      ##  ##    ####                    
##     ##          ##     ##    ##     ##     ##    ##   ##         ##        ##     ##         ##      ##      ##  ##     ###                  
##      ########    #######     ##      ##    ##     ##  ########  ##          ##    ##      ########    ########   ##      ##                             
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#UNDERSTAND CORRELATION STRUCTURE: kinome-binding vs MR-perturbations (across overlapping 88 drugs)




load("simple_PANGEA-HALLMARKS_KinomeBinding_overlap.RData")

CORRmatrix_kinomeBIND_effectorACTIVITY<-cor(kinase_inhibitors_Binding,kinase_inhibitors_avgHALLMARKS,method="spearman")

plot(CORRmatrix_kinomeBIND_effectorACTIVITY["MAP2K2",],CORRmatrix_kinomeBIND_effectorACTIVITY["BRAF",])



##################################################################
#PLOTTING FUNCTIONS:
##################################################################
library("gplots")
library("RColorBrewer")
library("devtools")

hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  dist_columns <- dist(t(matrix))
  dist_rows <- dist(matrix)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "-log10(Kd)",key.xlab = "-log10(Kd)")
}
hcluster_corr<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(15,15), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "NES",key.xlab = "NES")
}

na_idx<-which(is.na(rowMeans(CORRmatrix_kinomeBIND_effectorACTIVITY))==TRUE)
CORRmatrix_kinomeBIND_effectorACTIVITY_na<-CORRmatrix_kinomeBIND_effectorACTIVITY[-na_idx,]


hcluster_corr(CORRmatrix_kinomeBIND_effectorACTIVITY_na,col_size = 0.5,row_size=0.4)

ranked_genes<-names(sort(apply(CORRmatrix_kinomeBIND_effectorACTIVITY_na,2,var),decreasing=TRUE))
ranked_kinases<-names(sort(apply(CORRmatrix_kinomeBIND_effectorACTIVITY_na,1,var),decreasing=TRUE))


hcluster_corr(t(CORRmatrix_kinomeBIND_effectorACTIVITY_na[ranked_kinases[1:200],ranked_genes[1:40]]),col_size = 0.6,row_size=1,title="Kinome-N1-Correlation:  top-2000-var genes")



##################################################################
#GENOME & KINOME ANNOTATION:
##################################################################
library("gplots")
library("RColorBrewer")
library("devtools")
heatmap.3 <- function(x,Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,distfun = dist,hclustfun = hclust,dendrogram = c("both","row", "column", "none"),symm = FALSE,scale = c("none","row", "column"),na.rm = TRUE,revC = identical(Colv,"Rowv"),add.expr,
                      breaks,symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",col = "heat.colors",colsep,rowsep,sepcolor = "white",sepwidth = c(0.05, 0.05),cellnote,notecex = 1,notecol = "cyan",na.color = par("bg"),
                      trace = c("none", "column","row", "both"),tracecol = "cyan",hline = median(breaks),vline = median(breaks),linecol = tracecol,margins = c(5,5),ColSideColors,RowSideColors,side.height.fraction=0.3,cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,labCol = NULL,key = TRUE,keysize = 1.5,density.info = c("none", "histogram", "density"),denscol = tracecol,symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,main = NULL,xlab = NULL,ylab = NULL,lmat = NULL,lhei = NULL,lwid = NULL,ColSideColorsSize = 1,RowSideColorsSize = 1,KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

load("KINOME-classifications.RData")

#Symetric Key:  NES's
hcluster_corr_MULTI_sidebar_G2M_mod2<-function(matrix,database,DB_col,DB_class,class_names,class_colors,main_title,G2M_plus=TRUE,x_lab="Genes"){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  #0 RENAME MATRIX to names compatible with database 
  rownames(database)<-make.names(toupper(rownames(database)),unique = TRUE)
  
  rownames(matrix)<-make.names(toupper(rownames(matrix)),unique = TRUE)
  
  
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix<-matrix[rownames(database),]
  
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(class_names),ncol=length(rownames(matrix)))
  rownames(rlab)<-class_names
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(class_names)){
    #1 extract 1st rows DB-entry
    db_col_tmp<-DB_col[row]
    db_class_tmp<-DB_class[row]
    
    #2 extract class-indices in rlab
    db_class_idx<-which(database[,db_col_tmp] == db_class_tmp)
    class_names_tmp<-rownames(database[db_class_idx,])
    rlab_class_idx<-which(colnames(rlab) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-class_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  
  
  
  #########################################
  #XX G2M TOP-bar
  #########################################
  
  ###1 extract relevant G2M genes:
  G2M_chckpt<-c("TOP2A","CENPF","CHEK1","FOXM1","CENPI","BLM","BRIP1","MCM6","MCM4","E2F1","E2F2","E2F3","E2F8","BUB1B","MYBL2","MCM7","FEN1","AURKA","AURKB","CENPU","FANCA","TYMS","RAD51","HMGA1","DNMT3B","E2F7","BRCA1","BRCA2","PRKDC","CENPK","MCM3","DMNT1","FANCD2","RAN","RANBP1","MCM2","CENPZ","RACGAP1","PCNA","H2AFX","MCM5")
  G2M_chckpt_plus<-unique(c("TOP2A","CENPF","CHEK1","FOXM1","CENPI","BLM","BRIP1","MCM6","MCM4","E2F1","E2F2","E2F3","E2F8","BUB1B","MYBL2","MCM7","FEN1","AURKA","AURKB","CENPU","FANCA","TYMS","RAD51","HMGA1","DNMT3B","E2F7","BRCA1","BRCA2","PRKDC","CENPK","MCM3","DMNT1","FANCD2","RAN","RANBP1","MCM2","CENPZ","RACGAP1","PCNA","H2AFX","MCM5","TRIP13","ASF1B","UHRF1","ATAD2","ECT2","TTK","WHSC1","ARHGAP11A","ZNF367","TCF19","HELLS","RACGAP1","MKI67","SPAG5","UBE2C","TONSL","RUVBLT","ZWINT","HMGB2","DEPDC1","CBX3","RFC4","CHAF1A","GTSE1","HMGN1","PSMC3IP","HDGF","CKS1B","TCEB1","CKS2","CCNA2","ENY2","GGCT","ILE2","HSPE1","PSRC1","NDC80","TMPO","CDK2","GMPS","EIF4EBP1","PA2G4","GMNN","ZC3H15","PTTG1","ZC3H1","CYCS","GTF2E2","HLTF","TCF19","WHSC1","TTF2","SHOX2","EIF2AK2","ARNTL2","UHRF1","GTSE1","DEPDC1","SUV39H2","ZWINT","YEASTS2","LRPPRC","NAA15","POLD1","RFC4","RBL1","MCM3","ILF3","CDK7","CIT","GLRA2","PHF19","CDCA7","ZNF695","EZH2","PCNA","HNRNPAB","ZCRB1","SUB1","FH","HSPE1","RAB1A","CCDC47","COPS3","H2AFX","DAXX","HDGF","APH1A","PUF60","POLR3K","HSBP1","IDH2","VPS72","PRMT1","UBE2L3","PPP2CA","HSPA5","MCM2","HDAC2","CEBPZ","TIMELESS","UBEC2","PTTG1","TRAIP","MNX1","ARHGEF39","CHEK2","IQGAP3","HSPD1"))
  if (G2M_plus == TRUE){
    G2M_chckpt_filtered<-G2M_chckpt_plus[which(G2M_chckpt_plus %in% colnames(matrix))]
  }
  else {
    G2M_chckpt_filtered<-G2M_chckpt[which(G2M_chckpt %in% colnames(matrix))]
  }
  
  ###2 MAKE column color-key:
  clab<-matrix("white",nrow=length(colnames(matrix)),ncol=1)
  
  matrix_G2M_idx<-which(colnames(matrix) %in% G2M_chckpt_filtered)
  
  clab[matrix_G2M_idx,1]<-"darkblue"
  
  colnames(clab)=c("G2M")
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  heatmap.3(matrix, 
            Rowv=as.dendrogram(hclust_rows),
            Colv=as.dendrogram(hclust_columns), 
            na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,6),
            ColSideColors=clab, RowSideColors=rlab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
            density.info="none", trace="none", main=main_title, cexRow=0.6,cexCol=0.1,col=colorRampPalette(c("blue",  "white",  "red"))(n = 20),
            ColSideColorsSize=1, RowSideColorsSize=5, KeyValueName="spearman",keysize = 0.8,xlab=x_lab,ylab="kinases")
  
}

hcluster_corr_MULTI_sidebar_G2M_mod2_transpose<-function(matrix,database,DB_col,DB_class,class_names,class_colors,main_title,G2M_plus=TRUE,x_lab="Genes"){
  #########################################
  #1 MAKE DATABASE & DRUG-SCREEN COMPATIBLE:
  #########################################
  
  #0 RENAME MATRIX to names compatible with database 
  rownames(database)<-make.names(toupper(rownames(database)),unique = TRUE)
  
  rownames(matrix)<-make.names(toupper(rownames(matrix)),unique = TRUE)
  
  
  
  
  #1 slice database by drugs in matrix (and then vice versa)
  DB_drug_idx<-which(rownames(database) %in% rownames(matrix))
  database<-database[DB_drug_idx,]
  
  matrix<-matrix[rownames(database),]
  
  
  #########################################
  #1 MAKE Key-matrix for each pair of DB_col and DB_class
  #########################################
  #ROW-KEY MATRIC:     rows = drug_class          cols = drug
  
  ###1 make blank ("white") drug_class matrix
  rlab<-matrix("white",nrow=length(class_names),ncol=length(rownames(matrix)))
  rownames(rlab)<-class_names
  colnames(rlab)<-rownames(matrix)
  
  #fill in each row with labels for each class listed
  for (row in 1:length(class_names)){
    #1 extract 1st rows DB-entry
    db_col_tmp<-DB_col[row]
    db_class_tmp<-DB_class[row]
    
    #2 extract class-indices in rlab
    db_class_idx<-which(database[,db_col_tmp] == db_class_tmp)
    class_names_tmp<-rownames(database[db_class_idx,])
    rlab_class_idx<-which(colnames(rlab) %in% class_names_tmp)
    
    #3 replace "white" w/ color assigned to each class:
    class_color_temp<-class_colors[row]
    rlab[row,rlab_class_idx]<-class_color_temp
  }
  
  
  clab<-t(rlab)
  matrix<-t(matrix)
  
  
  #########################################
  #XX DEFINE distance/clustering functions:
  #########################################
  #mydist=function(c) {as.dist(1-cor(c,method="spearman"))}
  #myclust=function(c) {hclust(c,method="average")}
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  #########################################
  #XX PLOT HEATMAP
  #########################################
  heatmap.3(matrix, 
            Rowv=as.dendrogram(hclust_rows),
            Colv=as.dendrogram(hclust_columns), 
            na.rm = TRUE, scale="none", dendrogram="both", margins=c(3,20),
            ColSideColors=clab, symbreaks=TRUE, key=TRUE, symkey=TRUE,
            density.info="none", trace="none", main=main_title, cexRow=0.9,cexCol=0.3,col=colorRampPalette(c("blue",  "white",  "red"))(n = 20),
            ColSideColorsSize=5, KeyValueName="spearman",keysize = 0.8,xlab="kinases",)
  
}




hcluster_corr_MULTI_sidebar_G2M_mod2(CORRmatrix_kinomeBIND_effectorACTIVITY_na,kinome_classification_hugo_filtered,
                                     DB_col = c(3,3,3,3,3,3,3),
                                     DB_class = c("AGC","CAMK","CK1","CMGC","STE","TK","TKL"),
                                     class_names =c("AGC","CAMK","CK1","CMGC","STE","TK","TKL"),
                                     class_colors = c("darkred","red","darkorange","darkgreen","cyan4","darkblue","blueviolet"),
                                     main_title = "KINOME-binding/N1-screen Correlation Matrix:")



hcluster_corr_MULTI_sidebar_G2M_mod2_transpose(CORRmatrix_kinomeBIND_effectorACTIVITY_na[,ranked_genes[1:40]],kinome_classification_hugo_filtered,
                                     DB_col = c(3,3,3,3,3,3,3),
                                     DB_class = c("AGC","CAMK","CK1","CMGC","STE","TK","TKL"),
                                     class_names =c("AGC","CAMK","CK1","CMGC","STE","TK","TKL"),
                                     class_colors = c("darkred","red","darkorange","darkgreen","cyan4","darkblue","blueviolet"),
                                     main_title = "KINOME-binding/N1-screen Correlation Matrix:")





####################################################################################################################################
#TRANSFORM KINASES INTO KEGG PATHWAYS
####################################################################################################################################

library(viper)
load("kegg_pathways_regulon.RData")



KinasePathway_mRNAHallmark_corr<-aREA(CORRmatrix_kinomeBIND_effectorACTIVITY_na,kegg_pathways_regulon,minsize = 2)


hcluster_corr(KinasePathway_mRNAHallmark_corr$nes,col_size = 0.6,row_size=0.5,title="Kinome-N1-Correlation:  top-2000-var genes")


sort(rownames(KinasePathway_mRNAHallmark_corr$nes))

KinasePathway_mRNAHallmark_corr$nes



                        
kegg_relevant<-c("APOPTOSIS","CALCIUM_SIGNALING_PATHWAY","CELL_CYCLE","CHEMOKINE_SIGNALING_PATHWAY","CYTOSOLIC_DNA_SENSING_PATHWAY",
                 "ERBB_SIGNALING_PATHWAY","INSULIN_SIGNALING_PATHWAY","JAK_STAT_SIGNALING_PATHWAY","MAPK_SIGNALING_PATHWAY",
                 "MTOR_SIGNALING_PATHWAY","NUCLEOTIDE_EXCISION_REPAIR","P53_SIGNALING_PATHWAY","PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM",
                 "REGULATION_OF_AUTOPHAGY","TGF_BETA_SIGNALING_PATHWAY","TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","VEGF_SIGNALING_PATHWAY","WNT_SIGNALING_PATHWAY")

hallmark_relevant<-c("APOPTOSIS","DNA_REPAIR","E2F_TARGETS","FATTY_ACID_METABOLISM","G2M_CHECKPOINT","GLYCOLYSIS","HEDGEHOG_SIGNALING",
  "HYPOXIA","IL2_STAT5_SIGNALING","IL6_JAK_STAT3_SIGNALING","INFLAMMATORY_RESPONSE",
  "INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","KRAS_SIGNALING_DN","KRAS_SIGNALING_UP","MITOTIC_SPINDLE",
  "MTORC1_SIGNALING","MYC_TARGETS_V1","MYC_TARGETS_V2","NOTCH_SIGNALING",
  "OXIDATIVE_PHOSPHORYLATION","P53_PATHWAY","PANCREAS_BETA_CELLS","PEROXISOME","PI3K_AKT_MTOR_SIGNALING",
  "PROTEIN_SECRETION","REACTIVE_OXIGEN_SPECIES_PATHWAY","TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB",
  "UNFOLDED_PROTEIN_RESPONSE","WNT_BETA_CATENIN_SIGNALING")

hcluster_corr(KinasePathway_mRNAHallmark_corr$nes[kegg_relevant,hallmark_relevant],col_size = 1,row_size=1,title="KEGG-Pathways vs mRNA-Programs")

sort(KinasePathway_mRNAHallmark_corr$nes[kegg_relevant,"MYC_TARGETS_V1"])
sort(KinasePathway_mRNAHallmark_corr$nes[kegg_relevant,"APOPTOSIS"])


kegg_relevant2<-c("APOPTOSIS","CALCIUM_SIGNALING_PATHWAY","CELL_CYCLE",
                 "ERBB_SIGNALING_PATHWAY","INSULIN_SIGNALING_PATHWAY","JAK_STAT_SIGNALING_PATHWAY","MAPK_SIGNALING_PATHWAY",
                 "MTOR_SIGNALING_PATHWAY","P53_SIGNALING_PATHWAY","PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM",
                 "TGF_BETA_SIGNALING_PATHWAY","VEGF_SIGNALING_PATHWAY")

kegg_relevant3<-c("CALCIUM_SIGNALING_PATHWAY","CELL_CYCLE",
                  "ERBB_SIGNALING_PATHWAY","INSULIN_SIGNALING_PATHWAY","JAK_STAT_SIGNALING_PATHWAY",
                  "P53_SIGNALING_PATHWAY",
                  "TGF_BETA_SIGNALING_PATHWAY","VEGF_SIGNALING_PATHWAY")

hcluster(KinasePathway_mRNAHallmark_corr$nes[kegg_relevant3,c("MYC_TARGETS_V1","MYC_TARGETS_V2","E2F_TARGETS","G2M_CHECKPOINT","TGF_BETA_SIGNALING","APOPTOSIS")],col_size = 1,row_size=1,title="KEGG-Pathways vs mRNA-Programs")













####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##    ########     ########    ##     ##     ########             ########       ##      ###     ##   ##     ## 
##    ##     ##    ##     ##   ##     ##    ##                    ##     ##     ####     ####    ##   ##    ##               
##    ##      ##   ##     ##   ##     ##   ##                     ##     ##    ##  ##    ## ##   ##   ##   ##                 
##    ##      ##   ########    ##     ##   ##    ####    ######   ########    ##    ##   ##  ##  ##   ######         
##    ##      ##   ##    ##    ##     ##   ##       ##            ##     ##  ##########  ##   ## ##   ##   ##                                  
##    ##     ##    ##     ##   ##     ##    ##      ##            ##     ##  ##      ##  ##    ####   ##    ##
##    ########     ##      ##   #######      ########             ########   ##      ##  ##     ###   ##     ##                                    
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########    #######      ####    ##    ########  ##     ##   #######    ##    #######    ###     ##                                                
##     ##          ##     ##     ## ##   ##    ##        ##     ##  ##          ##   ##     ##   ####    ##                                                  
##    ##           ##     ##     ##  ##  ##    ##        ##     ##  ##          ##   ##     ##   ## ##   ##                                                           
##    ##           ##     ##     ##   ## ##    #######   ##     ##   #######    ##   ##     ##   ##  ##  ##                                                            
##    ##           ##     ##     ##    ####    ##        ##     ##         ##   ##   ##     ##   ##   ## ##                                                     
##     ##          ##     ##     ##     ###    ##        ##     ##         ##   ##   ##     ##   ##    ####                                                   
##      ########    #######      ##      ##    ##         #######    #######    ##    #######    ##     ###                                                                         
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#NOT ENOUGH TEAMS TO DO DEEP DIVE IN METHODS OR DATA BASES AS:
#         -> data sources quite diverse (each person had a unique data-set and so hard to deconvolue NORM & Data biases)
#GOAL:  #1 analyze drugbank vs kinome gold-standards
#       #2 see if non-drug-bank targets picked up in submissions...




#GOLD STANDARD:
load("OVERLAP_n1-Kinome-Screens.RData")
rm(kinome_selectivity_scores_OVERLAP)


load("final_results_processed.RData")


####################################################################################################################################
#1 RECOVERY OF DRUG-BANK TARGETS
#####################################################################################################################################
kinome_GoldStandard<-kinome_inhibitor_matrix_pKa_OVERLAP[rownames(final_results_processed$SBNB),which(colnames(kinome_inhibitor_matrix_pKa_OVERLAP) %in% colnames(final_results_processed$SBNB))]

max(kinome_GoldStandard)

threshold_matrix<-matrix(ncol=length(seq(4.6,8.5,0.01)),nrow=length(rownames(final_results_processed$SBNB)))
rownames(threshold_matrix)<-rownames(final_results_processed$SBNB)
colnames(threshold_matrix)<-seq(4.6,8.5,0.01)

#drug="SORAFENIB"
for (drug in rownames(threshold_matrix)){
  cutoff="8.5"
  for (cutoff in colnames(threshold_matrix)){
    #ID DrugBank Kinases (for this drug):
    targets_list<-strsplit(n1druglibrary_annotation_CTEP_OVERLAP[drug,"target"],",")[[1]]
    KINASES_drugbank<-targets_list[which(targets_list %in% colnames(kinome_GoldStandard))]
    
    #ID KINOME Drugs (for this kinase)
    KINASES_kinome<-names(which(kinome_GoldStandard[drug,]>as.numeric(cutoff)))
    
    #FRACTION DRUG-BANK KINASES Recovered
    drugbank_recovered<-length(intersect(KINASES_kinome,KINASES_drugbank))/length(KINASES_drugbank)
    
    #STORE:
    threshold_matrix[drug,cutoff]<-drugbank_recovered
  }
}

plot(as.numeric(colnames(threshold_matrix)),colMeans(threshold_matrix),type="n",main="Recovered DrugBank Targets Based on Kd Threshold",
     ylab="Fraction Drug-Bank Kinase Targets",xlab="Affinity Threshold Used (-log10(Kd))",cex.main=1.5,cex.lab=1.5,cex.axis=1.5,xlim=c(8.5,4.5))
lines(as.numeric(colnames(threshold_matrix)),colMeans(threshold_matrix),lwd=3)
abline(v=6,col="red",lty=2)

####################################################################################################################################
#  Kinases-Drug <uM events NOT in DrugBank:
#####################################################################################################################################
kinome_inhibitor_matrix_pKa_OVERLAP_filtered<-kinome_inhibitor_matrix_pKa_OVERLAP[,which(colnames(kinome_inhibitor_matrix_pKa_OVERLAP) %in% colnames(final_results_processed$SBNB))]

drug_unique_num<-c()
drug_overlap_num<-c()

Bootstrapped_targets<-matrix(nrow=nrow(final_results_processed$SBNB),ncol=6)
rownames(Bootstrapped_targets)<-rownames(final_results_processed$SBNB)
colnames(Bootstrapped_targets)<-c("Atom_overlap","Atom_kinome","netphar_overlap","netphar_kinome","SBNB_overlap","SBNB_kinome")

#drug="DASATINIB"
for (drug in rownames(final_results_processed$SBNB)){
  #################################
  #1 Characterize Gold-Standard Overlap:
  #################################
  
  #ID DrugBank Kinases (for this drug):
  targets_list<-strsplit(n1druglibrary_annotation_CTEP_OVERLAP[drug,"target"],",")[[1]]
  KINASES_drugbank<-targets_list[which(targets_list %in% colnames(kinome_inhibitor_matrix_pKa_OVERLAP_filtered))]

  #ID KINOME Drugs (for this kinase)
  KINASES_kinome<-names(which(kinome_inhibitor_matrix_pKa_OVERLAP_filtered[drug,]>6))
  
  #PERCENT OF DrugBank drugs w/n kinome:

  tmp_num<-(length(KINASES_kinome)-length(which(KINASES_drugbank %in% KINASES_kinome)))
  drug_unique_num<-append(drug_unique_num,tmp_num)
  
  tmp_overlap<-length(which(KINASES_drugbank %in% KINASES_kinome))
  drug_overlap_num<-append(drug_overlap_num,tmp_overlap)
  
  
  #################################
  #2 POPULATE RANK-MATRIX OF PREDICTIONS
  #################################
  
  Bootstrapped_targets[drug,"Atom_overlap"]<-mean(final_results_processed$Atom[drug,KINASES_drugbank])/mean(final_results_processed$Atom[drug,])
  Bootstrapped_targets[drug,"Atom_kinome"]<-mean(final_results_processed$Atom[drug,KINASES_kinome])/mean(final_results_processed$Atom[drug,])
  
  
  Bootstrapped_targets[drug,"SBNB_overlap"]<-mean(final_results_processed$SBNB[drug,KINASES_drugbank])/mean(final_results_processed$SBNB[drug,])
  Bootstrapped_targets[drug,"SBNB_kinome"]<-mean(final_results_processed$SBNB[drug,KINASES_kinome])/mean(final_results_processed$SBNB[drug,])
  
  
  Bootstrapped_targets[drug,"netphar_overlap"]<-mean(final_results_processed$netphar[drug,KINASES_drugbank])/mean(final_results_processed$netphar[drug,])
  Bootstrapped_targets[drug,"netphar_kinome"]<-mean(final_results_processed$netphar[drug,KINASES_kinome])/mean(final_results_processed$netphar[drug,])
  
}

##########################
#TARGET-DATA OVERLAP: bar-plot
##########################
names(drug_unique_num)<-rownames(final_results_processed$SBNB)
names(drug_overlap_num)<-rownames(final_results_processed$SBNB)


DrugBank_Kinome_compare_drugs<-cbind(drug_overlap_num,drug_unique_num)
names_order<-names(sort(rowSums(DrugBank_Kinome_compare_drugs),decreasing=T))

names_order2<-names(sort(drug_overlap_num,decreasing=T))


par(las=2) # make label text perpendicular to axis
par(mar=c(9,5,4,5),cex=1.2) # increase y-axis margin.
#PLOT:
barplot(t(DrugBank_Kinome_compare_drugs[names_order,]),horiz=F,cex.lab=2,
        ylab="# of Kinase-Targets",xlab="",main="Gold-Standard Comparisons:",
        cex.axis = 1.5,col=c("black","red"),cex.main=2)
legend("topright",c("DrugBank-Kinome Overlap","Unique Kinome-Targets"),
       fill=c("black","red"),text.col=c("black","red"),cex=1.7)

save(list=c("DrugBank_Kinome_compare_drugs","Bootstrapped_targets"),file="KINOME-DrugBank_library-compare.RData")

##########################
#PREDICTIONS COMPARISONS:
##########################

plot(Bootstrapped_targets[,"Atom_overlap"],Bootstrapped_targets[,"Atom_kinome"], col="white"
     ,xlab="DrugBank-Targets Average-Score (average)",ylab="Unique Kinome-Targets Avg-Score")
text(Bootstrapped_targets[,"Atom_overlap"],Bootstrapped_targets[,"Atom_kinome"], labels=rownames(Bootstrapped_targets))


plot(Bootstrapped_targets[,"netphar_overlap"],Bootstrapped_targets[,"netphar_kinome"], col="white"
     ,xlab="DrugBank-Targets Average-Score (average)",ylab="Unique Kinome-Targets Avg-Score")
text(Bootstrapped_targets[,"netphar_overlap"],Bootstrapped_targets[,"netphar_kinome"], labels=rownames(Bootstrapped_targets))

plot(Bootstrapped_targets[,"SBNB_overlap"],Bootstrapped_targets[,"SBNB_kinome"], col="white"
     ,xlab="DrugBank-Targets Average-Score (average)",ylab="Unique Kinome-Targets Avg-Score")
text(Bootstrapped_targets[,"SBNB_overlap"],Bootstrapped_targets[,"SBNB_kinome"], labels=rownames(Bootstrapped_targets))

#CONCENSUS:
plot(rowMeans(apply(Bootstrapped_targets[,grep("overlap",colnames(Bootstrapped_targets))],2,rank)),
     rowMeans(apply(Bootstrapped_targets[,grep("kinome",colnames(Bootstrapped_targets))],2,rank)), 
     col="white",xlab="DrugBank-Targets Average-Score (Rank)",ylab="Unique Kinome-Targets Avg-Score (Rank)")
text(rowMeans(apply(Bootstrapped_targets[,grep("overlap",colnames(Bootstrapped_targets))],2,rank)),
     rowMeans(apply(Bootstrapped_targets[,grep("kinome",colnames(Bootstrapped_targets))],2,rank)),
     labels=rownames(Bootstrapped_targets))











####################################################################################################################################
#  Kinases-Drug <uM events NOT in DrugBank:
#####################################################################################################################################

kinome_unique_fract<-c()
kinome_unique_num<-c()
overlap_num<-c()


#target="EGFR"
for (target in colnames(kinome_inhibitor_matrix_pKa_OVERLAP)){
  #ID DrugBank Drugs (for this kinase):
  targets_list<-strsplit(n1druglibrary_annotation_CTEP_OVERLAP$target,",")
  target_drug_idx<-which(unlist(lapply(targets_list,function(x) (target %in% x)))==TRUE)
  DRUGS_drugbank<-rownames(n1druglibrary_annotation_CTEP_OVERLAP[target_drug_idx,])
  
  #ID KINOME Drugs (for this kinase)
  DRUGS_kinome<-names(which(kinome_inhibitor_matrix_pKa_OVERLAP[,target]>6))
  
  #PERCENT OF DrugBank drugs w/n kinome:
  tmp_fract<-(length(DRUGS_kinome)-length(which(DRUGS_drugbank %in% DRUGS_kinome)))/length(DRUGS_kinome)
  kinome_unique_fract<-append(kinome_unique_fract,tmp_fract)
  
  tmp_num<-(length(DRUGS_kinome)-length(which(DRUGS_drugbank %in% DRUGS_kinome)))
  kinome_unique_num<-append(kinome_unique_num,tmp_num)
  
  tmp_overlap<-length(which(DRUGS_drugbank %in% DRUGS_kinome))
  overlap_num<-append(overlap_num,tmp_overlap)
  
}
names(kinome_unique_num)<-colnames(kinome_inhibitor_matrix_pKa_OVERLAP)
names(kinome_unique_fract)<-colnames(kinome_inhibitor_matrix_pKa_OVERLAP)
names(overlap_num)<-colnames(kinome_inhibitor_matrix_pKa_OVERLAP)


DrugBank_Kinome_compare<-cbind(overlap_num,kinome_unique_num)
names_order<-names(sort(rowSums(DrugBank_Kinome_compare),decreasing=T))



par(las=2) # make label text perpendicular to axis
par(mar=c(5,5,4,5),cex=1) # increase y-axis margin.
#PLOT:
barplot(t(DrugBank_Kinome_compare[names_order[1:40],]),horiz=F,cex.main=1.5,cex.lab=1.5,
        ylab="Correlation",xlab="135 drugs",main="PANACEA-L1000 Correlation:  Drug-Signatures",
        cex.axis = 1.5,col=c("black","red"))


legend("topright",legend_info,cex=1,fill="red",text.col="darkred")












