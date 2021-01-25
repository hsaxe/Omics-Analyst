#' Takes an OmicsAnalyst formatted data.frame, replaces any missing value(s) with zeros, calculates the interquartile ranges (IQRs) for all variables and sorts them by IQR. Only top varbiables(rows) in this rank will be retained in the output. The default is to retain the 5000 highest IQR variables to conserve computing resources
#' @param x a data.frame
#' @param Filter logical argument. If TRUE (T), data will be IQR filtered by the specified threshold. If FALSE (F), data will not be filtered
#' @param Threshold a value indicating the maximum amount of IQR ranked variables to be included in the output.
#' @import DescTools
#' @export
clean_data <- function(x, Filter = T, Threshold = 5000){
  iqrs <- NULL
  a <- x[-1,]

  ifelse(any(is.na(a)), print(paste(sum(is.na(a)), "missing value(s) replaced with zeros")), print("Data has no missing values"))

  if(any(is.na(a))){a[is.na(a)] <- 0}

  if(Filter %in% T){
    a$iqrs = as.numeric(matrixStats::rowIQRs(as.matrix(sapply(a[,-1], as.numeric))))

    filtered <- dplyr::arrange(a[-1,], desc(iqrs))[1:Threshold,]

    if(nrow(a) <= Threshold) {print(paste("No features removed. Data only has", nrow(x), "features"))}

    filtered[1,] <- x[1,]

    final <- filtered[,-ncol(filtered)]

    print(paste("Removed", (nrow(x)-1) - (nrow(final)-1), "features based on specified threshold"))

    return(final)
  }

}

#' Takes an OmicsAnalyst formatted data.frame and converts it to long format for plotting descriptive statistics and individual features
#' @param x a data.frame
#' @export

format_to_plot <- function(x){

  a <- reshape2::melt(x[-1,-ncol(x)], id.vars = 1, variable.name = "Sample", value.name = "value")

  a$Class = gsub("\\.\\d$", "", a$Sample)

  a$value = as.numeric(a$value)

  return(a)

}


#' Takes a plotting data.frame (typically long format) and returns the data.frame/data.table with a new column grouped by "group", transformed by "fun", and named by "name" arguments.
#' @param x a data.frame
#' @param group column name of groups to group by
#' @param variable column name of variable that will be grouped by "groups"
#' @param fun function to apply to grouped variable
#' @param name desired name for new grouped variable produced by "fun"
#' @export
group_stats <- function(x, group = "", variable = "", fun, name = ""){

  d <- data.table::data.table(x)

  a <- d[, fun(get(variable)), by = group]

  colnames(a)[ncol(a)] <- name

  m <- merge(d, a, by = group)

  return(m)
}

#' Takes an OmicsAnalyst formatted data.frame and performs general logarithm transformation (glog) and/or stdandardization (autoscale) returning a normalized data.frame
#' @param x a data.frame
#' @param glog logical argument. If TRUE (T), data will be glog normalized. If FALSE (F), data will not be glog normalized
#' @param autoscale logical argument. If TRUE (T), data will be autoscaled. If FALSE (F), data will not be autoscaled
#' @export
norm_Omics_df <- function(x, glog = T, autoscale = T){

  rawT <- data.table::transpose(x, keep.names = , make.names = 1)

  rawT[, c(2:ncol(rawT))] <- sapply(rawT[, c(2:ncol(rawT))], FUN = as.numeric)

  if(glog %in% T){
    rawT[, c(2:ncol(rawT))] <- sapply(rawT[, c(2:ncol(rawT))], FUN = FitAR::glog)
  }

  if(autoscale %in% T){
    rawT[, c(2:ncol(rawT))] <- sapply(rawT[, c(2:ncol(rawT))], FUN = BBmisc::normalize, method = "standardize")
  }

  rawT[, c(2:ncol(rawT))] <- sapply(rawT[, c(2:ncol(rawT))], FUN = format, digits = 3, width = 3)

  raw <- data.table::transpose(rawT, keep.names = "Gene.ID")

  colnames(raw) = colnames(x)

  return(raw)

}

#' Takes an OmicsAnalyst formatted data.frame and returns a data.fram prepared for PCA analysis
#' @param x a data.frame
#' @param legend.hjust adjusts horizontal justification of PCA plot legend
#' @import stats
#' @import ggplot2
#' @import grDevices
#' @import viridis
#' @export
PCA_Omics <- function(x, legend.hjust = 0.84){
  Group <- PC1 <- PC2 <- NULL
  forPCA <- data.table::transpose(x, keep.names = , make.names = 1)

  forPCA[2:ncol(forPCA)] <- sapply(forPCA[2:ncol(forPCA)], as.numeric)

  rownames(forPCA) <- colnames(x)[-1]

  PCA <- prcomp(forPCA[2:ncol(forPCA)], center = F, scale. = F)

  pca <- as.data.frame(PCA$x)

  proportion <- round((PCA$sdev^2)/ sum(PCA$sdev^2) * 100, 2)

  label <- paste(colnames(pca)[1:13], "(", as.character(proportion), "%)", sep = "")

  pca$Group = as.factor(gsub("\\.\\d$", "", rownames(pca)))

  print(summary(PCA))

  tiff("PCA.tiff", 8.5, 8.5, units = "cm", compression = "lzw", res = 800)
  g <- ggplot2::ggplot(pca, aes(x = PC1, y = PC2, fill = Group))+
    scale_fill_viridis(begin = 0.2, discrete = T)+
    scale_color_viridis(begin = 0.2, discrete = T)+
    geom_point(aes(fill = Group), size = 1, shape = 21, colour = "black")+
    stat_ellipse(geom = "polygon", alpha = 0.25, level = 0.95)+
    theme_bw()+
    theme(text = element_text(size = 7), legend.title = element_blank(), legend.text = element_text(size = 6), legend.key.size = unit(0.4,"cm"), legend.position = c(legend.hjust, 0.9), legend.background = element_rect(linetype = "solid", colour = "black", size = 0.3), legend.margin = margin(r=1, l=1,t=-2,b=1))+
    labs(x = label[1], y = label[2])
  print(g)
  dev.off()

  print(g)
}


