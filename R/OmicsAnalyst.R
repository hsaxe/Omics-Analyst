#' Takes an OmicsAnalyst formatted data.frame, replaces any missing value(s) with zeros, calculates the interquartile ranges (IQRs) for all variables, sorts them by IQR, and filters the data by either IQR.Rank or IQR.
#' @param x a data.frame
#' @param Filter If "IQR.Rank", data will be IQR filtered by the specified Rank.threshold (default 5000). If IQR, data will be IQR filtered by the specified iqr.threshold (defaut 0.5).
#' @param Rank.Threshold a value indicating the maximum amount of IQR ranked variables to be included in the output.
#' @param iqr.threshold a value indicating the minimum IQR variables to be included in the output.
#' @import DescTools
#' @export
clean_data <- function(x, CPMfilter = T, CPMfilterFUN = mean, CPMfilter.Threshold = 20, IQRfilter = "", Rank.Threshold = 5000, iqr.threshold = 0.5){
  iqrs <- NULL
  a <- as.data.frame(x)

  ifelse(any(is.na(a)), print(paste(sum(is.na(a)), "missing value(s) replaced with zeros")), print("Data has no missing values"))

  if(any(is.na(a))){a[is.na(a)] <- 0}

  CPMfilterFUN = enquo(CPMfilterFUN)

  if(CPMfilter == T){
    final = a %>%
      edgeR::cpm() %>%
      as.data.frame() %>%
      mutate(!!CPMfilterFUN := apply(., 1, !!CPMfilterFUN)) %>%
      arrange(desc(!!CPMfilterFUN)) %>%
      filter(!!CPMfilterFUN >= CPMfilter.Threshold)

    if(nrow(a) <= Rank.Threshold) {print(paste("No features removed. Data only has", nrow(x), "features"))}

    print(paste("Removed", (nrow(x)) - (nrow(final)), "features based on CPMfilter.Threshold threshold"))

    return(final)

  }

  if(IQRfilter %in% "IQR.Rank"){

    final = a %>%
      mutate(iqrs = apply(., 1, IQR)) %>%
      arrange(desc(iqrs)) %>%
      head(n = Rank.Threshold)

    if(nrow(a) <= Rank.Threshold) {print(paste("No features removed. Data only has", nrow(x), "features"))}

    print(paste("Removed", (nrow(x)) - (nrow(final)), "features based on IQR rank threshold"))

    Cleaned = list(dat = select(final, !iqrs), IQRs = final %>% select(iqrs))

    return(Cleaned)

  }

  if(IQRfilter %in% "IQR"){

    final = a %>%
      mutate(iqrs = apply(., 1, IQR)) %>%
      filter(iqrs >= iqr.threshold)

    print(paste("Removed", (nrow(x)-1) - (nrow(final)-1), "features based on IQR threshold"))

    Cleaned = list(dat = select(final, !iqrs), IQRs = final %>% select(iqrs))

    return(Cleaned)
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

#' Takes an OmicsAnalyst formatted data.frame and returns a data.frame prepared for PCA analysis
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


